import 'dart:math' as math;

import 'package:vector_math/vector_math_64.dart';

import '../ijo_log.dart';
import '../measure/bounds.dart';
import '../mesh/mesh_kit.dart';
import 'boundary_boolean.dart';
import 'classify.dart';
import 'geom.dart';
import 'heal.dart';
import 'intersect.dart';
import 'provenance.dart';
import 'solid.dart';
import 'contact.dart';
import 'cap_join.dart';
import 'tessellate.dart';
import 'xform.dart';

/// Imprint / classify / collect boolean. Planar solids are general; curved
/// faces split when the intersection is a circle or axis-parallel generator.
class BrepBoolean {
  BrepBoolean._();

  static const unionOp = 0;
  static const subtractOp = 1;
  static const intersectOp = 2;

  /// Runs the boolean and carries face and edge provenance across it, so a
  /// face split by the cut still points at the face it came from and the cut's
  /// own new faces are marked as such. Set [sharedBoundary] to retain finite
  /// face arrangements directly when clipping along existing planar sides.
  static BrepSolid run(
    BrepSolid a,
    BrepSolid b,
    int op, {
    bool sharedBoundary = false,
  }) => BrepProvenance.inherit(
    BrepHeal.splitPlanarIslands(
      sharedBoundary ? _closedSharedBoundary(a, b, op) : _run(a, b, op),
    ),
    [a, b],
    op: 'boolean ${opName(op)}',
  );

  /// Preserve finite face regions when clipping along an existing face.
  /// The independent planar collector can mistake a concave cap's holes for
  /// outside regions when the cutter shares a long side of that cap.
  static BrepSolid _closedSharedBoundary(BrepSolid a, BrepSolid b, int op) {
    final result = _disconnectNonManifold(
      BrepBoundaryBoolean.run(a, b, op),
      splitPinches: false,
    );
    if (result.faces.isEmpty) throw StateError('Boolean produced no geometry');
    if (!BrepHeal.topologyIsClosed(result)) {
      throw BrepUnsupported(
        'Boolean could not construct a closed shared boundary',
      );
    }
    return result;
  }

  /// Sew a smaller cap into a convex planar face when the bodies occupy
  /// opposite half-spaces. This is an exact union, with no surface splitting:
  /// the removed cap becomes the shoulder's inner boundary. Useful for heads,
  /// bosses and threaded shafts, in either operand order.
  static BrepSolid? _joinContainedPlanarCap(BrepSolid a, BrepSolid b) {
    bool onNegativeSide(BrepSolid solid, Vector3 origin, Vector3 normal) {
      // In these surfaces the axial coordinate has no interior maximum;
      // supported boundary curves attain their extrema at their endpoints.
      bool parallel(Vector3 axis) => axis.cross(normal).length < 1e-8;
      for (final f in solid.faces) {
        final s = f.surface;
        if (s is BrepCylinder && !parallel(s.axis)) return false;
        if (s is BrepRuledSurface) {
          if (s.lower is! BrepHelix ||
              s.upper is! BrepHelix ||
              !parallel((s.lower as BrepHelix).axis) ||
              !parallel((s.upper as BrepHelix).axis))
            return false;
        } else if (s is! BrepPlane && s is! BrepCylinder)
          return false;
        for (final l in f.loops) {
          for (final co in l.coedges) {
            final c = solid.edges[co.edge].curve;
            if (c is BrepCircle && !parallel(c.normal)) return false;
            if (c is BrepHelix && !parallel(c.axis)) return false;
            if (c is BrepPlaneSectionCurve && !parallel(c.plane.normal))
              return false;
            if (c is! BrepLine &&
                c is! BrepCircle &&
                c is! BrepHelix &&
                c is! BrepPlaneSectionCurve)
              return false;
            if ((c.evaluate(c.t0) - origin).dot(normal) > 1e-7 ||
                (c.evaluate(c.t1) - origin).dot(normal) > 1e-7)
              return false;
          }
        }
      }
      return true;
    }

    for (var ai = 0; ai < a.faces.length; ai++) {
      final fa = a.faces[ai];
      if (fa.surface is! BrepPlane || fa.loops.length != 1) continue;
      final plane = fa.surface as BrepPlane, n = fa.outwardNormal(0, 0);
      for (var bi = 0; bi < b.faces.length; bi++) {
        final fb = b.faces[bi];
        if (fb.surface is! BrepPlane ||
            fb.loops.length != 1 ||
            fb.outwardNormal(0, 0).dot(n) > -1 + 1e-10)
          continue;
        final pb = fb.surface as BrepPlane;
        if ((pb.origin - plane.origin).dot(n).abs() > 1e-7) continue;
        final loop = fb.loops.single;
        if (loop.coedges.length < 3 ||
            loop.coedges.any((c) => b.edges[c.edge].curve is! BrepLine))
          continue;
        final points = [
          for (final co in loop.coedges)
            b.vertices[co.reversed ? b.edges[co.edge].v1 : b.edges[co.edge].v0],
        ];
        var center = Vector3.zero();
        for (final p in points) center += p;
        center /= points.length.toDouble();
        final halfSpaces = <(Vector3, Vector3)>[];
        var convex = true;
        for (var i = 0; i < points.length; i++) {
          final p = points[i], d = points[(i + 1) % points.length] - p;
          if (d.length < 1e-9) {
            convex = false;
            break;
          }
          var inward = n.cross(d).normalized();
          if ((center - p).dot(inward) < 0) inward = -inward;
          if (points.any((q) => (q - p).dot(inward) < -1e-7)) {
            convex = false;
            break;
          }
          halfSpaces.add((p, inward));
        }
        if (!convex) continue;
        final box = conservativeSolidBounds(
          BrepSolid(vertices: a.vertices, edges: a.edges, faces: [fa]),
        );
        if (box == null) continue;
        final corners = [
          for (final x in [box.min.x, box.max.x])
            for (final y in [box.min.y, box.max.y])
              for (final z in [box.min.z, box.max.z]) Vector3(x, y, z),
        ];
        if (halfSpaces.any(
          (h) => corners.any((p) => (p - h.$1).dot(h.$2) <= 1e-7),
        ))
          continue;
        if (!onNegativeSide(a, plane.origin, n) ||
            !onNegativeSide(b, plane.origin, -n))
          continue;
        final compound = BrepSolid.compound([a, b]);
        final shoulder = compound.faces[a.faces.length + bi];
        final faces = [
          for (var i = 0; i < compound.faces.length; i++)
            if (i != ai)
              i == a.faces.length + bi
                  ? BrepFace(
                      surface: shoulder.surface,
                      loops: [
                        ...shoulder.loops,
                        compound.faces[ai].loops.single,
                      ],
                      reversed: shoulder.reversed,
                      origin: shoulder.origin,
                    )
                  : compound.faces[i],
        ];
        final result = BrepSolid(
          vertices: compound.vertices,
          edges: compound.edges,
          faces: faces,
        );
        if (BrepHeal.topologyIsClosed(result)) return result;
      }
    }
    return null;
  }

  static BrepSolid _run(BrepSolid a, BrepSolid b, int op) {
    if (identical(a, b)) {
      if (op == subtractOp) throw StateError('Boolean produced no geometry');
      if (op == unionOp || op == intersectOp) return a.clone();
    }
    BrepSolid ruledInput(BrepSolid s) {
      if (!s.faces.any(
        (f) =>
            f.surface is BrepSplineSurface &&
            (f.surface as BrepSplineSurface).ruled != null,
      ))
        return s;
      return BrepSolid(
        vertices: s.vertices,
        edges: s.edges,
        faces: [
          for (final f in s.faces)
            BrepFace(
              surface: f.surface is BrepSplineSurface
                  ? (f.surface as BrepSplineSurface).ruled ?? f.surface
                  : f.surface,
              loops: f.loops,
              reversed: f.reversed,
              origin: f.origin,
            ),
        ],
      );
    }

    a = ruledInput(
      BrepHeal.alignCircularTrims(
        BrepHeal.restorePrismaticWalls(BrepXform.materializeSplineBody(a)),
      ),
    );
    b = ruledInput(
      BrepHeal.alignCircularTrims(
        BrepHeal.restorePrismaticWalls(BrepXform.materializeSplineBody(b)),
      ),
    );
    final boundA = conservativeSolidBounds(a),
        boundB = conservativeSolidBounds(b);
    if (boundA != null &&
        boundB != null &&
        !boundA.intersects(boundB, tolerance: 1e-7)) {
      if (op == subtractOp) return a.clone();
      if (op == intersectOp) throw StateError('Boolean produced no geometry');
      final builder = BrepBuilder();
      for (final s in [a, b]) {
        final vo = builder.vertices.length, eo = builder.edges.length;
        for (final v in s.vertices) builder.vertex(v);
        for (final e in s.edges)
          builder.edge(e.v0 + vo, e.v1 + vo, e.curve, origin: e.origin);
        for (final f in s.faces)
          builder.face(
            f.surface,
            [
              for (final l in f.loops)
                BrepLoop([
                  for (final c in l.coedges)
                    BrepCoedge(c.edge + eo, reversed: c.reversed),
                ]),
            ],
            reversed: f.reversed,
            origin: f.origin,
          );
      }
      return builder.build();
    }
    final ta = BrepXform.unwrapAffine(a), tb = BrepXform.unwrapAffine(b);
    if (ta != null || tb != null) {
      if (op == unionOp) {
        final joined = joinMatchingPlanarCaps(a, b);
        if (joined != null) return joined;
      }
      // Merge coplanar patches in the native frame before wrapping them.
      // The healer cannot recognize affine-wrapped planes; leaving the
      // subdivisions here exposes old operand outlines in X-ray and selection.
      if (ta != null && tb != null && ta.transform.sameAs(tb.transform)) {
        return BrepXform.affine(
          BrepHeal.run(_run(ta.solid, tb.solid, op)),
          ta.transform,
        );
      }
      // A similarity preserves analytic cylinders/cones/etc. Move the analytic
      // operand into the other body's native frame before trimming. An ordinary
      // (unwrapped) operand is already in world coordinates; it needs the same
      // treatment as an affine operand with an identity transform.
      BrepSolid? analyticInFrame(
        BrepSolid source,
        BrepAffineTransform relative,
      ) {
        final lengths = [
              relative.x.length,
              relative.y.length,
              relative.z.length,
            ],
            scale = lengths.first;
        final orthogonal =
            relative.x.dot(relative.y).abs() < 1e-10 * scale * scale &&
            relative.x.dot(relative.z).abs() < 1e-10 * scale * scale &&
            relative.y.dot(relative.z).abs() < 1e-10 * scale * scale;
        bool rationalRail(BrepCurve c) =>
            c is BrepLine ||
            c is BrepCircle ||
            c is BrepEllipse ||
            c is BrepSplineCurve ||
            c is BrepAffineCurve && rationalRail(c.base);
        if (orthogonal &&
            lengths.every((v) => (v - 1).abs() < 1e-10) &&
            relative.determinant > 0 &&
            source.faces.every(
              (f) =>
                  f.surface is BrepPlane ||
                  f.surface is BrepCylinder ||
                  f.surface is BrepCone ||
                  f.surface is BrepSphere ||
                  f.surface is BrepTorus ||
                  f.surface is BrepSplineSurface ||
                  f.surface is BrepRuledSurface &&
                      rationalRail((f.surface as BrepRuledSurface).lower) &&
                      rationalRail((f.surface as BrepRuledSurface).upper),
            )) {
          // Rotation/translation preserves rational and ruled surfaces exactly.
          // Materialize the placement so differing body frames do not recurse
          // through generic affine wrappers or reject ordinary sketch sweeps.
          return BrepXform.rigid(source, relative);
        }
        if (!orthogonal ||
            !lengths.every((v) => (v - scale).abs() < 1e-10 * scale) ||
            !source.faces.every(
              (f) =>
                  f.surface is BrepPlane ||
                  f.surface is BrepCylinder ||
                  f.surface is BrepCone ||
                  f.surface is BrepSphere ||
                  f.surface is BrepTorus,
            )) {
          return null;
        }
        return BrepXform.map(
          source,
          relative.point,
          reverseWinding: relative.determinant < 0,
          flipFaceNormals: false,
        );
      }

      if (ta != null) {
        final inverse = BrepAffineTransform.fromFunction(
          ta.transform.inversePoint,
        );
        final localB = analyticInFrame(
          tb?.solid ?? b,
          tb?.transform.then(inverse) ?? inverse,
        );
        if (localB != null) {
          return BrepXform.affine(
            BrepHeal.run(_run(ta.solid, localB, op)),
            ta.transform,
          );
        }
      }
      if (tb != null) {
        final inverse = BrepAffineTransform.fromFunction(
          tb.transform.inversePoint,
        );
        final localA = analyticInFrame(
          ta?.solid ?? a,
          ta?.transform.then(inverse) ?? inverse,
        );
        if (localA != null) {
          return BrepXform.affine(
            BrepHeal.run(_run(localA, tb.solid, op)),
            tb.transform,
          );
        }
      }
      if (ta != null &&
          b.faces.every((f) => f.surface is BrepPlane) &&
          b.edges.every((e) => e.curve is BrepLine)) {
        final localB = BrepXform.map(
          b,
          ta.transform.inversePoint,
          reverseWinding: ta.transform.determinant < 0,
          flipFaceNormals: false,
        );
        return BrepXform.affine(
          BrepHeal.run(_run(ta.solid, localB, op)),
          ta.transform,
        );
      }
      if (tb != null &&
          a.faces.every((f) => f.surface is BrepPlane) &&
          a.edges.every((e) => e.curve is BrepLine)) {
        final localA = BrepXform.map(
          a,
          tb.transform.inversePoint,
          reverseWinding: tb.transform.determinant < 0,
          flipFaceNormals: false,
        );
        return BrepXform.affine(
          BrepHeal.run(_run(localA, tb.solid, op)),
          tb.transform,
        );
      }
      // Differently scaled translational walls can still be represented
      // exactly as ruled curves. Keep the cheap common-frame paths above,
      // then expose supported world geometry for the general BREP trimmer.
      final worldA = BrepXform.materializeForEditing(a);
      final worldB = BrepXform.materializeForEditing(b);
      if (worldA.faces.every((f) => f.surface is! BrepAffineSurface) &&
          worldB.faces.every((f) => f.surface is! BrepAffineSurface)) {
        return _run(worldA, worldB, op);
      }
      throw BrepUnsupported(
        'Boolean between different affine curved surfaces needs a general surface trimmer',
      );
    }
    if (op == unionOp) {
      final joined =
          _joinContainedPlanarCap(a, b) ?? _joinContainedPlanarCap(b, a);
      if (joined != null) return joined;
    }
    final sw = Stopwatch()..start();
    IjoLog.line(
      'brep boolean ${opName(op)} '
      'a=${a.faces.length}f b=${b.faces.length}f',
    );
    try {
      if (op == subtractOp) {
        final boundsA = conservativeSolidBounds(a);
        final boundsB = conservativeSolidBounds(b);
        if (boundsA != null &&
            boundsB != null &&
            !boundsA.intersects(boundsB, tolerance: 1e-6)) {
          // A proven no-op must not rebuild trims or split existing blends.
          // Clone because provenance is assigned to the result in place.
          IjoLog.line('brep boolean subtract disjoint: preserving input');
          return a.clone();
        }
      }
      if (hasPairedPlanarContacts(a) || hasPairedPlanarContacts(b)) {
        final result = _disconnectNonManifold(
          BrepBoundaryBoolean.run(a, b, op),
          splitPinches: false,
        );
        if (!BrepHeal.topologyIsClosed(result)) {
          throw BrepUnsupported(
            'Boolean could not preserve closed contact boundaries',
          );
        }
        return result;
      }
      if (_needsTorusBoundary(a, b) ||
          _needsSphereBoundary(a, b) ||
          [...a.faces, ...b.faces].any(
            (f) =>
                f.surface is BrepRuledSurface || f.surface is BrepSplineSurface,
          )) {
        final result = BrepHeal.sewClosedRims(
          BrepBoundaryBoolean.run(a, b, op),
        );
        if (!BrepHeal.topologyIsClosed(result)) {
          if (IjoLog.verbose) {
            final uses = <int, List<int>>{};
            for (var fi = 0; fi < result.faces.length; fi++) {
              for (final loop in result.faces[fi].loops)
                for (final co in loop.coedges) (uses[co.edge] ??= []).add(fi);
            }
            for (final entry in uses.entries.where(
              (e) => e.value.length != 2,
            )) {
              final e = result.edges[entry.key];
              IjoLog.line(
                'Unpaired BREP edge ${entry.key} faces=${entry.value} ${e.curve.runtimeType} ${e.curve.evaluate(e.curve.t0)} -> ${e.curve.evaluate(e.curve.t1)}',
              );
            }
          }
          throw BrepUnsupported(
            'Curved boolean could not construct a closed boundary',
          );
        }
        return result;
      }
      if (op == subtractOp && _needsSharedBoundary(a, b)) {
        try {
          final candidate = BrepBoundaryBoolean.subtract(a, b);
          if (BrepHeal.topologyIsClosed(candidate) &&
              (MeshKit.isClosedSolid(BrepTessellate.run(candidate)) ||
                  _axialAndPlanarOnly(a, b) &&
                      MeshKit.isClosedSolid(
                        BrepHeal.tessellateClosed(candidate, op: 'Boolean'),
                      ))) {
            IjoLog.line(
              'brep shared-boundary subtract done ${sw.elapsedMilliseconds}ms out=${candidate.faces.length}f',
            );
            return candidate;
          }
          IjoLog.line(
            'brep shared-boundary validation failed; trying existing path',
          );
        } on BrepUnsupported catch (e) {
          IjoLog.line('brep shared-boundary unsupported: $e');
        }
      }
      final pieces = <_Piece>[];
      _collect(a, b, op, fromA: true, out: pieces);
      _collect(b, a, op, fromA: false, out: pieces);
      if (pieces.isEmpty) {
        throw StateError('Boolean produced no geometry');
      }
      final assembled = _pairOpenTwins(_assemble(pieces));
      final result = BrepHeal.sewClosedRims(
        _pairOpenTwins(
          _disconnectNonManifold(
            BrepHeal.run(_disconnectNonManifold(assembled)),
          ),
        ),
      );
      if (result.faces.isEmpty) {
        throw StateError('Boolean produced no geometry');
      }
      if (!BrepHeal.topologyIsClosed(result) &&
          [...a.faces, ...b.faces].every(
            (f) =>
                f.surface is BrepPlane ||
                f.surface is BrepCylinder ||
                f.surface is BrepCone ||
                f.surface is BrepSphere,
          )) {
        // A sphere can cross several edges of a finite planar face. The
        // independent face collector may then leave circular sphere trims
        // paired with straight planar chords. Rebuild both sides from their
        // shared intersections, just as for intersecting axial surfaces.
        final shared = _disconnectNonManifold(
          BrepBoundaryBoolean.run(a, b, op),
          splitPinches: false,
        );
        if (!BrepHeal.topologyIsClosed(shared)) {
          if (IjoLog.verbose) {
            final uses = <int, List<int>>{};
            for (var i = 0; i < shared.faces.length; i++)
              for (final l in shared.faces[i].loops)
                for (final c in l.coedges) (uses[c.edge] ??= []).add(i);
            for (final e in uses.entries.where((e) => e.value.length != 2))
              IjoLog.line(
                'Unpaired shared edge faces=${e.value} ${shared.edges[e.key].curve.toJson()}',
              );
          }
          throw BrepUnsupported(
            'Boolean could not construct a closed shared boundary',
          );
        }
        IjoLog.line(
          'brep shared-boundary ${opName(op)} done '
          '${sw.elapsedMilliseconds}ms out=${shared.faces.length}f',
        );
        return shared;
      }
      IjoLog.line(
        'brep boolean done ${sw.elapsedMilliseconds}ms '
        'out=${result.faces.length}f',
      );
      return result;
    } catch (e, st) {
      IjoLog.fail('brep boolean', e, st);
      rethrow;
    }
  }

  static bool _needsSphereBoundary(BrepSolid a, BrepSolid b) {
    if ([...a.edges, ...b.edges].any(
      (e) =>
          e.curve is BrepIntersectionCurve &&
          (e.curve as BrepIntersectionCurve).target is BrepSphere,
    ))
      return true;
    for (final pair in [(a, b), (b, a)]) {
      for (final fa in pair.$1.faces) {
        final c = fa.surface;
        if (c is BrepPlane) {
          final polygon = _planePoly(pair.$1, fa);
          for (final fb in pair.$2.faces) {
            if (fb.surface is! BrepSphere) continue;
            if (!_sphereMayMeetPlanarFace(pair.$1, fa, pair.$2, fb)) continue;
            for (final curve in BrepIntersect.surfaces(c, fb.surface)) {
              if (curve is BrepCircle &&
                  !_circleIsInteriorHole(polygon, curve, c.normal) &&
                  _circleHitsPoly(curve, polygon, c.normal)) {
                // Clipped spherical arcs can enclose the large complement.
                // Build shared, oriented patches before choosing either side;
                // edge counts alone cannot detect a wrongly filled sphere.
                return true;
              }
            }
          }
        }
        if (c is! BrepCylinder) continue;
        for (final fb in pair.$2.faces) {
          final s = fb.surface;
          if (s is! BrepSphere) continue;
          final d = s.center - c.origin;
          if ((d - c.axis * d.dot(c.axis)).length > 1e-6 &&
              BrepBoundaryBoolean.facesIntersect(pair.$1, fa, pair.$2, fb))
            return true;
        }
      }
    }
    return false;
  }

  static bool _needsSharedBoundary(BrepSolid a, BrepSolid b) {
    if ([...a.edges, ...b.edges].any((e) => e.curve is BrepIntersectionCurve))
      return true;
    // A planar face with several trim loops needs their shared arrangement
    // when a cut reaches it too. Selecting this path only for cylinder
    // intersections lets a later planar notch rebuild a multiply-holed cap
    // through the legacy polygon path and lose parts of its trim boundary.
    for (final face in _axialAndPlanarOnly(a, b) ? a.faces : <BrepFace>[]) {
      if (face.surface is BrepPlane &&
          face.loops.length > 1 &&
          b.faces.any(
            (other) => BrepBoundaryBoolean.facesIntersect(a, face, b, other),
          )) {
        return true;
      }
    }
    for (final fa in a.faces) {
      if (fa.surface is! BrepCylinder) continue;
      for (final fb in b.faces) {
        if (fb.surface is! BrepCylinder && fb.surface is! BrepPlane) continue;
        if (BrepBoundaryBoolean.facesIntersect(a, fa, b, fb)) return true;
      }
    }
    return false;
  }

  static bool _axialAndPlanarOnly(BrepSolid a, BrepSolid b) =>
      [...a.faces, ...b.faces].every(
        (f) =>
            f.surface is BrepPlane ||
            f.surface is BrepCylinder ||
            f.surface is BrepCone,
      );

  static String opName(int op) {
    switch (op) {
      case subtractOp:
        return 'subtract';
      case intersectOp:
        return 'intersect';
      default:
        return 'union';
    }
  }

  static void _collect(
    BrepSolid self,
    BrepSolid other,
    int op, {
    required bool fromA,
    required List<_Piece> out,
  }) {
    for (final face in self.faces) {
      final surface = face.surface;
      if (surface is BrepPlane) {
        final sourcePoly = _planePoly(self, face);
        // Intersect finite faces, not the sphere hidden inside a union.
        final sphereFaces = [
          for (final otherFace in other.faces)
            if (otherFace.surface is BrepSphere &&
                _sphereMayMeetPlanarFace(self, face, other, otherFace))
              otherFace,
        ];
        var polys = [sourcePoly];
        if (polys.first.length < 3) continue;
        // Prism cut: imprint the other solid's profile on this cap
        // instead of slicing it with infinite wall planes. Subtract
        // imprints the tool onto the base; union imprints a simple
        // prism onto already-complex caps (holes / arcs) so a boss
        // sews instead of shattering the D-notch face. Through-cuts
        // and slice half-spaces outspan with line walls and fall
        // through to plane splits.
        final priorHoles = _faceHoles(self, face);
        if ((op == subtractOp && fromA) ||
            ((op == unionOp || op == intersectOp) &&
                _faceNeedsPrismImprint(self, face))) {
          final imprint = _imprintPrismOnPlane(
            self,
            face,
            other,
            op,
            priorHoles,
          );
          if (imprint != null) {
            if (!imprint.drop) {
              final n = face.outwardNormal(0, 0);
              final nn = surface.normal;
              for (final leftover in imprint.outers) {
                var poly = leftover.pts;
                final extraHoles = <List<Vector3>>[];
                final extraCurves = <List<BrepCurve>>[];
                final extraCircles = <BrepCircle>[];
                final notchArcs = <BrepCircle>[];
                for (var hi = 0; hi < priorHoles.holes.length; hi++) {
                  if (_holeSamplesInPoly(poly, priorHoles.holes[hi], nn)) {
                    extraHoles.add(priorHoles.holes[hi]);
                    extraCurves.add(
                      hi < priorHoles.curves.length
                          ? priorHoles.curves[hi]
                          : const <BrepCurve>[],
                    );
                  }
                }
                for (final hole in priorHoles.circles) {
                  if (_circleIsInteriorHole(poly, hole, nn)) {
                    extraCircles.add(hole);
                    continue;
                  }
                  final notched = _notchPolyWithCircle(poly, hole, nn);
                  if (notched == null) continue;
                  poly = notched.poly;
                  notchArcs.add(notched.arc);
                }
                final hasHole =
                    leftover.holes.isNotEmpty ||
                    leftover.holeCircles.isNotEmpty ||
                    extraHoles.isNotEmpty ||
                    extraCircles.isNotEmpty;
                if (!_keepNeighborhood(
                  op,
                  fromA: fromA,
                  other: other,
                  centroid: _probeOnOuter(poly),
                  outward: n,
                  hasHole: hasHole,
                )) {
                  continue;
                }
                out.add(
                  _Piece(
                    surface,
                    poly,
                    face.reversed,
                    holes: [...extraHoles, ...leftover.holes],
                    holeCircles: [...extraCircles, ...leftover.holeCircles],
                    holeCurves: [...extraCurves, ...leftover.holeCurves],
                    outerArcs: [...notchArcs, ...leftover.outerArcs],
                    outerCurves: leftover.curves,
                  ),
                );
              }
            }
            continue;
          }
        }
        for (final of in other.faces) {
          if (of.surface is! BrepPlane) continue;
          final pl = of.surface as BrepPlane;
          final otherPoly = _planePoly(other, of);
          final next = <List<Vector3>>[];
          for (final poly in polys) {
            if (!_aabbsOverlap(poly, otherPoly) ||
                !_finitePlanesHit(poly, surface, otherPoly, pl)) {
              next.add(poly);
              continue;
            }
            next.addAll(_splitPolyByPlane(poly, pl));
          }
          polys = next;
        }
        for (final of in other.faces) {
          if (of.surface is! BrepCylinder) continue;
          final cyl = of.surface as BrepCylinder;
          if (surface.normal.dot(cyl.axis).abs() > 1e-3) continue;
          final next = <List<Vector3>>[];
          for (final poly in polys) {
            next.addAll(_splitPolyByCylinderGenerators(poly, cyl, surface));
          }
          polys = next;
        }
        final holes = <List<List<Vector3>>>[];
        final holeCircles = <List<BrepCircle>>[];
        final arcs = <List<BrepCircle?>>[];
        for (final of in other.faces) {
          if (of.surface is! BrepCylinder) continue;
          final cyl = of.surface as BrepCylinder;
          final an = surface.normal.dot(cyl.axis).abs();
          if (an < 0.999) continue;
          final domain = BrepTessellate.cylinderFaceDomain(other, of);
          final zCap = (surface.origin - cyl.origin).dot(cyl.axis);
          if (zCap < domain.$3 - 1e-4 || zCap > domain.$4 + 1e-4) continue;
          final next = <List<Vector3>>[];
          final nextHoles = <List<List<Vector3>>>[];
          final nextCircles = <List<BrepCircle>>[];
          final nextArcs = <List<BrepCircle?>>[];
          for (var pi = 0; pi < polys.length; pi++) {
            final poly = polys[pi];
            final cut = _cutPlaneByCylinder(
              poly,
              cyl,
              surface,
              op: op,
              fromA: fromA,
              u0: domain.$1,
              u1: domain.$2,
              other: other,
            );
            for (final piece in cut.outers) {
              if (piece.poly.length < 2) continue;
              if (piece.poly.length < 3 && piece.arc12 == null) continue;
              next.add(piece.poly);
              nextHoles.add([
                ...((pi < holes.length) ? holes[pi] : const <List<Vector3>>[]),
                if (cut.hole != null) cut.hole!,
              ]);
              nextCircles.add([
                ...((pi < holeCircles.length)
                    ? holeCircles[pi]
                    : const <BrepCircle>[]),
                if (cut.holeCircle != null) cut.holeCircle!,
              ]);
              nextArcs.add([piece.arc12]);
            }
          }
          polys = next;
          holes
            ..clear()
            ..addAll(nextHoles);
          holeCircles
            ..clear()
            ..addAll(nextCircles);
          arcs
            ..clear()
            ..addAll(nextArcs);
          // Hanging nodes: a partial-arc split inserts T-junctions on one
          // piece only. Subdivide coincident edges on sibling faces so the
          // diameter sews (lens stays attached to the rest of the cap).
          final imprint = <Vector3>[];
          for (final list in arcs) {
            for (final arc in list) {
              if (arc == null) continue;
              imprint.add(arc.evaluate(arc.t0));
              imprint.add(arc.evaluate(arc.t1));
            }
          }
          if (imprint.isNotEmpty) {
            final nn = surface.normal;
            for (var i = 0; i < polys.length; i++) {
              var p = polys[i];
              for (final q in imprint) {
                p = _insertOnBoundary(p, q, nn) ?? p;
              }
              polys[i] = p;
            }
          }
        }
        for (final of in other.faces) {
          if (of.surface is! BrepCone) continue;
          final cone = of.surface as BrepCone;
          final an = surface.normal.dot(cone.axis).abs();
          if (an < 0.999) continue;
          final zCap = (surface.origin - cone.origin).dot(cone.axis);
          if (zCap < -1e-4 || zCap > cone.height + 1e-4) continue;
          final r = cone.radiusAt(zCap);
          if (r < 1e-6) continue;
          final dummy = BrepCylinder(
            origin: cone.origin,
            axis: cone.axis,
            radius: r,
            xDir: cone.xDir,
          );
          final next = <List<Vector3>>[];
          final nextHoles = <List<List<Vector3>>>[];
          final nextCircles = <List<BrepCircle>>[];
          final nextArcs = <List<BrepCircle?>>[];
          for (var pi = 0; pi < polys.length; pi++) {
            final poly = polys[pi];
            final cut = _cutPlaneByCylinder(
              poly,
              dummy,
              surface,
              op: op,
              fromA: fromA,
              other: other,
            );
            for (final piece in cut.outers) {
              if (piece.poly.length < 2) continue;
              if (piece.poly.length < 3 && piece.arc12 == null) continue;
              next.add(piece.poly);
              nextHoles.add([
                ...((pi < holes.length) ? holes[pi] : const <List<Vector3>>[]),
                if (cut.hole != null) cut.hole!,
              ]);
              nextCircles.add([
                ...((pi < holeCircles.length)
                    ? holeCircles[pi]
                    : const <BrepCircle>[]),
                if (cut.holeCircle != null) cut.holeCircle!,
              ]);
              nextArcs.add([piece.arc12]);
            }
          }
          if (next.isNotEmpty) {
            polys = next;
            holes
              ..clear()
              ..addAll(nextHoles);
            holeCircles
              ..clear()
              ..addAll(nextCircles);
            arcs
              ..clear()
              ..addAll(nextArcs);
          }
        }
        for (final of in sphereFaces) {
          final sph = of.surface as BrepSphere;
          for (final curve in BrepIntersect.surfaces(surface, sph)) {
            if (curve is! BrepCircle) continue;
            while (holeCircles.length < polys.length) {
              holeCircles.add(const []);
            }
            while (holes.length < polys.length) {
              holes.add(const []);
            }
            while (arcs.length < polys.length) {
              arcs.add(const [null]);
            }
            for (var pi = 0; pi < polys.length; pi++) {
              if (_circleIsInteriorHole(polys[pi], curve, surface.normal)) {
                holeCircles[pi] = [...holeCircles[pi], curve];
              }
            }
          }
        }
        for (final of in other.faces) {
          if (of.surface is! BrepTorus) continue;
          final tor = of.surface as BrepTorus;
          // A tangent plane touches the tube along a circle but removes no
          // planar area. Turning that contact into a hole opens the box floor.
          if (surface.normal.cross(tor.axis).length < 1e-8 &&
              (surface.height(tor.center).abs() - tor.minorRadius).abs() <
                  1e-8) {
            continue;
          }
          for (final curve in BrepIntersect.surfaces(surface, tor)) {
            if (curve is! BrepCircle) continue;
            while (holeCircles.length < polys.length) {
              holeCircles.add(const []);
            }
            while (holes.length < polys.length) {
              holes.add(const []);
            }
            while (arcs.length < polys.length) {
              arcs.add(const [null]);
            }
            for (var pi = 0; pi < polys.length; pi++) {
              if (_circleIsInteriorHole(polys[pi], curve, surface.normal)) {
                holeCircles[pi] = [...holeCircles[pi], curve];
              }
            }
          }
        }
        final n = face.outwardNormal(0, 0);
        for (var pi = 0; pi < polys.length; pi++) {
          var poly = polys[pi];
          final polyHoles = pi < holes.length
              ? holes[pi]
              : const <List<Vector3>>[];
          final circles = pi < holeCircles.length
              ? holeCircles[pi]
              : const <BrepCircle>[];
          final arc = pi < arcs.length && arcs[pi].isNotEmpty
              ? arcs[pi].first
              : null;
          if (poly.length < 2) continue;
          if (poly.length < 3 && arc == null) continue;
          final nn = surface.normal;
          final ellipses = [
            for (final of in other.faces)
              if (of.surface is! BrepPlane)
                for (final curve in BrepIntersect.surfaces(surface, of.surface))
                  if (curve is BrepEllipse &&
                      _ellipseIsInteriorHole(poly, curve, nn))
                    curve,
          ];
          final hasInner =
              polyHoles.isNotEmpty || circles.isNotEmpty || ellipses.isNotEmpty;
          final c = hasInner
              ? _probeOnOuter(poly)
              : arc != null
              ? (poly.length == 2 ? _arcProbe(poly, arc) : poly.first)
              : _centroid(poly);
          if (!_keepNeighborhood(
            op,
            fromA: fromA,
            other: other,
            centroid: c,
            outward: n,
            hasHole: hasInner,
          )) {
            continue;
          }
          var rev = face.reversed;
          if (!fromA && op == subtractOp) rev = !rev;
          // Keep analytic loops when the outer boundary was not changed.
          // A sampled arc is suitable for classification, not a replacement
          // for the shared curve on the adjacent fillet face. Sphere notches
          // are handled below; retain their existing interior-hole path only.
          final keepSource =
              polys.length == 1 &&
              identical(poly, sourcePoly) &&
              polyHoles.isEmpty &&
              ellipses.isEmpty &&
              arc == null;
          final extraHoles = keepSource
              ? const <List<Vector3>>[]
              : [
                  for (var hi = 0; hi < priorHoles.holes.length; hi++)
                    if (_holeSamplesInPoly(poly, priorHoles.holes[hi], nn))
                      priorHoles.holes[hi],
                ];
          final extraCurves = keepSource
              ? const <List<BrepCurve>>[]
              : [
                  for (var hi = 0; hi < priorHoles.holes.length; hi++)
                    if (_holeSamplesInPoly(poly, priorHoles.holes[hi], nn))
                      hi < priorHoles.curves.length
                          ? priorHoles.curves[hi]
                          : const <BrepCurve>[],
                ];
          final extraCircles = <BrepCircle>[];
          final notchArcs = <BrepCircle>[];
          if (!keepSource) {
            for (final hole in priorHoles.circles) {
              if (_circleIsInteriorHole(poly, hole, nn)) {
                extraCircles.add(hole);
                continue;
              }
              final notched = _notchPolyWithCircle(poly, hole, nn);
              if (notched == null) continue;
              poly = notched.poly;
              notchArcs.add(notched.arc);
            }
          }
          for (final of in sphereFaces) {
            for (final curve in BrepIntersect.surfaces(surface, of.surface)) {
              if (curve is! BrepCircle) continue;
              if (_circleIsInteriorHole(poly, curve, nn)) continue;
              final notched = _notchPolyWithCircle(poly, curve, nn);
              if (notched == null) continue;
              poly = notched.poly;
              notchArcs.add(notched.arc);
            }
          }
          out.add(
            _Piece(
              surface,
              poly,
              rev,
              holes: [...extraHoles, ...polyHoles],
              holeCircles: [...extraCircles, ...circles],
              holeCurves: extraCurves,
              arc12: arc,
              outerArcs: notchArcs,
              innerClosed: ellipses,
              sourceSolid: keepSource && notchArcs.isEmpty ? self : null,
              sourceFace: keepSource && notchArcs.isEmpty ? face : null,
            ),
          );
        }
      } else if (surface is BrepCylinder) {
        _collectCylinder(self, face, other, op, fromA: fromA, out: out);
      } else if (surface is BrepSphere) {
        _collectSphere(self, face, other, op, fromA: fromA, out: out);
      } else if (surface is BrepTorus) {
        _collectTorus(self, face, other, op, fromA: fromA, out: out);
      } else if (surface is BrepCone) {
        _collectCone(self, face, other, op, fromA: fromA, out: out);
      } else {
        throw BrepUnsupported(
          'Boolean on ${surface.runtimeType} is not implemented yet',
        );
      }
    }
  }

  /// Reject a planar face wholly outside a retained spherical patch. Circular
  /// trims define half-spaces in world coordinates, so the test still works
  /// after arbitrary rotations where the faces' axis-aligned bounds overlap.
  /// Bound analytic curves conservatively; sampled arcs can miss an extremum.
  static bool _sphereMayMeetPlanarFace(
    BrepSolid planarSolid,
    BrepFace planarFace,
    BrepSolid sphereSolid,
    BrepFace sphereFace,
  ) {
    final planarBounds = conservativeSolidBounds(
      BrepSolid(
        vertices: planarSolid.vertices,
        edges: planarSolid.edges,
        faces: [planarFace],
      ),
    );
    if (planarBounds != null &&
        !conservativeSphereFaceBounds(
          sphereSolid,
          sphereFace,
        ).intersects(planarBounds, tolerance: 1e-6)) {
      return false;
    }
    if (!sphereFace.loops.every(
      (loop) =>
          loop.coedges.length == 1 &&
          sphereSolid.edges[loop.coedges.single.edge].curve is BrepCircle &&
          (sphereSolid.edges[loop.coedges.single.edge].curve as BrepCircle)
              .isClosed,
    )) {
      return true;
    }
    for (final loop in sphereFace.loops) {
      final co = loop.coedges.single;
      final rim = sphereSolid.edges[co.edge].curve as BrepCircle;
      var n = co.reversed ? -rim.normal : rim.normal;
      // The retained side is <= 0, matching circular sphere classification.
      if (sphereFace.loops.length == 1) n = -n;
      var minimum = double.infinity;
      for (final boundary in planarFace.loops) {
        for (final use in boundary.coedges) {
          final edge = planarSolid.edges[use.edge];
          for (final p in [
            planarSolid.vertices[edge.v0],
            planarSolid.vertices[edge.v1],
          ]) {
            minimum = math.min(minimum, (p - rim.center).dot(n));
          }
          final curve = edge.curve;
          final double lower;
          if (curve is BrepLine) {
            lower = math.min(
              (curve.a - rim.center).dot(n),
              (curve.b - rim.center).dot(n),
            );
          } else if (curve is BrepCircle) {
            final x = curve.radius * curve.xDir.dot(n);
            final y = curve.radius * curve.yDir.dot(n);
            lower =
                (curve.center - rim.center).dot(n) - math.sqrt(x * x + y * y);
          } else if (curve is BrepEllipse) {
            final x = curve.rMajor * curve.xDir.normalized().dot(n);
            final y = curve.rMinor * curve.yDir.dot(n);
            lower =
                (curve.center - rim.center).dot(n) - math.sqrt(x * x + y * y);
          } else {
            return true;
          }
          minimum = math.min(minimum, lower);
        }
      }
      if (minimum.isFinite && minimum > 1e-6) return false;
    }
    return true;
  }

  static bool _sameCircle(BrepCircle a, BrepCircle b) {
    if ((a.center - b.center).length > 1e-4) return false;
    if ((a.radius - b.radius).abs() > 1e-4) return false;
    return a.normal.cross(b.normal).length < 2e-2;
  }

  static void _addUniqueCircle(List<BrepCircle> circles, BrepCircle curve) {
    for (final c in circles) {
      if (_sameCircle(c, curve)) return;
    }
    circles.add(curve);
  }

  static bool _sphereRemainderKept(
    BrepSphere sph,
    BrepFace face,
    BrepSolid other,
    int op, {
    required bool fromA,
  }) {
    const dirs = [
      [1.0, 0.0, 0.0],
      [-1.0, 0.0, 0.0],
      [0.0, 1.0, 0.0],
      [0.0, -1.0, 0.0],
      [0.0, 0.0, 1.0],
      [0.0, 0.0, -1.0],
      [0.7071, 0.7071, 0.0],
      [0.7071, -0.7071, 0.0],
      [-0.7071, 0.7071, 0.0],
      [-0.7071, -0.7071, 0.0],
      [0.7071, 0.0, 0.7071],
      [0.7071, 0.0, -0.7071],
    ];
    for (final d in dirs) {
      final p = sph.center + Vector3(d[0], d[1], d[2]) * sph.radius;
      var n = p - sph.center;
      if (n.length < 1e-12) continue;
      n.normalize();
      if (face.reversed) n = -n;
      if (_keepNeighborhood(
        op,
        fromA: fromA,
        other: other,
        centroid: p,
        outward: n,
        hasHole: true,
      )) {
        return true;
      }
    }
    return false;
  }

  static void _collectSphere(
    BrepSolid self,
    BrepFace face,
    BrepSolid other,
    int op, {
    required bool fromA,
    required List<_Piece> out,
  }) {
    final sph = face.surface as BrepSphere;
    final otherBounds = conservativeSolidBounds(other);
    if (otherBounds != null &&
        !conservativeSphereFaceBounds(
          self,
          face,
        ).intersects(otherBounds, tolerance: 1e-6)) {
      // The trimmed face is outside the other operand even when its full
      // supporting sphere intersects an infinite cutter surface.
      if (op == unionOp || (op == subtractOp && fromA)) {
        out.add(
          _Piece(
            sph,
            const [],
            face.reversed,
            sourceSolid: self,
            sourceFace: face,
          ),
        );
      }
      return;
    }

    final circles = <BrepCircle>[];
    final clipped = <BrepCircle>[];
    var crossesEdge = false;
    final existingHoles = <List<BrepCurve>>[];
    for (final loop in face.loops) {
      if (loop.coedges.length == 1) {
        final ei = loop.coedges.first.edge;
        if (ei >= 0 && ei < self.edges.length) {
          final curve = self.edges[ei].curve;
          if (curve is BrepCircle && curve.isClosed) {
            _addUniqueCircle(circles, curve);
            continue;
          }
        }
      }
      if (loop.coedges.isEmpty) continue;
      final curves = <BrepCurve>[];
      var ok = true;
      for (final c in loop.coedges) {
        if (c.edge < 0 || c.edge >= self.edges.length) {
          ok = false;
          break;
        }
        curves.add(self.edges[c.edge].curve);
      }
      if (ok && curves.isNotEmpty) existingHoles.add(curves);
    }
    final priorCircles = circles.length;
    for (final of in other.faces) {
      if (of.surface is BrepPlane) {
        final pl = of.surface as BrepPlane;
        final poly = _planePoly(other, of);
        if (poly.length < 3) continue;
        for (final curve in BrepIntersect.surfaces(pl, sph)) {
          if (curve is! BrepCircle) continue;
          if (_circleIsInteriorHole(poly, curve, pl.normal)) {
            _addUniqueCircle(circles, curve);
          } else if (_circleHitsPoly(curve, poly, pl.normal)) {
            crossesEdge = true;
            clipped.addAll(_clipCircleToPoly(curve, poly, pl.normal));
          }
        }
        continue;
      }
      if (of.surface is BrepSphere) {
        for (final curve in BrepIntersect.surfaces(sph, of.surface)) {
          if (curve is BrepCircle) _addUniqueCircle(circles, curve);
        }
        continue;
      }
      if (of.surface is BrepCylinder) {
        for (final curve in BrepIntersect.surfaces(sph, of.surface)) {
          if (curve is BrepCircle) _addUniqueCircle(circles, curve);
        }
      }
    }
    final noNew = clipped.isEmpty && circles.length == priorCircles;
    if (noNew) {
      if (crossesEdge) {
        throw BrepUnsupported(
          'Sphere boolean that splits the sphere across a face edge '
          'is not implemented yet',
        );
      }
      final probe = sph.center + Vector3(0, 0, sph.radius * 0.5);
      var n = probe - sph.center;
      if (n.length < 1e-12) n = Vector3(0, 0, 1);
      n.normalize();
      if (face.reversed) n = -n;
      if (!_keepNeighborhood(
        op,
        fromA: fromA,
        other: other,
        centroid: probe,
        outward: n,
      )) {
        return;
      }
      out.add(
        _Piece(
          sph,
          const [],
          face.reversed,
          sourceSolid: self,
          sourceFace: face,
        ),
      );
      return;
    }
    if (circles.isEmpty && clipped.isEmpty && existingHoles.isEmpty) {
      if (crossesEdge) {
        throw BrepUnsupported(
          'Sphere boolean that splits the sphere across a face edge '
          'is not implemented yet',
        );
      }
      final probe = sph.center + Vector3(0, 0, sph.radius * 0.5);
      var n = probe - sph.center;
      if (n.length < 1e-12) n = Vector3(0, 0, 1);
      n.normalize();
      if (face.reversed) n = -n;
      if (!_keepNeighborhood(
        op,
        fromA: fromA,
        other: other,
        centroid: probe,
        outward: n,
      )) {
        return;
      }
      out.add(
        _Piece(
          sph,
          const [],
          face.reversed,
          sourceSolid: self,
          sourceFace: face,
        ),
      );
      return;
    }
    if (clipped.isNotEmpty) {
      final loops = _chainCircleArcs(clipped);
      if (loops.isEmpty) {
        throw BrepUnsupported(
          'Sphere boolean that splits the sphere across a face edge '
          'is not implemented yet',
        );
      }
      final hasPrior = circles.isNotEmpty || existingHoles.isNotEmpty;
      if (hasPrior) {
        for (final loop in loops) {
          if (loop.isEmpty) continue;
          if (!_clipLoopOnSphereFace(sph, circles, loop)) continue;
          existingHoles.add([for (final s in loop) s.curve]);
        }
      } else {
        for (final loop in loops) {
          final pts = [for (final s in loop) s.a.clone()];
          if (pts.length < 2) continue;
          var probe = Vector3.zero();
          for (final p in pts) {
            probe += p;
          }
          probe = sph.closestPoint(probe / pts.length.toDouble());
          final small = probe;
          final large = sph.center * 2.0 - small;
          Vector3? kept;
          for (final pole in [small, large]) {
            var n = pole - sph.center;
            if (n.length < 1e-12) continue;
            n.normalize();
            if (face.reversed) n = -n;
            if (!_keepNeighborhood(
              op,
              fromA: fromA,
              other: other,
              centroid: pole,
              outward: n,
            )) {
              continue;
            }
            kept = pole;
            break;
          }
          if (kept == null) continue;
          var rev = face.reversed;
          if (!fromA && op == subtractOp) rev = !rev;
          out.add(
            _Piece(
              sph,
              [kept, ...pts],
              rev,
              outerArcs: [
                for (final s in loop)
                  if (s.curve is BrepCircle) s.curve as BrepCircle,
              ],
              outerCurves: [for (final s in loop) s.curve],
            ),
          );
        }
        if (circles.isEmpty) return;
      }
    }
    if (circles.length > 1 || existingHoles.isNotEmpty) {
      final remainderKept = _sphereRemainderKept(
        sph,
        face,
        other,
        op,
        fromA: fromA,
      );
      if (remainderKept) {
        var rev = face.reversed;
        if (!fromA && op == subtractOp) rev = !rev;
        out.add(
          _Piece(
            sph,
            const [],
            rev,
            holeCircles: [...circles],
            holeCurves: existingHoles,
          ),
        );
        return;
      }
    }
    for (final circle in circles) {
      final polesPair = _sphereCapPoles(sph, circle);
      final small = polesPair.$1;
      final large = polesPair.$2;
      final poles = circles.length > 1 ? [small] : [small, large];
      for (final pole in poles) {
        var n = pole - sph.center;
        if (n.length < 1e-12) continue;
        n.normalize();
        if (face.reversed) n = -n;
        final keep = _keepNeighborhood(
          op,
          fromA: fromA,
          other: other,
          centroid: pole,
          outward: n,
        );
        if (!keep) {
          continue;
        }
        var rev = face.reversed;
        if (!fromA && op == subtractOp) rev = !rev;
        out.add(_Piece(sph, [pole], rev, holeCircles: [circle]));
      }
    }
  }

  static bool _clipLoopOnSphereFace(
    BrepSphere sph,
    List<BrepCircle> circles,
    List<_ImprintSeg> loop,
  ) {
    if (circles.length < 2) return true;
    var p = Vector3.zero();
    for (final s in loop) {
      p += s.a;
    }
    if (loop.isEmpty) return false;
    p = sph.closestPoint(p / loop.length.toDouble());
    BrepCircle? a;
    BrepCircle? b;
    for (var i = 0; i < circles.length; i++) {
      for (var j = i + 1; j < circles.length; j++) {
        if (circles[i].normal.dot(circles[j].normal).abs() > 0.9) {
          a = circles[i];
          b = circles[j];
          break;
        }
      }
      if (a != null) break;
    }
    if (a == null || b == null) return true;
    final n = a.normal.normalized();
    final za = a.center.dot(n);
    final zb = b.center.dot(n);
    final zp = p.dot(n);
    final lo = math.min(za, zb) - 1e-3;
    final hi = math.max(za, zb) + 1e-3;
    return zp >= lo && zp <= hi;
  }

  /// Poles of the two spherical caps bounded by [circle]. A great circle
  /// (centre at the sphere centre) has no "small" cap, so the poles are
  /// taken along the circle normal.
  static (Vector3, Vector3) _sphereCapPoles(BrepSphere sph, BrepCircle circle) {
    final toC = circle.center - sph.center;
    if (toC.length < sph.radius * 1e-6) {
      var n = circle.normal.clone();
      if (n.length < 1e-12) n = Vector3(0, 0, 1);
      n.normalize();
      final a = sph.center + n * sph.radius;
      final b = sph.center - n * sph.radius;
      return (a, b);
    }
    final small = sph.closestPoint(circle.center);
    return (small, sph.center * 2 - small);
  }

  /// Arcs of [circle] whose interior lies on [poly]. Empty when the circle
  /// misses the face or is a full interior hole (caller handles those).
  static List<BrepCircle> _clipCircleToPoly(
    BrepCircle circle,
    List<Vector3> poly,
    Vector3 n,
  ) {
    if (_circleIsInteriorHole(poly, circle, n)) return const [];
    final hits = <(double, Vector3)>[];
    void addHit(Vector3 p) {
      final t = circle.closestParam(p);
      for (final h in hits) {
        if ((h.$1 - t).abs() < 1e-6) return;
      }
      hits.add((t, p.clone()));
    }

    for (var i = 0; i < poly.length; i++) {
      for (final p in _circleSegmentHits(
        circle,
        poly[i],
        poly[(i + 1) % poly.length],
      )) {
        addHit(p);
      }
    }
    if (hits.length < 2) return const [];
    hits.sort((a, b) => a.$1.compareTo(b.$1));
    final out = <BrepCircle>[];
    for (var i = 0; i < hits.length; i++) {
      final a = hits[i];
      final b = hits[(i + 1) % hits.length];
      var t0 = a.$1;
      var t1 = b.$1;
      if (i == hits.length - 1) t1 += math.pi * 2;
      if (t1 - t0 < 1e-6) continue;
      final arc = BrepCircle(
        center: circle.center,
        normal: circle.normal,
        radius: circle.radius,
        xDir: circle.xDir,
        t0: t0,
        t1: t1,
      );
      final mid = arc.evaluate(0.5 * (t0 + t1));
      if (!_pointInPoly(poly, mid, n) && _minDistToEdges(mid, poly, n) > 1e-4) {
        continue;
      }
      out.add(arc);
    }
    return out;
  }

  static List<List<_ImprintSeg>> _chainCircleArcs(List<BrepCircle> arcs) {
    final segs = [
      for (final a in arcs) _ImprintSeg(a.evaluate(a.t0), a.evaluate(a.t1), a),
    ];
    return _chainImprintSegs(segs);
  }

  static bool _circleHitsPoly(
    BrepCircle circle,
    List<Vector3> poly,
    Vector3 n,
  ) {
    if (poly.length < 3) return false;
    const samples = 24;
    for (var i = 0; i < samples; i++) {
      final p = circle.evaluate(
        circle.t0 + (circle.t1 - circle.t0) * (i / samples),
      );
      if (_pointInPoly(poly, p, n) || _minDistToEdges(p, poly, n) <= 1e-4) {
        return true;
      }
    }
    return _pointInPoly(poly, circle.center, n) &&
        _minDistToEdges(circle.center, poly, n) <= circle.radius + 1e-4;
  }

  static void _collectCylinder(
    BrepSolid self,
    BrepFace face,
    BrepSolid other,
    int op, {
    required bool fromA,
    required List<_Piece> out,
  }) {
    final cyl = face.surface as BrepCylinder;
    final domain = BrepTessellate.cylinderFaceDomain(self, face);
    var z0 = domain.$3;
    var z1 = domain.$4;
    if (z1 <= z0) {
      z0 = 0;
      z1 = 1;
    }
    // Split the wall in UV by planes whose *faces* actually meet the cylinder.
    var rects = [(domain.$1, domain.$2, z0, z1)];
    var cutAtSeam = false;
    for (final of in other.faces) {
      if (of.surface is! BrepPlane) continue;
      final pl = of.surface as BrepPlane;
      final otherPoly = _planePoly(other, of);
      if (!_planeFaceHitsCylinder(pl, otherPoly, cyl, z0, z1)) continue;
      final an = pl.normal.dot(cyl.axis).abs();
      if (an > 0.999) {
        final zhit = (pl.origin - cyl.origin).dot(cyl.axis);
        if (zhit <= z0 + 1e-8 || zhit >= z1 - 1e-8) continue;
        final next = <(double, double, double, double)>[];
        for (final r in rects) {
          if (zhit > r.$3 + 1e-8 && zhit < r.$4 - 1e-8) {
            next.add((r.$1, r.$2, r.$3, zhit));
            next.add((r.$1, r.$2, zhit, r.$4));
          } else {
            next.add(r);
          }
        }
        rects = next;
      } else if (an < 1e-3) {
        final curves = BrepIntersect.surfaces(pl, cyl);
        for (final curve in curves) {
          if (curve is! BrepLine) continue;
          // The supporting plane can meet the cylinder along two generators.
          // A finite cutter face hitting one does not authorize splitting at
          // the other (which may run through a distant fillet patch).
          if (!_generatorHitsPoly(curve, otherPoly, pl.normal, cyl, z0, z1)) {
            continue;
          }
          final uv = cyl.parameters(curve.evaluate(0.5));
          if (uv == null) continue;
          var uCut = uv.$1;
          if (uCut < 0) uCut += math.pi * 2;
          uCut = uCut % (math.pi * 2);
          if (uCut < 1e-5 || (math.pi * 2 - uCut) < 1e-5) cutAtSeam = true;
          final next = <(double, double, double, double)>[];
          for (final r in rects) {
            // A patch domain may run past 2pi, so lift the cut into this
            // rect's turn instead of comparing it in [0, 2pi).
            var u = uCut;
            while (u < r.$1 - 1e-9) {
              u += math.pi * 2;
            }
            while (u >= r.$1 + math.pi * 2 - 1e-9) {
              u -= math.pi * 2;
            }
            if (u > r.$1 + 1e-6 && u < r.$2 - 1e-6) {
              next.add((r.$1, u, r.$3, r.$4));
              next.add((u, r.$2, r.$3, r.$4));
            } else {
              next.add(r);
            }
          }
          rects = next;
        }
      } else {
        // Oblique ellipse: plane faces imprint the curve. Keep this wall
        // whole so classification can drop or keep it.
      }
    }
    for (final of in other.faces) {
      if (of.surface is! BrepCylinder) continue;
      final otherCyl = of.surface as BrepCylinder;
      if (cyl.axis.cross(otherCyl.axis).length > 1e-6) continue;
      // Intersect finite walls, not their infinite supporting cylinders. A
      // distant boss must not split an older wall at a (possibly near-tangent)
      // generator. Express both trimmed axial spans in the same frame.
      final od = BrepTessellate.cylinderFaceDomain(other, of);
      final oz0 = (otherCyl.evaluate(0, od.$3) - cyl.origin).dot(cyl.axis);
      final oz1 = (otherCyl.evaluate(0, od.$4) - cyl.origin).dot(cyl.axis);
      if (math.max(oz0, oz1) < z0 - 1e-7 || math.min(oz0, oz1) > z1 + 1e-7) {
        continue;
      }
      final curves = BrepIntersect.surfaces(cyl, otherCyl);
      for (final curve in curves) {
        if (curve is! BrepLine) continue;
        final uv = cyl.parameters(curve.evaluate(0.5));
        if (uv == null) continue;
        var uCut = uv.$1;
        if (uCut < 0) uCut += math.pi * 2;
        uCut = uCut % (math.pi * 2);
        final next = <(double, double, double, double)>[];
        for (final r in rects) {
          var u = uCut;
          while (u < r.$1 - 1e-9) {
            u += math.pi * 2;
          }
          while (u >= r.$1 + math.pi * 2 - 1e-9) {
            u -= math.pi * 2;
          }
          if (u > r.$1 + 1e-6 && u < r.$2 - 1e-6) {
            next.add((r.$1, u, r.$3, r.$4));
            next.add((u, r.$2, r.$3, r.$4));
          } else {
            next.add(r);
          }
        }
        rects = next;
      }
    }
    for (final of in other.faces) {
      if (of.surface is! BrepSphere) continue;
      final selfBounds = conservativeSolidBounds(self);
      if (selfBounds != null &&
          !conservativeSphereFaceBounds(
            other,
            of,
          ).intersects(selfBounds, tolerance: 1e-6)) {
        continue;
      }
      for (final curve in BrepIntersect.surfaces(cyl, of.surface)) {
        if (curve is! BrepCircle) continue;
        final zhit = (curve.center - cyl.origin).dot(cyl.axis);
        if (zhit <= z0 + 1e-8 || zhit >= z1 - 1e-8) continue;
        final next = <(double, double, double, double)>[];
        for (final r in rects) {
          if (zhit > r.$3 + 1e-8 && zhit < r.$4 - 1e-8) {
            next.add((r.$1, r.$2, r.$3, zhit));
            next.add((r.$1, r.$2, zhit, r.$4));
          } else {
            next.add(r);
          }
        }
        rects = next;
      }
    }
    // Subtract keeps only the tool wall that sits inside A. A through-cut
    // taller than A must not leave the overhanging band as a dangling face.
    if (!fromA && op == subtractOp) {
      final span = _axisSpan(other, cyl);
      if (span != null) {
        final clipped = <(double, double, double, double)>[];
        for (final r in rects) {
          final z0c = math.max(r.$3, span.$1);
          final z1c = math.min(r.$4, span.$2);
          if (z1c - z0c > 1e-8) clipped.add((r.$1, r.$2, z0c, z1c));
        }
        rects = clipped;
      }
    }
    final keptRects = _mergeCylinderRects(
      rects.where((r) {
        final u = 0.5 * (r.$1 + r.$2);
        final v = 0.5 * (r.$3 + r.$4);
        final p = cyl.evaluate(u, v);
        var n = cyl.normalAt(u, v);
        if (face.reversed) n = -n;
        return _keepNeighborhood(
          op,
          fromA: fromA,
          other: other,
          centroid: p,
          outward: n,
          cylinderWall: true,
        );
      }),
      wrapSeam: !cutAtSeam,
    );
    if (keptRects.length == 1) {
      final r = keptRects.single;
      if ((r.$1 - domain.$1).abs() < 1e-9 &&
          (r.$2 - domain.$2).abs() < 1e-9 &&
          (r.$3 - domain.$3).abs() < 1e-9 &&
          (r.$4 - domain.$4).abs() < 1e-9) {
        out.add(
          _Piece(
            cyl,
            const [],
            !fromA && op == subtractOp ? !face.reversed : face.reversed,
            sourceSolid: self,
            sourceFace: face,
          ),
        );
        return;
      }
    }
    for (final r in keptRects) {
      var rev = face.reversed;
      if (!fromA && op == subtractOp) rev = !rev;
      final poly = [
        cyl.evaluate(r.$1, r.$3),
        cyl.evaluate(r.$2, r.$3),
        cyl.evaluate(r.$2, r.$4),
        cyl.evaluate(r.$1, r.$4),
      ];
      out.add(_Piece(cyl, poly, rev, uv: r));
    }
  }

  /// Join cylinder patches that wrap across the u=0 seam. Adjacent patches
  /// split by a cutting generator stay separate so that generator still sews.
  static List<(double, double, double, double)> _mergeCylinderRects(
    Iterable<(double, double, double, double)> rects, {
    required bool wrapSeam,
  }) {
    final list = rects.toList();
    if (list.length <= 1) return list;
    final byZ = <String, List<(double, double, double, double)>>{};
    for (final r in list) {
      final k = '${r.$3.toStringAsFixed(6)}_${r.$4.toStringAsFixed(6)}';
      (byZ[k] ??= []).add(r);
    }
    final out = <(double, double, double, double)>[];
    for (final group in byZ.values) {
      group.sort((a, b) => a.$1.compareTo(b.$1));
      if (wrapSeam &&
          group.length >= 2 &&
          group.first.$1.abs() < 1e-5 &&
          (group.last.$2 - math.pi * 2).abs() < 1e-5) {
        final first = group.removeAt(0);
        final last = group.removeLast();
        group.add((last.$1, first.$2 + math.pi * 2, first.$3, first.$4));
      }
      out.addAll(group);
    }
    return out;
  }

  static bool _needsTorusBoundary(BrepSolid a, BrepSolid b) {
    for (final pair in [(a, b), (b, a)]) {
      for (final face in pair.$1.faces) {
        final tor = face.surface;
        if (tor is! BrepTorus) continue;
        for (final other in pair.$2.faces) {
          final plane = other.surface;
          if (plane is! BrepPlane) continue;
          final dot = plane.normal.dot(tor.axis).abs();
          if (dot > 0.999 ||
              dot < 1e-3 && plane.height(tor.center).abs() <= 1e-6)
            continue;
          if (_planeHitsTorus(plane, tor)) return true;
        }
      }
    }
    return false;
  }

  static void _collectTorus(
    BrepSolid self,
    BrepFace face,
    BrepSolid other,
    int op, {
    required bool fromA,
    required List<_Piece> out,
  }) {
    final tor = face.surface as BrepTorus;
    for (final of in other.faces) {
      if (of.surface is! BrepPlane) continue;
      final pl = of.surface as BrepPlane;
      if (pl.normal.dot(tor.axis).abs() > 0.999) continue;
      // A plane that contains the axis splits the tube at two meridians.
      if (pl.normal.dot(tor.axis).abs() < 1e-3 &&
          pl.height(tor.center).abs() <= 1e-4) {
        continue;
      }
      if (_planeHitsTorus(pl, tor)) {
        throw BrepUnsupported(
          'Torus–plane intersection (non-principal) is not implemented yet',
        );
      }
    }
    // Fillets occupy only a trimmed band of their supporting torus. Starting
    // from the full surface resurrects discarded tube quadrants, including
    // coincident rims that later look like collapsed disconnected shells.
    var rects = [
      BrepTessellate.torusFaceDomain(self, tor, face) ??
          (0.0, math.pi * 2, 0.0, math.pi * 2),
    ];
    for (final of in other.faces) {
      if (of.surface is! BrepPlane) continue;
      final pl = of.surface as BrepPlane;
      if (pl.normal.dot(tor.axis).abs() < 0.999) continue;
      for (final curve in BrepIntersect.surfaces(pl, tor)) {
        if (curve is! BrepCircle) continue;
        final uv = tor.parameters(curve.evaluate(0));
        if (uv == null) continue;
        var vCut = uv.$2;
        while (vCut < 0) {
          vCut += math.pi * 2;
        }
        vCut %= math.pi * 2;
        if (vCut < 1e-5 || (math.pi * 2 - vCut) < 1e-5) continue;
        final next = <(double, double, double, double)>[];
        for (final r in rects) {
          var v = vCut;
          while (v < r.$3 - 1e-9) {
            v += math.pi * 2;
          }
          while (v >= r.$3 + math.pi * 2 - 1e-9) {
            v -= math.pi * 2;
          }
          if (v > r.$3 + 1e-6 && v < r.$4 - 1e-6) {
            next.add((r.$1, r.$2, r.$3, v));
            next.add((r.$1, r.$2, v, r.$4));
          } else {
            next.add(r);
          }
        }
        rects = next;
      }
    }
    for (final of in other.faces) {
      if (of.surface is! BrepPlane) continue;
      final pl = of.surface as BrepPlane;
      if (pl.normal.dot(tor.axis).abs() > 1e-3) continue;
      if (pl.height(tor.center).abs() > 1e-4) continue;
      for (final curve in BrepIntersect.surfaces(pl, tor)) {
        if (curve is! BrepCircle) continue;
        final uv = tor.parameters(curve.center);
        if (uv == null) continue;
        var uCut = uv.$1;
        if (uCut < 0) uCut += math.pi * 2;
        uCut %= math.pi * 2;
        if (uCut < 1e-5 || (math.pi * 2 - uCut) < 1e-5) continue;
        final next = <(double, double, double, double)>[];
        for (final r in rects) {
          var u = uCut;
          while (u < r.$1 - 1e-9) {
            u += math.pi * 2;
          }
          while (u >= r.$1 + math.pi * 2 - 1e-9) {
            u -= math.pi * 2;
          }
          if (u > r.$1 + 1e-6 && u < r.$2 - 1e-6) {
            next.add((r.$1, u, r.$3, r.$4));
            next.add((u, r.$2, r.$3, r.$4));
          } else {
            next.add(r);
          }
        }
        rects = next;
      }
    }
    for (final r in rects) {
      final u = 0.5 * (r.$1 + r.$2);
      final v = 0.5 * (r.$3 + r.$4);
      final p = tor.evaluate(u, v);
      var n = tor.normalAt(u, v);
      if (face.reversed) n = -n;
      if (!_keepNeighborhood(
        op,
        fromA: fromA,
        other: other,
        centroid: p,
        outward: n,
      )) {
        continue;
      }
      var rev = face.reversed;
      if (!fromA && op == subtractOp) rev = !rev;
      out.add(
        _Piece(
          tor,
          [
            tor.evaluate(r.$1, r.$3),
            tor.evaluate(r.$2, r.$3),
            tor.evaluate(r.$2, r.$4),
            tor.evaluate(r.$1, r.$4),
          ],
          rev,
          uv: r,
        ),
      );
    }
  }

  static bool _planeHitsTorus(BrepPlane pl, BrepTorus tor) {
    var pos = false;
    var neg = false;
    const n = 16;
    for (var i = 0; i < n; i++) {
      for (var j = 0; j < n; j++) {
        final h = pl.height(
          tor.evaluate(i / n * math.pi * 2, j / n * math.pi * 2),
        );
        if (h > 1e-4) pos = true;
        if (h < -1e-4) neg = true;
        if (pos && neg) return true;
      }
    }
    return false;
  }

  static bool _planeHitsCone(BrepPlane pl, BrepCone cone) {
    var pos = false;
    var neg = false;
    const n = 16;
    for (var i = 0; i < n; i++) {
      for (var j = 0; j <= n; j++) {
        final z = cone.height * (j / n);
        final h = pl.height(cone.evaluate(i / n * math.pi * 2, z));
        if (h > 1e-4) pos = true;
        if (h < -1e-4) neg = true;
        if (pos && neg) return true;
      }
    }
    return false;
  }

  static void _collectCone(
    BrepSolid self,
    BrepFace face,
    BrepSolid other,
    int op, {
    required bool fromA,
    required List<_Piece> out,
  }) {
    final cone = face.surface as BrepCone;
    var z0 = 0.0;
    var z1 = cone.height;
    var rects = [(0.0, math.pi * 2, z0, z1)];
    for (final of in other.faces) {
      if (of.surface is! BrepPlane) continue;
      // A finite polygon entirely inside the supporting cone cannot cut its
      // lateral wall. Extending a socket's side faces to infinite planes would
      // invent oblique intersections far outside the actual cutting tool.
      // The cone interior is convex, so testing every straight trim endpoint
      // proves containment of the whole polygon (including concave polygons).
      if (of.loops.isNotEmpty &&
          of.loops.every((l) => l.coedges.isNotEmpty) &&
          of.loops.every(
            (l) => l.coedges.every((co) {
              final e = other.edges[co.edge];
              if (e.curve is! BrepLine) return false;
              return [other.vertices[e.v0], other.vertices[e.v1]].every((p) {
                final d = p - cone.origin, z = d.dot(cone.axis);
                return (d - cone.axis * z).length < cone.radiusAt(z) - 1e-7;
              });
            }),
          )) {
        continue;
      }
      final pl = of.surface as BrepPlane;
      final an = pl.normal.dot(cone.axis).abs();
      if (an > 0.999) {
        final den = pl.normal.dot(cone.axis);
        if (den.abs() < 1e-12) continue;
        final zhit = -pl.height(cone.origin) / den;
        if (zhit <= z0 + 1e-8 || zhit >= z1 - 1e-8) continue;
        final next = <(double, double, double, double)>[];
        for (final r in rects) {
          if (zhit > r.$3 + 1e-8 && zhit < r.$4 - 1e-8) {
            next.add((r.$1, r.$2, r.$3, zhit));
            next.add((r.$1, r.$2, zhit, r.$4));
          } else {
            next.add(r);
          }
        }
        rects = next;
      } else {
        final curves = BrepIntersect.surfaces(pl, cone);
        final ellipses = [
          for (final c in curves)
            if (c is BrepEllipse) c,
        ];
        if (ellipses.isEmpty && _planeHitsCone(pl, cone)) {
          throw BrepUnsupported(
            'Cone–plane intersection (oblique) is not implemented yet',
          );
        }
        for (final ellipse in ellipses) {
          final apex = cone.origin + cone.axis * cone.height;
          final towardApex = cone.closestPoint(
            ellipse.center + (apex - ellipse.center) * 0.15,
          );
          final towardBase = cone.closestPoint(
            ellipse.center + (cone.origin - ellipse.center) * 0.15,
          );
          for (final probe in [towardApex, towardBase]) {
            var nrm = cone.normalAt(
              cone.parameters(probe)?.$1 ?? 0,
              cone.parameters(probe)?.$2 ?? 0,
            );
            if (face.reversed) nrm = -nrm;
            if (!_keepNeighborhood(
              op,
              fromA: fromA,
              other: other,
              centroid: probe,
              outward: nrm,
            )) {
              continue;
            }
            var rev = face.reversed;
            if (!fromA && op == subtractOp) rev = !rev;
            final wantApex =
                (probe - apex).length < (probe - cone.origin).length;
            out.add(
              _Piece(
                cone,
                const [],
                rev,
                outerCurves: [ellipse],
                holeCircles: wantApex
                    ? const []
                    : [
                        BrepCircle(
                          center: cone.origin,
                          normal: cone.axis,
                          radius: cone.radius,
                          xDir: cone.xDir,
                        ),
                      ],
              ),
            );
          }
        }
        if (ellipses.isNotEmpty) {
          rects = [];
        }
      }
    }
    if (!fromA && op == subtractOp) {
      // Clip hanging cone wall to the other solid's axis span.
      var lo = double.infinity;
      var hi = -double.infinity;
      for (final v in other.vertices) {
        final z = (v - cone.origin).dot(cone.axis);
        lo = math.min(lo, z);
        hi = math.max(hi, z);
      }
      if (lo < hi) {
        final clipped = <(double, double, double, double)>[];
        for (final r in rects) {
          final z0c = math.max(r.$3, lo);
          final z1c = math.min(r.$4, hi);
          if (z1c - z0c > 1e-8) clipped.add((r.$1, r.$2, z0c, z1c));
        }
        rects = clipped;
      }
    }
    for (final r in rects) {
      final u = 0.5 * (r.$1 + r.$2);
      final v = 0.5 * (r.$3 + r.$4);
      final p = cone.evaluate(u, v);
      var n = cone.normalAt(u, v);
      if (face.reversed) n = -n;
      if (!_keepNeighborhood(
        op,
        fromA: fromA,
        other: other,
        centroid: p,
        outward: n,
      )) {
        continue;
      }
      var rev = face.reversed;
      if (!fromA && op == subtractOp) rev = !rev;
      out.add(
        _Piece(
          cone,
          [
            cone.evaluate(r.$1, r.$3),
            cone.evaluate(r.$2, r.$3),
            cone.evaluate(r.$2, r.$4),
            cone.evaluate(r.$1, r.$4),
          ],
          rev,
          uv: r,
        ),
      );
    }
  }

  static bool _inOther(BrepSide side) => side == BrepSide.inside;

  static bool _inResult(
    int op, {
    required bool fromA,
    required bool inSelf,
    required BrepSide other,
  }) {
    final inO = _inOther(other);
    switch (op) {
      case subtractOp:
        return fromA ? inSelf && !inO : inSelf && inO;
      case intersectOp:
        return inSelf && inO;
      default:
        return inSelf || inO;
    }
  }

  static bool _keepNeighborhood(
    int op, {
    required bool fromA,
    required BrepSolid other,
    required Vector3 centroid,
    required Vector3 outward,
    bool hasHole = false,
    bool cylinderWall = false,
  }) {
    if (op == subtractOp) {
      final side = BrepClassify.side(other, centroid);
      if (fromA) {
        if (side == BrepSide.inside) return false;
        if (side == BrepSide.outside) return true;
        // Coplanar with the tool: keep a punched ring, drop a fully
        // overwritten face (e.g. a cylinder cap sitting on a box top).
        return hasHole;
      }
      if (side == BrepSide.inside) return true;
      // Through-hole wall can lie on A's face (inscribed cylinder).
      return side == BrepSide.on && cylinderWall;
    }
    if (op == intersectOp) {
      final side = BrepClassify.side(other, centroid);
      return side == BrepSide.inside || (fromA && side == BrepSide.on);
    }
    var n = outward.clone();
    if (n.length < 1e-12) {
      final side = BrepClassify.side(other, centroid);
      switch (op) {
        case subtractOp:
          return fromA ? side != BrepSide.inside : side == BrepSide.inside;
        case intersectOp:
          return side == BrepSide.inside || (fromA && side == BrepSide.on);
        default:
          return side == BrepSide.outside;
      }
    }
    n.normalize();
    const eps = 1e-4;
    final inward = BrepClassify.side(other, centroid - n * eps);
    final outS = BrepClassify.side(other, centroid + n * eps);
    final keep =
        _inResult(op, fromA: fromA, inSelf: true, other: inward) !=
        _inResult(op, fromA: fromA, inSelf: false, other: outS);
    if (!keep) return false;
    // Prefer A's copy of a coplanar envelope face so union does not double-skin.
    if (!fromA &&
        op != subtractOp &&
        BrepClassify.side(other, centroid) == BrepSide.on) {
      return false;
    }
    return true;
  }

  /// Inner loops already on [face] (through-holes from earlier cuts).
  static ({
    List<List<Vector3>> holes,
    List<List<BrepCurve>> curves,
    List<BrepCircle> circles,
  })
  _faceHoles(BrepSolid s, BrepFace face) {
    final holes = <List<Vector3>>[];
    final curves = <List<BrepCurve>>[];
    final circles = <BrepCircle>[];
    for (var i = 1; i < face.loops.length; i++) {
      final loop = face.loops[i];
      if (loop.coedges.isEmpty) continue;
      if (loop.coedges.length == 1) {
        final e = s.edges[loop.coedges.first.edge];
        final c = e.curve;
        if (c is BrepCircle && c.isClosed) {
          circles.add(c);
          continue;
        }
      }
      final pts = <Vector3>[];
      final cs = <BrepCurve>[];
      for (final co in loop.coedges) {
        final e = s.edges[co.edge];
        pts.add(s.vertices[co.reversed ? e.v1 : e.v0].clone());
        cs.add(e.curve);
      }
      if (pts.length >= 2) {
        holes.add(pts);
        curves.add(cs);
      }
    }
    return (holes: holes, curves: curves, circles: circles);
  }

  static bool _holeSamplesInPoly(
    List<Vector3> poly,
    List<Vector3> hole,
    Vector3 n,
  ) {
    if (hole.isEmpty) return false;
    var c = Vector3.zero();
    for (final p in hole) {
      c += p;
    }
    c /= hole.length.toDouble();
    return _pointInPoly(poly, c, n);
  }

  /// True when [circle] sits strictly inside [poly] as an inner loop.
  static bool _circleIsInteriorHole(
    List<Vector3> poly,
    BrepCircle circle,
    Vector3 n,
  ) {
    if (_minDistToEdges(circle.center, poly, n) <= circle.radius + 1e-4) {
      return false;
    }
    if (!_pointInPoly(poly, circle.center, n)) return false;
    const samples = 16;
    for (var i = 0; i < samples; i++) {
      final p = circle.evaluate(
        circle.t0 + (circle.t1 - circle.t0) * (i / samples),
      );
      if (!_pointInPoly(poly, p, n)) return false;
    }
    return true;
  }

  static bool _ellipseIsInteriorHole(
    List<Vector3> poly,
    BrepEllipse e,
    Vector3 n,
  ) {
    final reach = math.max(e.rMajor, e.rMinor);
    if (_minDistToEdges(e.center, poly, n) <= reach + 1e-4) {
      return false;
    }
    if (!_pointInPoly(poly, e.center, n)) return false;
    const samples = 16;
    for (var i = 0; i < samples; i++) {
      final p = e.evaluate(e.t0 + (e.t1 - e.t0) * (i / samples));
      if (!_pointInPoly(poly, p, n)) return false;
    }
    return true;
  }

  /// Open a circular hole that the split boundary cut through: replace the
  /// chord across the disk with the same analytic arc the cylinder wall uses.
  static ({List<Vector3> poly, BrepCircle arc})? _notchPolyWithCircle(
    List<Vector3> poly,
    BrepCircle circle,
    Vector3 n,
  ) {
    final hits = <Vector3>[];
    void addHit(Vector3 p) {
      for (final h in hits) {
        if ((h - p).length < 1e-4) return;
      }
      hits.add(p.clone());
    }

    for (var i = 0; i < poly.length; i++) {
      for (final p in _circleSegmentHits(
        circle,
        poly[i],
        poly[(i + 1) % poly.length],
      )) {
        addHit(p);
      }
    }
    if (hits.length != 2) return null;

    BrepCircle arcBetween(double ta, double tb) {
      var a = ta;
      var b = tb;
      while (b < a) {
        b += math.pi * 2;
      }
      return BrepCircle(
        center: circle.center,
        normal: circle.normal,
        radius: circle.radius,
        xDir: circle.xDir,
        t0: a,
        t1: b,
      );
    }

    final tA = circle.closestParam(hits[0]);
    final tB = circle.closestParam(hits[1]);
    final a1 = arcBetween(tA, tB);
    final a2 = arcBetween(tB, tA);
    bool interiorMid(BrepCircle arc) {
      final mid = arc.evaluate(0.5 * (arc.t0 + arc.t1));
      return _pointInPoly(poly, mid, n) && _minDistToEdges(mid, poly, n) > 1e-4;
    }

    final BrepCircle arc;
    final in1 = interiorMid(a1);
    final in2 = interiorMid(a2);
    if (in1 && !in2) {
      arc = a1;
    } else if (in2 && !in1) {
      arc = a2;
    } else {
      // Both hits on one edge: the in-face test can see both midpoints as
      // outside (center sits off the face). Keep the arc nearer the face.
      final c = _centroid(poly);
      final d1 = (a1.evaluate(0.5 * (a1.t0 + a1.t1)) - c).length;
      final d2 = (a2.evaluate(0.5 * (a2.t0 + a2.t1)) - c).length;
      arc = d1 <= d2 ? a1 : a2;
    }

    var withHits = _insertOnBoundary(poly, hits[0], n);
    if (withHits == null) return null;
    withHits = _insertOnBoundary(withHits, hits[1], n);
    if (withHits == null) return null;
    final verts = withHits;
    final i0 = _nearestVert(verts, hits[0]);
    final i1 = _nearestVert(verts, hits[1]);
    if (i0 == i1) return null;

    List<Vector3> walk(int from, int to) {
      final out = <Vector3>[];
      var i = from;
      for (var k = 0; k <= verts.length; k++) {
        out.add(verts[i].clone());
        if (i == to) break;
        i = (i + 1) % verts.length;
      }
      return out;
    }

    final pathA = walk(i0, i1);
    final pathB = walk(i1, i0);
    bool throughDisk(List<Vector3> path) {
      for (var i = 0; i < path.length - 1; i++) {
        final m = (path[i] + path[i + 1]) * 0.5;
        if ((m - circle.center).length < circle.radius - 1e-4) return true;
      }
      return false;
    }

    final List<Vector3> keep;
    final aHit = throughDisk(pathA);
    final bHit = throughDisk(pathB);
    if (aHit && !bHit) {
      keep = pathB;
    } else if (bHit && !aHit) {
      keep = pathA;
    } else {
      return null;
    }
    if (keep.length < 3) return null;
    return (poly: keep, arc: arc);
  }

  static List<Vector3> _circleSegmentHits(
    BrepCircle circle,
    Vector3 a,
    Vector3 b,
  ) {
    final d = b - a;
    final len2 = d.length2;
    if (len2 < 1e-20) return const [];
    final f = a - circle.center;
    final A = d.dot(d);
    final bCoef = 2 * f.dot(d);
    final C = f.dot(f) - circle.radius * circle.radius;
    var disc = bCoef * bCoef - 4 * A * C;
    // Roundoff at tangency must not invent two nearby intersection points.
    // Those split a single contact into a microscopic gap after rotation.
    if (disc.abs() <= 1e-14 * (bCoef * bCoef + (4 * A * C).abs())) disc = 0;
    if (disc < 0) {
      if (disc > -1e-8) {
        disc = 0;
      } else {
        return const [];
      }
    }
    final s = math.sqrt(disc);
    final out = <Vector3>[];
    for (final t in [(-bCoef - s) / (2 * A), (-bCoef + s) / (2 * A)]) {
      if (t < -1e-6 || t > 1 + 1e-6) continue;
      final p = a + d * t.clamp(0.0, 1.0);
      if (((p - circle.center).length - circle.radius).abs() > 1e-3) continue;
      if ((p - circle.center).dot(circle.normal).abs() > 1e-3) continue;
      out.add(p);
    }
    return out;
  }

  /// Caps with holes or circular edges cannot be split as line polygons.
  /// Union/intersect only imprint onto those; subtract still imprints
  /// onto simple faces so a D-notch can cut a box.
  static bool _faceNeedsPrismImprint(BrepSolid s, BrepFace f) {
    if (f.loops.length > 1) return true;
    if (f.loops.isEmpty) return false;
    for (final c in f.loops.first.coedges) {
      if (s.edges[c.edge].curve is! BrepLine) return true;
    }
    return false;
  }

  static bool _segsHaveCurve(List<_ImprintSeg> segs) {
    for (final s in segs) {
      if (s.curve is! BrepLine) return true;
    }
    return false;
  }

  /// Imprint a prism tool onto a planar cap: pocket (inner loop), side
  /// notch (outer-loop splice), or a full wipe. Full cylinders stay on
  /// [_cutPlaneByCylinder] so inscribed / through-hole paths do not change.
  static _PrismImprint? _imprintPrismOnPlane(
    BrepSolid self,
    BrepFace face,
    BrepSolid tool,
    int op,
    ({
      List<List<Vector3>> holes,
      List<List<BrepCurve>> curves,
      List<BrepCircle> circles,
    })
    priorHoles,
  ) {
    final plane = face.surface as BrepPlane;
    final axis = plane.normal;
    if (!_prismMeetsPlane(tool, plane)) return null;
    final segs = <({Vector3 a, Vector3 b, BrepCurve curve})>[];
    final aligned = _isPrismAlong(tool, axis);
    if (aligned) {
      for (final tf in tool.faces) {
        final seg = _wallSectionOnPlane(tool, tf, plane, axis);
        if (seg != null) segs.add(seg);
      }
    } else {
      final section = _planarSolidSection(tool, plane);
      if (section == null) return null;
      segs.addAll(section);
    }
    final subject = _planeOuterSegs(self, face);
    if (subject.length < 3) return null;
    final List<_ImprintSeg> clip;
    if (segs.length == 1 &&
        segs.first.curve is BrepCircle &&
        (segs.first.curve as BrepCircle).isClosed) {
      // A cylinder straddling the face boundary. Split the rim where it
      // crosses, so the bite stays an analytic arc. Fully inside /
      // containing disks find no crossing and stay on
      // [_cutPlaneByCylinder].
      final arcs = _closedCircleToArcs(
        segs.first.curve as BrepCircle,
        subject,
        plane.normal,
      );
      if (arcs == null) return null;
      clip = arcs;
    } else {
      if (segs.length < 2) return null;
      final chain = _largestClosedChain(segs);
      if (chain == null) return null;
      // A section with multiple loops needs a multi-region overlay. Do not
      // silently discard smaller loops (including holes) in this path.
      if (!aligned && chain.length != segs.length) return null;
      clip = [
        for (final s in chain) _ImprintSeg(s.a.clone(), s.b.clone(), s.curve),
      ];
    }
    final poly = _segsToPoly(subject);
    final samples = <Vector3>[];
    for (final s in clip) {
      samples.addAll(_sampleSection(s.curve));
    }
    if (samples.length < 3 || poly.length < 3) return null;
    final n = plane.normal;
    final clipPoly = _segsToPoly(clip);
    if (!_profileOverlapsFace(poly, clipPoly, subject, clip, n, priorHoles)) {
      return _PrismImprint(
        outers: [
          _ImprintOuter(
            pts: [for (final s in subject) s.a.clone()],
            curves: [for (final s in subject) s.curve],
          ),
        ],
      );
    }
    // A box is a prism along every axis. Imprinting its profile onto every
    // simple face independently is not a 3D boolean: through-cuts,
    // overlapping boxes, and slice half-spaces all outspan the face with
    // line walls, and plane splits handle those. A face that carries arcs
    // or holes must still imprint: plane splits would sample its arcs into
    // polylines that no longer sew to the cylinder wall.
    if (_profileOutspansPoly(poly, samples, n) &&
        !_segsHaveCurve(clip) &&
        !_faceNeedsPrismImprint(self, face)) {
      return null;
    }
    final wholePts = [for (final s in subject) s.a.clone()];
    final wholeCurves = [for (final s in subject) s.curve];
    final covers = _profileCoversPoly(poly, samples, n);
    final inside = _strictlyInside(poly, samples, n);
    if (op == intersectOp) {
      if (covers) {
        return _PrismImprint(
          outers: [_ImprintOuter(pts: wholePts, curves: wholeCurves)],
        );
      }
      if (inside) {
        return _PrismImprint(
          outers: [
            _ImprintOuter(
              pts: [for (final s in clip) s.a.clone()],
              curves: [for (final s in clip) s.curve],
            ),
          ],
        );
      }
      final overlay = _overlayLoops(subject, clip, n, difference: false);
      if (overlay == null || overlay.isEmpty) return null;
      return _PrismImprint(outers: overlay);
    }
    if (covers) {
      // Half-space / slice tools swallow the face and stick out of it.
      // Dropping the face would delete a cap; fall through to plane splits.
      if (_profileOutspansPoly(poly, samples, n)) return null;
      return const _PrismImprint(drop: true);
    }
    if (inside) {
      return _PrismImprint(
        outers: [
          _ImprintOuter(
            pts: wholePts,
            curves: wholeCurves,
            holes: [
              [for (final s in clip) s.a.clone()],
            ],
            holeCurves: [
              [for (final s in clip) s.curve],
            ],
          ),
        ],
      );
    }
    final overlay = _overlayLoops(subject, clip, n, difference: true);
    if (overlay == null || overlay.isEmpty) return null;
    // A through-slot on a simple cap splits it into two outers. Plane
    // splits keep the walls; independent 2D leftovers do not.
    if (overlay.length > 1 && !_faceNeedsPrismImprint(self, face)) {
      return null;
    }
    if (op == unionOp) {
      // Both bodies imprint each other, so both would subtract the shared
      // region and nothing would cover it. Emit it as its own region;
      // [_keepNeighborhood] keeps A's copy and drops B's when the two
      // faces are coplanar, and drops both when the tool really covers it.
      final within = _overlayLoops(subject, clip, n, difference: false);
      if (within != null && within.isNotEmpty) {
        return _PrismImprint(outers: [...overlay, ...within]);
      }
    }
    return _PrismImprint(outers: overlay);
  }

  /// A closed rim split into arcs where it crosses [subject]. Null when the
  /// circle never crosses the boundary (through-hole, tangent, or covering).
  static List<_ImprintSeg>? _closedCircleToArcs(
    BrepCircle circle,
    List<_ImprintSeg> subject,
    Vector3 n,
  ) {
    const twoPi = math.pi * 2;
    final params = <double>[];
    for (final s in subject) {
      if (s.curve is! BrepLine) continue;
      for (final p in _circleSegmentHits(circle, s.a, s.b)) {
        if (!_pointOnSeg(s, p, n)) continue;
        var t = circle.closestParam(p);
        while (t < circle.t0 - 1e-9) {
          t += twoPi;
        }
        while (t >= circle.t0 + twoPi - 1e-9) {
          t -= twoPi;
        }
        var dup = false;
        for (final q in params) {
          if ((q - t).abs() < 1e-6) {
            dup = true;
            break;
          }
        }
        if (!dup) params.add(t);
      }
    }
    if (params.length < 2) return null;
    // Tangent rims (an inscribed cylinder touching all four walls) meet the
    // boundary without crossing it. Those keep the corner-cap handling in
    // [_cutPlaneByCylinder]; only a rim with material on both sides splits.
    final poly = _segsToPoly(subject);
    var anyIn = false;
    var anyOut = false;
    for (var i = 0; i < 32; i++) {
      final p = circle.evaluate(circle.t0 + twoPi * i / 32);
      if (_minDistToEdges(p, poly, n) < 1e-3) continue;
      if (_pointInPoly(poly, p, n)) {
        anyIn = true;
      } else {
        anyOut = true;
      }
    }
    if (!anyIn || !anyOut) return null;
    params.sort();
    final out = <_ImprintSeg>[];
    for (var i = 0; i < params.length; i++) {
      final a = params[i];
      final b = i == params.length - 1 ? params.first + twoPi : params[i + 1];
      if (b - a < 1e-6) continue;
      final arc = BrepCircle(
        center: circle.center,
        normal: circle.normal,
        radius: circle.radius,
        xDir: circle.xDir,
        t0: a,
        t1: b,
      );
      out.add(_ImprintSeg(arc.evaluate(a), arc.evaluate(b), arc));
    }
    return out.length < 2 ? null : out;
  }

  static List<_ImprintSeg> _planeOuterSegs(BrepSolid s, BrepFace f) {
    final segs = <_ImprintSeg>[];
    if (f.loops.isEmpty) return segs;
    for (final c in f.loops.first.coedges) {
      final e = s.edges[c.edge];
      final a = s.vertices[c.reversed ? e.v1 : e.v0];
      final b = s.vertices[c.reversed ? e.v0 : e.v1];
      if ((a - b).length < 1e-9 &&
          !(e.curve is BrepCircle && (e.curve as BrepCircle).isClosed)) {
        continue;
      }
      segs.add(_ImprintSeg(a.clone(), b.clone(), e.curve));
    }
    return segs;
  }

  static List<Vector3> _segsToPoly(List<_ImprintSeg> segs) {
    final pts = <Vector3>[];
    for (final s in segs) {
      pts.add(s.a.clone());
      if (s.curve is BrepCircle && !(s.curve as BrepCircle).isClosed) {
        final samples = _sampleSection(s.curve);
        for (var i = 1; i < samples.length - 1; i++) {
          pts.add(samples[i]);
        }
      }
    }
    return pts;
  }

  /// Subject minus clip (difference) or subject ∩ clip, keeping analytic
  /// line/circle edges.
  /// A point just inside [poly] alongside [s], or null when [s] is not on
  /// the boundary of [poly].
  static Vector3? _inwardProbe(_ImprintSeg s, List<Vector3> poly, Vector3 n) {
    final d = s.b - s.a;
    if (d.length2 < 1e-18) return null;
    var side = n.cross(d);
    if (side.length2 < 1e-18) return null;
    side = side.normalized() * 1e-3;
    final plus = s.mid + side;
    if (_pointInPoly(poly, plus, n)) return plus;
    final minus = s.mid - side;
    if (_pointInPoly(poly, minus, n)) return minus;
    return null;
  }

  static List<_ImprintOuter>? _overlayLoops(
    List<_ImprintSeg> subject,
    List<_ImprintSeg> clip,
    Vector3 n, {
    required bool difference,
  }) {
    final subPoly = _segsToPoly(subject);
    final clipPoly = _segsToPoly(clip);
    final hits = <Vector3>[];
    void addHit(Vector3 p) {
      for (final h in hits) {
        if ((h - p).length < 1e-4) return;
      }
      hits.add(p.clone());
    }

    for (final a in subject) {
      for (final b in clip) {
        for (final p in _segIntersect(a, b, n)) {
          addHit(p);
        }
      }
    }
    for (final s in subject) {
      addHit(s.a);
      addHit(s.b);
    }
    for (final s in clip) {
      addHit(s.a);
      addHit(s.b);
    }
    final splitA = _splitSegsAtHits(subject, hits, n);
    final splitB = _splitSegsAtHits(clip, hits, n);
    final kept = <_ImprintSeg>[];
    for (final s in splitA) {
      if (_midOnBoundary(s, clipPoly, n)) {
        // Both loops run along this segment, so both sides drop it and the
        // intersection can never close. Keep the subject's copy when the
        // material really is shared here.
        if (!difference) {
          final inw = _inwardProbe(s, subPoly, n);
          if (inw != null && _pointInPoly(clipPoly, inw, n)) kept.add(s);
        }
        continue;
      }
      final inClip = _pointInPoly(clipPoly, s.mid, n);
      if (difference ? inClip : !inClip) continue;
      kept.add(s);
    }
    for (final s in splitB) {
      if (_midOnBoundary(s, subPoly, n)) continue;
      if (!_pointInPoly(subPoly, s.mid, n)) continue;
      kept.add(difference ? _ImprintSeg(s.b.clone(), s.a.clone(), s.curve) : s);
    }
    if (kept.isEmpty) return null;
    final loops = _chainImprintSegs(kept);
    if (loops.isEmpty) return null;

    final sampled = [for (final loop in loops) _segsToPoly(loop)];
    final innerOf = List<int>.filled(loops.length, -1);
    for (var i = 0; i < loops.length; i++) {
      if (sampled[i].length < 3) continue;
      final c = _centroid(sampled[i]);
      var best = -1;
      var bestArea = double.infinity;
      for (var j = 0; j < loops.length; j++) {
        if (i == j || sampled[j].length < 3) continue;
        if (!_pointInPoly(sampled[j], c, n)) continue;
        final area = _loopArea(sampled[j], n).abs();
        if (area < bestArea) {
          bestArea = area;
          best = j;
        }
      }
      innerOf[i] = best;
    }
    final outers = <_ImprintOuter>[];
    for (var i = 0; i < loops.length; i++) {
      if (innerOf[i] >= 0) continue;
      final holes = <List<Vector3>>[];
      final holeCurves = <List<BrepCurve>>[];
      for (var j = 0; j < loops.length; j++) {
        if (innerOf[j] != i) continue;
        holes.add([for (final s in loops[j]) s.a.clone()]);
        holeCurves.add([for (final s in loops[j]) s.curve]);
      }
      final outerArcs = <BrepCircle>[
        for (final s in loops[i])
          if (s.curve is BrepCircle && !(s.curve as BrepCircle).isClosed)
            s.curve as BrepCircle,
      ];
      outers.add(
        _ImprintOuter(
          pts: [for (final s in loops[i]) s.a.clone()],
          curves: [for (final s in loops[i]) s.curve],
          holes: holes,
          holeCurves: holeCurves,
          outerArcs: outerArcs,
        ),
      );
    }
    return outers.isEmpty ? null : outers;
  }

  static bool _midOnBoundary(_ImprintSeg s, List<Vector3> poly, Vector3 n) =>
      _minDistToEdges(s.mid, poly, n) < 1e-4;

  static double _loopArea(List<Vector3> poly, Vector3 n) {
    if (poly.length < 3) return 0;
    var u = n.cross(Vector3(0, 0, 1));
    if (u.length < 1e-8) u = n.cross(Vector3(0, 1, 0));
    u.normalize();
    final v = n.cross(u);
    var area = 0.0;
    for (var i = 0; i < poly.length; i++) {
      final a = poly[i] - poly[0];
      final b = poly[(i + 1) % poly.length] - poly[0];
      area += a.dot(u) * b.dot(v) - a.dot(v) * b.dot(u);
    }
    return area * 0.5;
  }

  static List<Vector3> _segIntersect(_ImprintSeg a, _ImprintSeg b, Vector3 n) {
    final ac = a.curve;
    final bc = b.curve;
    if (ac is BrepCircle && bc is BrepLine) {
      return _hitsOnSeg(a, _circleSegmentHits(ac, b.a, b.b), n);
    }
    if (bc is BrepCircle && ac is BrepLine) {
      return _hitsOnSeg(b, _circleSegmentHits(bc, a.a, a.b), n);
    }
    if (ac is BrepLine && bc is BrepLine) {
      return _lineLineHit(a.a, a.b, b.a, b.b, n);
    }
    return const [];
  }

  static List<Vector3> _hitsOnSeg(_ImprintSeg s, List<Vector3> raw, Vector3 n) {
    final out = <Vector3>[];
    for (final p in raw) {
      if (_pointOnSeg(s, p, n)) out.add(p);
    }
    return out;
  }

  static bool _pointOnSeg(_ImprintSeg s, Vector3 p, Vector3 n) {
    if ((s.a - p).length < 1e-4 || (s.b - p).length < 1e-4) return true;
    if (s.curve is BrepCircle) {
      final c = s.curve as BrepCircle;
      if (((p - c.center).length - c.radius).abs() > 1e-3) return false;
      // Not closestParam: that clamps to t0/t1, which reports every point
      // on the full circle as lying on this arc and splits loops at points
      // outside the face.
      return _arcSpans(c, _circleAngle(c, p));
    }
    return _minDistToEdges(p, [s.a, s.b], n) < 1e-4;
  }

  static double _circleAngle(BrepCircle c, Vector3 p) {
    final d = p - c.center;
    final rad = d - c.normal * d.dot(c.normal);
    return math.atan2(rad.dot(c.yDir), rad.dot(c.xDir));
  }

  /// True when raw angle [t] lies within the arc's parameter range.
  static bool _arcSpans(BrepCircle c, double t) {
    const twoPi = math.pi * 2;
    final lo = math.min(c.t0, c.t1);
    final hi = math.max(c.t0, c.t1);
    if (hi - lo >= twoPi - 1e-6) return true;
    var u = t;
    while (u < lo - 1e-6) {
      u += twoPi;
    }
    while (u >= lo + twoPi - 1e-6) {
      u -= twoPi;
    }
    return u <= hi + 1e-4;
  }

  static List<Vector3> _lineLineHit(
    Vector3 a0,
    Vector3 a1,
    Vector3 b0,
    Vector3 b1,
    Vector3 n,
  ) {
    final da = a1 - a0;
    final db = b1 - b0;
    final denom = da.cross(db).dot(n);
    if (denom.abs() < 1e-12) return const [];
    final w = a0 - b0;
    final s = db.cross(w).dot(n) / denom;
    final t = da.cross(w).dot(n) / denom;
    if (s < -1e-6 || s > 1 + 1e-6 || t < -1e-6 || t > 1 + 1e-6) {
      return const [];
    }
    return [a0 + da * s.clamp(0.0, 1.0)];
  }

  static List<_ImprintSeg> _splitSegsAtHits(
    List<_ImprintSeg> segs,
    List<Vector3> hits,
    Vector3 n,
  ) {
    final out = <_ImprintSeg>[];
    for (final s in segs) {
      final on = <Vector3>[s.a.clone()];
      for (final h in hits) {
        if (_pointOnSeg(s, h, n)) {
          var dup = false;
          for (final p in on) {
            if ((p - h).length < 1e-4) {
              dup = true;
              break;
            }
          }
          if (!dup) on.add(h.clone());
        }
      }
      var endDup = false;
      for (final p in on) {
        if ((p - s.b).length < 1e-4) {
          endDup = true;
          break;
        }
      }
      if (!endDup) on.add(s.b.clone());
      if (on.length <= 2) {
        out.add(s);
        continue;
      }
      on.sort((p, q) => _segParam(s, p).compareTo(_segParam(s, q)));
      for (var i = 0; i < on.length - 1; i++) {
        if ((on[i] - on[i + 1]).length < 1e-8) continue;
        out.add(_sliceSeg(s, on[i], on[i + 1]));
      }
    }
    return out;
  }

  static double _segParam(_ImprintSeg s, Vector3 p) {
    if (s.curve is BrepCircle) {
      return (s.curve as BrepCircle).closestParam(p);
    }
    final d = s.b - s.a;
    final len2 = d.length2;
    if (len2 < 1e-20) return 0;
    return (p - s.a).dot(d) / len2;
  }

  static _ImprintSeg _sliceSeg(_ImprintSeg s, Vector3 a, Vector3 b) {
    final c = s.curve;
    if (c is BrepCircle) {
      var t0 = c.closestParam(a);
      var t1 = c.closestParam(b);
      if (c.t1 >= c.t0) {
        while (t1 < t0) {
          t1 += math.pi * 2;
        }
      } else {
        while (t0 < t1) {
          t0 += math.pi * 2;
        }
      }
      return _ImprintSeg(
        a,
        b,
        BrepCircle(
          center: c.center,
          normal: c.normal,
          radius: c.radius,
          xDir: c.xDir,
          t0: t0,
          t1: t1,
        ),
      );
    }
    return _ImprintSeg(a, b, BrepLine(a, b));
  }

  static List<List<_ImprintSeg>> _chainImprintSegs(List<_ImprintSeg> raw) {
    const tol = 1e-3;
    final unused = [...raw];
    final loops = <List<_ImprintSeg>>[];
    while (unused.isNotEmpty) {
      final loop = [unused.removeAt(0)];
      var guard = 0;
      while (guard++ < raw.length + 2) {
        final end = loop.last.b;
        var found = -1;
        var flip = false;
        for (var i = 0; i < unused.length; i++) {
          if ((unused[i].a - end).length < tol) {
            found = i;
            break;
          }
          if ((unused[i].b - end).length < tol) {
            found = i;
            flip = true;
            break;
          }
        }
        if (found < 0) break;
        var s = unused.removeAt(found);
        if (flip) s = _ImprintSeg(s.b.clone(), s.a.clone(), s.curve);
        loop.add(s);
      }
      if ((loop.last.b - loop.first.a).length > tol) return const [];
      if (loop.length >= 2 ||
          (loop.length == 1 &&
              loop.first.curve is BrepCircle &&
              (loop.first.curve as BrepCircle).isClosed)) {
        loops.add(loop);
      }
    }
    return loops;
  }

  static bool _isPrismAlong(BrepSolid s, Vector3 axis) {
    final ax = axis.normalized();
    var caps = 0;
    var walls = 0;
    for (final f in s.faces) {
      final surf = f.surface;
      if (surf is BrepPlane) {
        final an = surf.normal.dot(ax).abs();
        if (an > 0.999) {
          caps++;
        } else if (an < 1e-3) {
          walls++;
        } else {
          return false;
        }
      } else if (surf is BrepCylinder) {
        if (surf.axis.cross(ax).length > 1e-3) return false;
        walls++;
      } else {
        return false;
      }
    }
    // A cylinder primitive is a prism along its axis with a single wall.
    return caps >= 1 && walls >= 1;
  }

  /// Exact finite sections of planar-faced cutters, independent of sweep
  /// direction. Pair crossings by face membership, so concave faces and
  /// faces with holes can contribute more than one interval.
  static List<({Vector3 a, Vector3 b, BrepCurve curve})>? _planarSolidSection(
    BrepSolid solid,
    BrepPlane plane,
  ) {
    final result = <({Vector3 a, Vector3 b, BrepCurve curve})>[];
    const tol = 1e-8;
    for (final face in solid.faces) {
      final surface = face.surface;
      if (surface is! BrepPlane) return null;
      final direction = plane.normal.cross(surface.normal);
      if (direction.length < tol) {
        if (plane.height(surface.origin).abs() < tol) return null;
        continue;
      }
      direction.normalize();
      final hits = <Vector3>[];
      void add(Vector3 p) {
        if (!hits.any((q) => (p - q).length < tol)) hits.add(p);
      }

      for (final loop in face.loops) {
        for (final co in loop.coedges) {
          final curve = solid.edges[co.edge].curve;
          if (curve is! BrepLine) return null;
          final a = curve.a, b = curve.b;
          final ha = plane.height(a), hb = plane.height(b);
          if (ha.abs() < tol) add(a);
          if (hb.abs() < tol) add(b);
          if ((ha < -tol && hb > tol) || (hb < -tol && ha > tol)) {
            add(a + (b - a) * (ha / (ha - hb)));
          }
        }
      }
      hits.sort((a, b) => a.dot(direction).compareTo(b.dot(direction)));
      for (var i = 0; i + 1 < hits.length; i++) {
        final a = hits[i], b = hits[i + 1];
        final mid = (a + b) * .5;
        if ((BrepClassify.closestOnFace(solid, face, mid) - mid).length > tol) {
          continue;
        }
        result.add((a: a, b: b, curve: BrepLine(a, b)));
      }
    }
    return result;
  }

  static bool _prismMeetsPlane(BrepSolid s, BrepPlane plane) {
    if (s.vertices.isEmpty) return false;
    var lo = double.infinity;
    var hi = -double.infinity;
    for (final v in s.vertices) {
      final h = plane.height(v);
      lo = math.min(lo, h);
      hi = math.max(hi, h);
    }
    return lo <= 1e-4 && hi >= -1e-4;
  }

  static ({Vector3 a, Vector3 b, BrepCurve curve})? _wallSectionOnPlane(
    BrepSolid s,
    BrepFace face,
    BrepPlane plane,
    Vector3 axis,
  ) {
    final surf = face.surface;
    if (surf is BrepCylinder) {
      if (surf.axis.cross(axis).length > 1e-3) return null;
      final domain = BrepTessellate.cylinderFaceDomain(s, face);
      var u0 = domain.$1;
      var u1 = domain.$2;
      if (u1 <= u0 + 1e-9) u1 += math.pi * 2;
      final z = (plane.origin - surf.origin).dot(surf.axis);
      final center = surf.origin + surf.axis * z;
      if ((u1 - u0) >= math.pi * 2 - 1e-3) {
        final full = BrepCircle(
          center: center,
          normal: surf.axis,
          radius: surf.radius,
          xDir: surf.xDir,
        );
        final seam = full.evaluate(full.t0);
        return (a: seam, b: seam, curve: full);
      }
      final curve = BrepCircle(
        center: center,
        normal: surf.axis,
        radius: surf.radius,
        xDir: surf.xDir,
        t0: u0,
        t1: u1,
      );
      return (a: surf.evaluate(u0, z), b: surf.evaluate(u1, z), curve: curve);
    }
    if (surf is! BrepPlane) return null;
    if (surf.normal.dot(axis).abs() > 1e-3) return null;
    if (face.loops.isEmpty) return null;
    Vector3? a0;
    Vector3? b0;
    for (final c in face.loops.first.coedges) {
      final e = s.edges[c.edge];
      final a = s.vertices[c.reversed ? e.v1 : e.v0];
      final b = s.vertices[c.reversed ? e.v0 : e.v1];
      if ((b - a).dot(axis).abs() > 1e-3) continue;
      a0 = a;
      b0 = b;
      break;
    }
    if (a0 == null || b0 == null) return null;
    Vector3 onto(Vector3 p) => p - plane.normal * plane.height(p);
    final a = onto(a0);
    final b = onto(b0);
    return (a: a, b: b, curve: BrepLine(a, b));
  }

  static List<({Vector3 a, Vector3 b, BrepCurve curve})>? _largestClosedChain(
    List<({Vector3 a, Vector3 b, BrepCurve curve})> raw,
  ) {
    final all = _closedChains(raw);
    if (all.isEmpty) return null;
    if (all.length == 1) return all.first;
    var best = all.first;
    var bestArea = -1.0;
    for (final loop in all) {
      final samples = <Vector3>[];
      for (final s in loop) {
        samples.addAll(_sampleSection(s.curve));
      }
      if (samples.length < 3) continue;
      final n = (samples[1] - samples[0]).cross(samples[2] - samples[0]);
      if (n.length2 < 1e-20) continue;
      final area = _loopArea(samples, n.normalized()).abs();
      if (area > bestArea) {
        bestArea = area;
        best = loop;
      }
    }
    return best;
  }

  static List<List<({Vector3 a, Vector3 b, BrepCurve curve})>> _closedChains(
    List<({Vector3 a, Vector3 b, BrepCurve curve})> raw,
  ) {
    const tol = 1e-3;
    final unused = [...raw];
    final loops = <List<({Vector3 a, Vector3 b, BrepCurve curve})>>[];
    while (unused.isNotEmpty) {
      final out = [unused.removeAt(0)];
      var guard = 0;
      while (guard++ < raw.length + 2) {
        final end = out.last.b;
        var found = -1;
        var flip = false;
        for (var i = 0; i < unused.length; i++) {
          if ((unused[i].a - end).length < tol) {
            found = i;
            break;
          }
          if ((unused[i].b - end).length < tol) {
            found = i;
            flip = true;
            break;
          }
        }
        if (found < 0) break;
        var s = unused.removeAt(found);
        if (flip) s = (a: s.b, b: s.a, curve: s.curve);
        out.add(s);
      }
      if ((out.last.b - out.first.a).length > tol) continue;
      if (out.length >= 2 ||
          (out.length == 1 &&
              out.first.curve is BrepCircle &&
              (out.first.curve as BrepCircle).isClosed)) {
        loops.add(out);
      }
    }
    return loops;
  }

  static List<Vector3> _sampleSection(BrepCurve curve) {
    if (curve is BrepLine) {
      final mid = (curve.a + curve.b) * 0.5;
      return [curve.a.clone(), mid, curve.b.clone()];
    }
    if (curve is BrepCircle) {
      final n = curve.isClosed ? 16 : 8;
      return [
        for (var i = 0; i <= n; i++)
          curve.evaluate(curve.t0 + (curve.t1 - curve.t0) * (i / n)),
      ];
    }
    return [curve.evaluate(curve.t0), curve.evaluate(curve.t1)];
  }

  static bool _strictlyInside(
    List<Vector3> poly,
    List<Vector3> samples,
    Vector3 n,
  ) {
    for (final p in samples) {
      if (!_pointInPoly(poly, p, n)) return false;
      if (_minDistToEdges(p, poly, n) < 1e-4) return false;
    }
    return true;
  }

  static bool _pointInFaceMaterial(
    List<Vector3> poly,
    Vector3 p,
    Vector3 n, {
    required List<List<Vector3>> holes,
    required List<BrepCircle> circles,
  }) {
    if (!_pointInPoly(poly, p, n)) return false;
    if (_minDistToEdges(p, poly, n) < 1e-4) return true;
    for (final hole in holes) {
      if (_pointInPoly(hole, p, n) && _minDistToEdges(p, hole, n) > 1e-4) {
        return false;
      }
    }
    for (final c in circles) {
      final d = (p - c.center).length;
      if (d < c.radius - 1e-4) return false;
    }
    return true;
  }

  static bool _profileOverlapsFace(
    List<Vector3> poly,
    List<Vector3> clipPoly,
    List<_ImprintSeg> subject,
    List<_ImprintSeg> clip,
    Vector3 n,
    ({
      List<List<Vector3>> holes,
      List<List<BrepCurve>> curves,
      List<BrepCircle> circles,
    })
    prior,
  ) {
    bool onMaterial(Vector3 p) => _pointInFaceMaterial(
      poly,
      p,
      n,
      holes: prior.holes,
      circles: prior.circles,
    );
    for (final p in clipPoly) {
      if (onMaterial(p)) return true;
    }
    for (final s in clip) {
      if (onMaterial(s.mid)) return true;
    }
    for (final p in poly) {
      if (_pointInPoly(clipPoly, p, n) && onMaterial(p)) return true;
    }
    for (final a in subject) {
      for (final b in clip) {
        if (_segIntersect(a, b, n).isNotEmpty) return true;
      }
    }
    for (final hole in prior.holes) {
      if (hole.length < 2) continue;
      for (var i = 0; i < hole.length; i++) {
        final ha = hole[i];
        final hb = hole[(i + 1) % hole.length];
        final hs = _ImprintSeg(ha, hb, BrepLine(ha, hb));
        for (final c in clip) {
          if (_segIntersect(hs, c, n).isNotEmpty) return true;
        }
      }
    }
    return false;
  }

  static bool _profileOutspansPoly(
    List<Vector3> poly,
    List<Vector3> profile,
    Vector3 n,
  ) {
    for (final p in profile) {
      if (!_pointInPoly(poly, p, n) && _minDistToEdges(p, poly, n) > 1e-4) {
        return true;
      }
    }
    return false;
  }

  static bool _profileCoversPoly(
    List<Vector3> poly,
    List<Vector3> profile,
    Vector3 n,
  ) {
    if (profile.length < 3) return false;
    for (final p in poly) {
      if (!_pointInPoly(profile, p, n) &&
          _minDistToEdges(p, profile, n) > 1e-4) {
        return false;
      }
    }
    return true;
  }

  static List<Vector3> _planePoly(BrepSolid s, BrepFace f) {
    final pts = <Vector3>[];
    if (f.loops.isEmpty) return pts;
    for (final c in f.loops.first.coedges) {
      final e = s.edges[c.edge];
      if (e.curve is BrepLine) {
        pts.add(s.vertices[c.reversed ? e.v1 : e.v0]);
      } else {
        return BrepTessellate.sampleLoop(s, f.loops.first, 32);
      }
    }
    return pts;
  }

  static List<List<Vector3>> _splitPolyByPlane(
    List<Vector3> poly,
    BrepPlane plane,
  ) {
    var on = 0;
    var minH = double.infinity;
    var maxH = -double.infinity;
    for (final p in poly) {
      final h = plane.height(p);
      if (h.abs() < 1e-8) on++;
      minH = math.min(minH, h);
      maxH = math.max(maxH, h);
    }
    if (on == poly.length) return [poly];
    // Cap planes that only graze a boundary would emit a zero-thickness
    // sliver; merging that sliver later glues tangent wall-halves back.
    if (maxH <= 1e-8 || minH >= -1e-8) return [poly];
    final pos = _clip(poly, plane, true);
    final neg = _clip(poly, plane, false);
    final out = <List<Vector3>>[];
    if (pos.length >= 3) out.add(pos);
    if (neg.length >= 3) out.add(neg);
    return out.isEmpty ? [poly] : out;
  }

  static List<Vector3> _clip(
    List<Vector3> poly,
    BrepPlane plane,
    bool keepPos,
  ) {
    if (poly.isEmpty) return poly;
    final out = <Vector3>[];
    var prev = poly.last;
    var prevH = plane.height(prev);
    for (final cur in poly) {
      final h = plane.height(cur);
      final prevIn = keepPos ? prevH >= -1e-9 : prevH < -1e-9;
      final curIn = keepPos ? h >= -1e-9 : h < -1e-9;
      if (curIn) {
        if (!prevIn) {
          out.add(_lerpPlane(prev, cur, prevH, h));
        }
        out.add(cur.clone());
      } else if (prevIn) {
        out.add(_lerpPlane(prev, cur, prevH, h));
      }
      prev = cur;
      prevH = h;
    }
    return out;
  }

  static Vector3 _lerpPlane(Vector3 a, Vector3 b, double ha, double hb) {
    final denom = ha - hb;
    if (denom.abs() < 1e-18) return a.clone();
    final t = ha / denom;
    return a + (b - a) * t;
  }

  static double _radial(BrepCylinder cyl, Vector3 p) {
    final d = p - cyl.origin;
    return (d - cyl.axis * d.dot(cyl.axis)).length;
  }

  /// [s] extent along [cyl.axis], in the cylinder's V parameter.
  static (double, double)? _axisSpan(BrepSolid s, BrepCylinder cyl) {
    if (s.vertices.isEmpty) return null;
    var z0 = double.infinity;
    var z1 = -double.infinity;
    for (final v in s.vertices) {
      final z = (v - cyl.origin).dot(cyl.axis);
      z0 = math.min(z0, z);
      z1 = math.max(z1, z);
    }
    if (z1 <= z0) return null;
    return (z0, z1);
  }

  static bool _aabbsOverlap(
    List<Vector3> a,
    List<Vector3> b, {
    double pad = 1e-4,
  }) {
    if (a.isEmpty || b.isEmpty) return false;
    var amin = a.first.clone();
    var amax = a.first.clone();
    for (final p in a) {
      amin = Vector3(
        math.min(amin.x, p.x),
        math.min(amin.y, p.y),
        math.min(amin.z, p.z),
      );
      amax = Vector3(
        math.max(amax.x, p.x),
        math.max(amax.y, p.y),
        math.max(amax.z, p.z),
      );
    }
    var bmin = b.first.clone();
    var bmax = b.first.clone();
    for (final p in b) {
      bmin = Vector3(
        math.min(bmin.x, p.x),
        math.min(bmin.y, p.y),
        math.min(bmin.z, p.z),
      );
      bmax = Vector3(
        math.max(bmax.x, p.x),
        math.max(bmax.y, p.y),
        math.max(bmax.z, p.z),
      );
    }
    return amin.x <= bmax.x + pad &&
        amax.x >= bmin.x - pad &&
        amin.y <= bmax.y + pad &&
        amax.y >= bmin.y - pad &&
        amin.z <= bmax.z + pad &&
        amax.z >= bmin.z - pad;
  }

  /// True when two planar faces actually meet, not just their infinite planes.
  static bool _finitePlanesHit(
    List<Vector3> poly,
    BrepPlane selfPl,
    List<Vector3> otherPoly,
    BrepPlane otherPl,
  ) {
    if (poly.length < 3 || otherPoly.length < 3) return false;
    final an = selfPl.normal.dot(otherPl.normal).abs();
    if (an > 0.999) {
      return (otherPl.origin - selfPl.origin).dot(selfPl.normal).abs() < 1e-5;
    }
    bool edgeOn(List<Vector3> a, BrepPlane pl, List<Vector3> b, Vector3 bn) {
      for (var i = 0; i < a.length; i++) {
        final p = a[i];
        final q = a[(i + 1) % a.length];
        if (pl.height(p).abs() > 1e-4 || pl.height(q).abs() > 1e-4) continue;
        if ((p - q).length < 1e-8) continue;
        final mid = (p + q) * 0.5;
        if (_pointInPoly(b, mid, bn) || _minDistToEdges(mid, b, bn) <= 1e-4) {
          return true;
        }
        if (_pointInPoly(b, p, bn) || _pointInPoly(b, q, bn)) return true;
      }
      return false;
    }

    if (edgeOn(otherPoly, selfPl, poly, selfPl.normal)) return true;
    if (edgeOn(poly, otherPl, otherPoly, otherPl.normal)) return true;
    return _intersectionSegmentHitsBoth(poly, selfPl, otherPoly, otherPl);
  }

  /// True when the two finite faces share a positive-length intersection
  /// segment (not merely that their infinite planes cross).
  static bool _intersectionSegmentHitsBoth(
    List<Vector3> poly,
    BrepPlane selfPl,
    List<Vector3> otherPoly,
    BrepPlane otherPl,
  ) {
    final dir = selfPl.normal.cross(otherPl.normal);
    if (dir.length2 < 1e-20) return false;
    dir.normalize();
    final n1 = selfPl.normal;
    final n2 = otherPl.normal;
    final d1 = n1.dot(selfPl.origin);
    final d2 = n2.dot(otherPl.origin);
    final n12 = n1.cross(n2);
    if (n12.length2 < 1e-20) return false;
    final origin = (n2 * d1 - n1 * d2).cross(n12) / n12.length2;
    final a = _lineTsOnPoly(origin, dir, poly, selfPl.normal);
    final b = _lineTsOnPoly(origin, dir, otherPoly, otherPl.normal);
    if (a == null || b == null) return false;
    final lo = math.max(a.$1, b.$1);
    final hi = math.min(a.$2, b.$2);
    return hi - lo > 1e-6;
  }

  static (double, double)? _lineTsOnPoly(
    Vector3 origin,
    Vector3 dir,
    List<Vector3> poly,
    Vector3 n,
  ) {
    final ts = <double>[];
    for (var i = 0; i < poly.length; i++) {
      final a = poly[i];
      final b = poly[(i + 1) % poly.length];
      final ab = b - a;
      final denom = ab.cross(dir).dot(n);
      if (denom.abs() < 1e-12) continue;
      final w = a - origin;
      final s = dir.cross(w).dot(n) / denom;
      if (s < -1e-6 || s > 1 + 1e-6) continue;
      ts.add(ab.cross(w).dot(n) / denom);
    }
    if (ts.length < 2) return null;
    ts.sort();
    return (ts.first, ts.last);
  }

  /// True when [poly] on [pl] actually meets the cylinder wall, not just the
  /// infinite plane.
  static bool _planeFaceHitsCylinder(
    BrepPlane pl,
    List<Vector3> poly,
    BrepCylinder cyl,
    double z0,
    double z1,
  ) {
    if (poly.length < 3) return false;
    final an = pl.normal.dot(cyl.axis).abs();
    if (an > 0.999) {
      var minR = double.infinity;
      var maxR = 0.0;
      var zHit = false;
      for (final p in poly) {
        final z = (p - cyl.origin).dot(cyl.axis);
        if (z >= z0 - 1e-4 && z <= z1 + 1e-4) zHit = true;
        final r = _radial(cyl, p);
        minR = math.min(minR, r);
        maxR = math.max(maxR, r);
      }
      if (!zHit) return false;
      // Vertices all inside the circle: the face never meets the wall.
      if (maxR < cyl.radius - 1e-4) return false;
      final z = (pl.origin - cyl.origin).dot(cyl.axis);
      final center = cyl.origin + cyl.axis * z;
      // Through-hole: every vertex is outside the circle, but the disk
      // sits on this face and must clip the tool wall at this Z.
      if (_pointInPoly(poly, center, pl.normal)) return true;
      if (_minDistToEdges(center, poly, pl.normal) <= cyl.radius + 1e-4) {
        return true;
      }
      return minR <= cyl.radius + 1e-4;
    }
    if (an < 1e-3) {
      final curves = BrepIntersect.surfaces(pl, cyl);
      for (final curve in curves) {
        if (curve is! BrepLine) continue;
        if (_generatorHitsPoly(curve, poly, pl.normal, cyl, z0, z1)) {
          return true;
        }
      }
      return false;
    }
    return _polyMeetsCylinderWall(poly, cyl, z0, z1);
  }

  /// True when [poly] comes within a wall-thickness of the finite cylinder.
  static bool _polyMeetsCylinderWall(
    List<Vector3> poly,
    BrepCylinder cyl,
    double z0,
    double z1,
  ) {
    var minR = double.infinity;
    var maxR = 0.0;
    var zMin = double.infinity;
    var zMax = -double.infinity;
    void consider(Vector3 p) {
      final z = (p - cyl.origin).dot(cyl.axis);
      zMin = math.min(zMin, z);
      zMax = math.max(zMax, z);
      final r = _radial(cyl, p);
      minR = math.min(minR, r);
      maxR = math.max(maxR, r);
    }

    for (final p in poly) {
      consider(p);
    }
    for (var i = 0; i < poly.length; i++) {
      final a = poly[i];
      final b = poly[(i + 1) % poly.length];
      var best = a;
      var bestD = _radial(cyl, a);
      const samples = 8;
      for (var s = 1; s <= samples; s++) {
        final p = a + (b - a) * (s / samples);
        final d = _radial(cyl, p);
        if (d < bestD) {
          bestD = d;
          best = p;
        }
      }
      consider(best);
    }
    if (zMax < z0 - 1e-4 || zMin > z1 + 1e-4) return false;
    if (maxR < cyl.radius - 1e-4) return false;
    if (minR > cyl.radius + 1e-4) return false;
    return true;
  }

  static bool _generatorHitsPoly(
    BrepLine line,
    List<Vector3> poly,
    Vector3 n,
    BrepCylinder cyl,
    double z0,
    double z1,
  ) {
    final midZ = 0.5 * (z0 + z1);
    final axisPt = cyl.origin + cyl.axis * midZ;
    final q = line.closestPoint(axisPt);
    final z = (q - cyl.origin).dot(cyl.axis);
    if (z < z0 - 1e-4 || z > z1 + 1e-4) return false;
    return _pointInPoly(poly, _projectToPolyPlane(q, poly, n), n) ||
        _pointInPoly(poly, q, n);
  }

  static Vector3 _projectToPolyPlane(Vector3 q, List<Vector3> poly, Vector3 n) {
    final o = poly.first;
    return q - n * (q - o).dot(n);
  }

  static List<List<Vector3>> _splitPolyByCylinderGenerators(
    List<Vector3> poly,
    BrepCylinder cyl,
    BrepPlane plane,
  ) {
    var parts = [poly];
    final curves = BrepIntersect.surfaces(plane, cyl);
    for (final curve in curves) {
      if (curve is! BrepLine) continue;
      final a = curve.evaluate(0.5);
      var dir = curve.b - curve.a;
      if (dir.length2 < 1e-20) continue;
      dir.normalize();
      var sn = plane.normal.cross(dir);
      if (sn.length2 < 1e-20) continue;
      sn.normalize();
      final next = <List<Vector3>>[];
      final splitter = BrepPlane(origin: a, normal: sn);
      for (final p in parts) {
        final zRef = (_centroid(p) - cyl.origin).dot(cyl.axis);
        if (p.length < 3 ||
            !_generatorHitsPoly(curve, p, plane.normal, cyl, zRef, zRef)) {
          next.add(p);
          continue;
        }
        next.addAll(_splitPolyByPlane(p, splitter));
      }
      parts = next;
    }
    return parts;
  }

  static Vector3 _probeOnOuter(List<Vector3> poly) {
    if (poly.length < 2) return _centroid(poly);
    final mid = (poly[0] + poly[1]) * 0.5;
    if (poly.length < 3) return mid;
    final inward = _centroid(poly) - mid;
    if (inward.length < 1e-12) return mid;
    return mid + inward.normalized() * 1e-3;
  }

  static Vector3 _arcProbe(List<Vector3> poly, BrepCircle arc) {
    final mid = arc.evaluate(0.5 * (arc.t0 + arc.t1));
    if (poly.length == 2) {
      return (poly[0] + poly[1] + mid) / 3.0;
    }
    return _centroid(poly);
  }

  static ({
    List<({List<Vector3> poly, BrepCircle? arc12})> outers,
    List<Vector3>? hole,
    BrepCircle? holeCircle,
  })
  _cutPlaneByCylinder(
    List<Vector3> poly,
    BrepCylinder cyl,
    BrepPlane plane, {
    required int op,
    required bool fromA,
    double u0 = 0,
    double u1 = math.pi * 2,
    BrepSolid? other,
  }) {
    if (poly.length < 3) {
      return (
        outers: [(poly: poly, arc12: null)],
        hole: null,
        holeCircle: null,
      );
    }
    final z = (plane.origin - cyl.origin).dot(cyl.axis);
    final center = cyl.origin + cyl.axis * z;
    var ua = u0;
    var ub = u1;
    if (ub <= ua + 1e-9) ub += math.pi * 2;
    if ((ub - ua) < math.pi * 2 - 1e-3) {
      return _cutPlaneByCylinderArc(poly, cyl, plane, ua, ub, z, center);
    }
    const n = 32;
    final disk = [
      for (var i = 0; i < n; i++) cyl.evaluate(i / n * math.pi * 2, z),
    ];
    // Sequential half-space clips cannot represent "polygon minus disk".
    // Faces that miss the disk entirely (hollow opening lips, a nearby
    // coplanar patch) would be wiped by a tangent facing the other way.
    final centerIn = _pointInPoly(poly, center, plane.normal);
    if (!centerIn &&
        _minDistToEdges(center, poly, plane.normal) > cyl.radius + 1e-4) {
      return (
        outers: [(poly: poly, arc12: null)],
        hole: null,
        holeCircle: null,
      );
    }
    final vertsInside = poly.every((p) => _radial(cyl, p) <= cyl.radius + 1e-6);
    final containsDisk =
        centerIn && disk.every((p) => _pointInPoly(poly, p, plane.normal));
    final punch = op == subtractOp && fromA;
    if (vertsInside) {
      if (punch) {
        return (outers: const [], hole: null, holeCircle: null);
      }
      return (
        outers: [(poly: poly, arc12: null)],
        hole: null,
        holeCircle: null,
      );
    }
    if (containsDisk) {
      if (op == unionOp) {
        final circle = BrepCircle(
          center: center,
          normal: cyl.axis.dot(plane.normal) >= 0 ? cyl.axis : -cyl.axis,
          radius: cyl.radius,
          xDir: cyl.xDir,
        );
        return (
          outers: [(poly: poly, arc12: null)],
          // The analytic circle is the hole boundary. Adding its sampled
          // polygon as a second hole duplicates the cut and prevents the
          // untouched outer trim from being retained.
          hole: null,
          holeCircle: circle,
        );
      }
      if (punch) {
        final minEdge = _minDistToEdges(center, poly, plane.normal);
        // Inscribed / tangent: the circle pinches the outer boundary into
        // separate corner pieces instead of one face with a hole.
        if (minEdge <= cyl.radius + 1e-4) {
          return (
            outers: _cornerCaps(poly, cyl, plane, center),
            hole: null,
            holeCircle: null,
          );
        }
        final circle = BrepCircle(
          center: center,
          normal: cyl.axis.dot(plane.normal) >= 0 ? cyl.axis : -cyl.axis,
          radius: cyl.radius,
          xDir: cyl.xDir,
        );
        return (
          outers: [(poly: poly, arc12: null)],
          hole: null,
          holeCircle: circle,
        );
      }
      // Slice / intersect against a holed body: the disk is a void, not
      // the cutting-face material.
      final diskIsVoid =
          other != null && BrepClassify.side(other, center) == BrepSide.outside;
      if (diskIsVoid) {
        final circle = BrepCircle(
          center: center,
          normal: cyl.axis.dot(plane.normal) >= 0 ? cyl.axis : -cyl.axis,
          radius: cyl.radius,
          xDir: cyl.xDir,
        );
        return (
          outers: [(poly: poly, arc12: null)],
          hole: null,
          holeCircle: circle,
        );
      }
      return (
        outers: [(poly: disk, arc12: null)],
        hole: null,
        holeCircle: null,
      );
    }
    var cur = poly;
    for (var i = 0; i < disk.length; i++) {
      final a = disk[i];
      final b = disk[(i + 1) % disk.length];
      var edgeN = (b - a).cross(plane.normal);
      if (edgeN.length2 < 1e-20) continue;
      edgeN.normalize();
      if (edgeN.dot(center - a) < 0) edgeN = -edgeN;
      if (punch) edgeN = -edgeN;
      cur = _clip(cur, BrepPlane(origin: a, normal: edgeN), true);
      if (cur.length < 3) {
        return (
          outers: [(poly: cur, arc12: null)],
          hole: null,
          holeCircle: null,
        );
      }
    }
    return (outers: [(poly: cur, arc12: null)], hole: null, holeCircle: null);
  }

  /// Split [poly] along a cylinder *patch* arc (not a full disk).
  ///
  /// A half-cylinder D-notch must not punch a round hole through the leftover
  /// lens. Only the face whose interior contains the arc midpoint is split.
  static ({
    List<({List<Vector3> poly, BrepCircle? arc12})> outers,
    List<Vector3>? hole,
    BrepCircle? holeCircle,
  })
  _cutPlaneByCylinderArc(
    List<Vector3> poly,
    BrepCylinder cyl,
    BrepPlane plane,
    double ua,
    double ub,
    double z,
    Vector3 center,
  ) {
    final unchanged = (
      outers: [(poly: poly, arc12: null)],
      hole: null as List<Vector3>?,
      holeCircle: null as BrepCircle?,
    );
    final mid = cyl.evaluate(0.5 * (ua + ub), z);
    if (!_pointInPoly(poly, mid, plane.normal)) return unchanged;
    final p0 = cyl.evaluate(ua, z);
    final p1 = cyl.evaluate(ub, z);
    if (_minDistToEdges(p0, poly, plane.normal) > 1e-3 ||
        _minDistToEdges(p1, poly, plane.normal) > 1e-3) {
      return unchanged;
    }
    final split = _splitPolyAlongChord(poly, p0, p1, plane.normal);
    if (split == null) return unchanged;
    final arc = BrepCircle(
      center: center,
      // Same frame as the cylinder wall so cap arcs sew to the generators.
      normal: cyl.axis,
      radius: cyl.radius,
      xDir: cyl.xDir,
      t0: ua,
      t1: ub,
    );
    return (
      outers: [(poly: split.$1, arc12: arc), (poly: split.$2, arc12: arc)],
      hole: null,
      holeCircle: null,
    );
  }

  static (List<Vector3>, List<Vector3>)? _splitPolyAlongChord(
    List<Vector3> poly,
    Vector3 a,
    Vector3 b,
    Vector3 n,
  ) {
    final withA = _insertOnBoundary(poly, a, n);
    if (withA == null) return null;
    final withBoth = _insertOnBoundary(withA, b, n);
    if (withBoth == null) return null;
    final i0 = _nearestVert(withBoth, a);
    final i1 = _nearestVert(withBoth, b);
    if (i0 == i1) return null;

    List<Vector3> walk(int from, int to) {
      final out = <Vector3>[];
      var i = from;
      for (var k = 0; k <= withBoth.length; k++) {
        out.add(withBoth[i].clone());
        if (i == to) break;
        i = (i + 1) % withBoth.length;
      }
      return out;
    }

    var pathA = walk(i0, i1);
    var pathB = walk(i1, i0);
    if (pathA.length < 2 || pathB.length < 2) return null;
    return (pathA, pathB);
  }

  static List<Vector3>? _insertOnBoundary(
    List<Vector3> poly,
    Vector3 p,
    Vector3 n,
  ) {
    for (final v in poly) {
      if ((v - p).length < 1e-4) {
        return [for (final q in poly) q.clone()];
      }
    }
    var bestI = -1;
    var bestD = 1e-3;
    for (var i = 0; i < poly.length; i++) {
      final a = poly[i];
      final b = poly[(i + 1) % poly.length];
      final d = b - a;
      final len2 = d.length2;
      if (len2 < 1e-20) continue;
      final t = ((p - a).dot(d) / len2).clamp(0.0, 1.0);
      final q = a + d * t;
      final dist = (p - q).length;
      if (dist < bestD) {
        bestD = dist;
        bestI = i;
      }
    }
    if (bestI < 0) return null;
    final out = <Vector3>[for (final q in poly) q.clone()];
    out.insert(bestI + 1, p.clone());
    return out;
  }

  static int _nearestVert(List<Vector3> poly, Vector3 p) {
    var best = 0;
    var bestD = double.infinity;
    for (var i = 0; i < poly.length; i++) {
      final d = (poly[i] - p).length;
      if (d < bestD) {
        bestD = d;
        best = i;
      }
    }
    return best;
  }

  static double _minDistToEdges(Vector3 center, List<Vector3> poly, Vector3 n) {
    var best = double.infinity;
    for (var i = 0; i < poly.length; i++) {
      final a = poly[i];
      final b = poly[(i + 1) % poly.length];
      final d = b - a;
      final len2 = d.length2;
      final t = len2 < 1e-20
          ? 0.0
          : ((center - a).dot(d) / len2).clamp(0.0, 1.0);
      final q = a + d * t;
      final dist = (center - q).length;
      if (dist < best) best = dist;
    }
    return best;
  }

  static List<({List<Vector3> poly, BrepCircle? arc12})> _cornerCaps(
    List<Vector3> poly,
    BrepCylinder cyl,
    BrepPlane plane,
    Vector3 center,
  ) {
    final r = cyl.radius;
    final n = poly.length;
    final out = <({List<Vector3> poly, BrepCircle? arc12})>[];
    Vector3? hitToward(Vector3 from, Vector3 to) {
      final d = to - from;
      final len2 = d.length2;
      if (len2 < 1e-20) return null;
      final f = from - center;
      final a = d.dot(d);
      final b = 2 * f.dot(d);
      final c = f.dot(f) - r * r;
      var disc = b * b - 4 * a * c;
      if (disc.abs() <= 1e-14 * (b * b + (4 * a * c).abs())) disc = 0;
      if (disc < 0) {
        if (disc > -1e-8) {
          disc = 0;
        } else {
          return null;
        }
      }
      final s = math.sqrt(disc);
      final t0 = (-b - s) / (2 * a);
      final t1 = (-b + s) / (2 * a);
      double? best;
      for (final t in [t0, t1]) {
        if (t >= -1e-6 && t <= 1 + 1e-6) {
          if (best == null || t < best) best = t;
        }
      }
      if (best == null) return null;
      return from + d * best.clamp(0.0, 1.0);
    }

    for (var i = 0; i < n; i++) {
      final v = poly[i];
      if (_radial(cyl, v) <= r + 1e-5) continue;
      final prev = poly[(i - 1 + n) % n];
      final next = poly[(i + 1) % n];
      final hitNext = hitToward(v, next);
      final hitPrev = hitToward(v, prev);
      if (hitNext == null || hitPrev == null) continue;
      if ((hitNext - v).length < 1e-8 || (hitPrev - v).length < 1e-8) {
        continue;
      }
      var tA = cyl.parameters(hitNext)?.$1 ?? 0.0;
      var tB = cyl.parameters(hitPrev)?.$1 ?? 0.0;
      while (tB < tA) {
        tB += math.pi * 2;
      }
      if (tB - tA > math.pi) {
        tB -= math.pi * 2;
      }
      final arc = BrepCircle(
        center: center,
        // Match the cylinder frame so cap arcs sew to wall arcs.
        normal: cyl.axis,
        radius: r,
        xDir: cyl.xDir,
        t0: tA,
        t1: tB,
      );
      out.add((poly: [v, hitNext, hitPrev], arc12: arc));
    }
    return out;
  }

  static bool _pointInPoly(List<Vector3> poly, Vector3 p, Vector3 n) {
    if (poly.length < 3 || n.length2 < 1e-20) return false;
    final nn = n.normalized();
    if (_minDistToEdges(p, poly, nn) < 1e-6) return true;
    var xAxis = nn.cross(Vector3(0, 0, 1));
    if (xAxis.length2 < 1e-12) xAxis = nn.cross(Vector3(0, 1, 0));
    xAxis.normalize();
    final yAxis = nn.cross(xAxis);
    final o = poly.first;
    double xd(Vector3 q) => (q - o).dot(xAxis);
    double yd(Vector3 q) => (q - o).dot(yAxis);
    final px = xd(p);
    final py = yd(p);
    var inside = false;
    for (var i = 0, j = poly.length - 1; i < poly.length; j = i++) {
      final yi = yd(poly[i]);
      final yj = yd(poly[j]);
      if ((yi > py) == (yj > py)) continue;
      final xi = xd(poly[i]);
      final xj = xd(poly[j]);
      final xHit = (xj - xi) * (py - yi) / (yj - yi) + xi;
      if (px < xHit) inside = !inside;
    }
    return inside;
  }

  static Vector3 _centroid(List<Vector3> poly) {
    var c = Vector3.zero();
    for (final p in poly) {
      c += p;
    }
    return c / poly.length.toDouble();
  }

  static BrepCurve? _arcEdgeIfMatch(BrepCircle? arc, Vector3 a, Vector3 b) {
    if (arc == null) return null;
    final a0 = arc.evaluate(arc.t0);
    final a1 = arc.evaluate(arc.t1);
    const tol = 1e-3;
    if ((a - a0).length < tol && (b - a1).length < tol) return arc;
    if ((a - a1).length < tol && (b - a0).length < tol) return arc;
    return null;
  }

  static BrepCurve? _arcEdgeIfMatchAny(
    List<BrepCircle> arcs,
    Vector3 a,
    Vector3 b,
  ) {
    for (final arc in arcs) {
      final m = _arcEdgeIfMatch(arc, a, b);
      if (m != null) return m;
    }
    return null;
  }

  /// Two 1-use edges that share endpoints are the same curve assembled
  /// twice (line vs arc, or two arcs with different frames). Keep the
  /// analytic one so the solid sews.
  static BrepSolid _pairOpenTwins(BrepSolid s) {
    final uses = List<int>.filled(s.edges.length, 0);
    for (final f in s.faces) {
      for (final loop in f.loops) {
        for (final c in loop.coedges) {
          if (c.edge >= 0 && c.edge < uses.length) uses[c.edge]++;
        }
      }
    }
    final open = <int>[
      for (var i = 0; i < uses.length; i++)
        if (uses[i] == 1) i,
    ];
    if (open.length < 2) return s;
    final map = <int, (int, bool)>{};
    for (var i = 0; i < open.length; i++) {
      final a = open[i];
      if (map.containsKey(a)) continue;
      final ea = s.edges[a];
      final loA = math.min(ea.v0, ea.v1);
      final hiA = math.max(ea.v0, ea.v1);
      for (var j = i + 1; j < open.length; j++) {
        final b = open[j];
        if (map.containsKey(b)) continue;
        final eb = s.edges[b];
        final loB = math.min(eb.v0, eb.v1);
        final hiB = math.max(eb.v0, eb.v1);
        if (loA != loB || hiA != hiB) continue;
        final ca = ea.curve;
        final cb = eb.curve;
        // Equal endpoints do not identify a curve: complementary semicircles
        // share endpoints but bound opposite sides of a hole.
        bool follows(BrepCurve x, BrepCurve y) => [0.25, 0.5, 0.75].every((t) {
          final p = x.evaluate(x.t0 + (x.t1 - x.t0) * t);
          return (p - y.closestPoint(p)).length < 1e-7;
        });
        if (!follows(ca, cb) || !follows(cb, ca)) continue;
        final keep = ca is BrepCircle ? a : (cb is BrepCircle ? b : a);
        final drop = keep == a ? b : a;
        map[drop] = (
          keep,
          s.edges[drop].v0 == s.edges[drop].v1
              ? _oppositeCurveTravel(s.edges[drop].curve, s.edges[keep].curve)
              : s.edges[drop].v0 != s.edges[keep].v0,
        );
        break;
      }
    }
    if (map.isEmpty) return s;
    final faces = <BrepFace>[
      for (final f in s.faces)
        BrepFace(
          surface: f.surface,
          reversed: f.reversed,
          loops: [
            for (final loop in f.loops)
              BrepLoop([
                for (final c in loop.coedges)
                  BrepCoedge(
                    map[c.edge]?.$1 ?? c.edge,
                    reversed: c.reversed != (map[c.edge]?.$2 ?? false),
                  ),
              ]),
          ],
        ),
    ];
    return BrepSolid(vertices: s.vertices, edges: s.edges, faces: faces);
  }

  /// Orientation of the same geometric edge in two parameterizations. The
  /// sign of the parameter interval matters too (a circle can run 2*pi to 0).
  static bool _oppositeCurveTravel(BrepCurve a, BrepCurve b) {
    final t = a.t0 + (a.t1 - a.t0) * .37;
    final u = b.closestParam(a.evaluate(t));
    return a.tangent(t).dot(b.tangent(u)) * (a.t1 - a.t0) * (b.t1 - b.t0) < 0;
  }

  static BrepSolid _assemble(List<_Piece> pieces) {
    final b = BrepBuilder();
    int weld(Vector3 p) {
      for (var i = 0; i < b.vertices.length; i++) {
        if ((b.vertices[i] - p).length < 1e-6) return i;
      }
      return b.vertex(p);
    }

    final edgeOf = <(int, int), List<int>>{};
    (int, bool) sew(int i, int j, BrepCurve curve) {
      final key = (math.min(i, j), math.max(i, j));
      final candidates = edgeOf.putIfAbsent(key, () => []);
      for (final existing in candidates) {
        final e = b.edges[existing], stored = e.curve;
        // Use geometric agreement, not rounded midpoint strings: equivalent
        // rims can fall on opposite sides of a rounding boundary after a move.
        // Closed rims also require travel direction; their endpoints coincide.
        final reversed = i == j
            ? _oppositeCurveTravel(curve, stored)
            : e.v0 != i;
        bool follows(BrepCurve a, BrepCurve b) => [.25, .5, .75].every((f) {
          final p = a.evaluate(a.t0 + (a.t1 - a.t0) * f);
          return p.distanceTo(b.closestPoint(p)) < 1e-7;
        });
        if (!follows(curve, stored) || !follows(stored, curve)) continue;
        return (existing, reversed);
      }
      final id = b.edge(i, j, curve);
      candidates.add(id);
      return (id, false);
    }

    // Two intersecting round holes have one boundary, not two full circles.
    // Keep the portions of each circle outside the other disk analytically.
    void mergeRoundHoles(List<BrepLoop> loops) {
      final ids = <int>[];
      for (var i = 1; i < loops.length; i++) {
        final co = loops[i].coedges;
        if (co.length == 1 &&
            b.edges[co.single.edge].curve is BrepCircle &&
            (b.edges[co.single.edge].curve as BrepCircle).isClosed) {
          ids.add(i);
        }
      }
      if (ids.length != 2) return;
      final coA = loops[ids[0]].coedges.single;
      final coB = loops[ids[1]].coedges.single;
      final a = b.edges[coA.edge].curve as BrepCircle;
      final c = b.edges[coB.edge].curve as BrepCircle;
      final delta = c.center - a.center;
      final d = delta.length;
      if (a.normal.cross(c.normal).length > 1e-8 ||
          delta.dot(a.normal).abs() > 1e-6 ||
          d >= a.radius + c.radius - 1e-6 ||
          d <= (a.radius - c.radius).abs() + 1e-6) {
        return;
      }
      final x = (d * d + a.radius * a.radius - c.radius * c.radius) / (2 * d);
      final h = math.sqrt(math.max(0, a.radius * a.radius - x * x));
      final mid = a.center + delta * (x / d);
      final side = a.normal.cross(delta / d) * h;
      final p = mid + side, q = mid - side;
      BrepCoedge outsideArc(
        BrepCircle circle,
        BrepCircle other,
        BrepCoedge co,
      ) {
        double angle(Vector3 v) {
          final dv = v - circle.center;
          return math.atan2(dv.dot(circle.yDir), dv.dot(circle.xDir));
        }

        var start = angle(p), end = angle(q);
        while (end <= start) {
          end += 2 * math.pi;
        }
        if ((circle.evaluate((start + end) / 2) - other.center).length <
            other.radius) {
          final old = start;
          start = end;
          end = old + 2 * math.pi;
        }
        final arc = BrepCircle(
          center: circle.center,
          normal: circle.normal,
          radius: circle.radius,
          xDir: circle.xDir,
          t0: start,
          t1: end,
        );
        final edge = sew(
          weld(arc.evaluate(start)),
          weld(arc.evaluate(end)),
          arc,
        );
        return BrepCoedge(edge.$1, reversed: co.reversed != edge.$2);
      }

      final arcs = [outsideArc(a, c, coA), outsideArc(c, a, coB)];
      loops[ids[0]] = BrepLoop(arcs);
      loops.removeAt(ids[1]);
    }

    for (final piece in pieces) {
      if (piece.sourceFace != null && piece.sourceSolid != null) {
        final src = piece.sourceSolid!;
        final face = piece.sourceFace!;
        final loops = <BrepLoop>[];
        for (final loop in face.loops) {
          final es = <BrepCoedge>[];
          for (final c in loop.coedges) {
            final e = src.edges[c.edge];
            final sewn = sew(
              weld(src.vertices[e.v0]),
              weld(src.vertices[e.v1]),
              e.curve,
            );
            es.add(BrepCoedge(sewn.$1, reversed: c.reversed != sewn.$2));
          }
          loops.add(BrepLoop(es));
        }
        for (final circle in piece.holeCircles) {
          final seam = weld(circle.evaluate(0));
          final e = sew(seam, seam, circle);
          loops.add(BrepLoop([BrepCoedge(e.$1, reversed: !e.$2)]));
        }
        for (final curve in piece.innerClosed) {
          final seam = weld(curve.evaluate(curve.t0));
          final e = sew(seam, seam, curve);
          loops.add(BrepLoop([BrepCoedge(e.$1, reversed: !e.$2)]));
        }
        if (face.surface is BrepPlane) mergeRoundHoles(loops);
        b.face(face.surface, loops, reversed: piece.reversed);
        continue;
      }
      if (piece.surface is BrepCone &&
          piece.poly.length < 2 &&
          piece.outerCurves.isNotEmpty) {
        final loops = <BrepLoop>[];
        for (final curve in piece.outerCurves) {
          final seam = weld(curve.evaluate(curve.t0));
          final e = sew(seam, seam, curve);
          loops.add(BrepLoop([BrepCoedge(e.$1, reversed: !e.$2)]));
        }
        for (final circle in piece.holeCircles) {
          final seam = weld(circle.evaluate(0));
          final e = sew(seam, seam, circle);
          loops.add(BrepLoop([BrepCoedge(e.$1, reversed: !e.$2)]));
        }
        b.face(piece.surface, loops, reversed: piece.reversed);
        continue;
      }
      if (piece.surface is BrepSphere) {
        if (piece.holeCircles.isNotEmpty || piece.holeCurves.isNotEmpty) {
          final loops = <BrepLoop>[];
          for (final circle in piece.holeCircles) {
            final seam = weld(circle.evaluate(0));
            final e = sew(seam, seam, circle);
            var rev = !e.$2;
            final stored = b.edges[e.$1].curve;
            if (piece.poly.isNotEmpty && stored is BrepCircle) {
              final sph = piece.surface as BrepSphere;
              final want = piece.poly.first;
              final pole = sph.circularCapPole(stored, coedgeReversed: rev);
              final flipped = sph.circularCapPole(stored, coedgeReversed: !rev);
              if ((flipped - want).length < (pole - want).length) rev = !rev;
            }
            loops.add(BrepLoop([BrepCoedge(e.$1, reversed: rev)]));
          }
          for (final chain in piece.holeCurves) {
            if (chain.isEmpty) continue;
            if (chain.length == 1 &&
                chain.first is BrepCircle &&
                (chain.first as BrepCircle).isClosed) {
              final circle = chain.first as BrepCircle;
              final seam = weld(circle.evaluate(0));
              final e = sew(seam, seam, circle);
              loops.add(BrepLoop([BrepCoedge(e.$1, reversed: !e.$2)]));
              continue;
            }
            final es = <BrepCoedge>[];
            for (final curve in chain) {
              final ia = weld(curve.evaluate(curve.t0));
              final ib = weld(curve.evaluate(curve.t1));
              final e = sew(ia, ib, curve);
              es.add(BrepCoedge(e.$1, reversed: e.$2));
            }
            if (es.isNotEmpty) loops.add(BrepLoop(es));
          }
          b.face(piece.surface, loops, reversed: piece.reversed);
          continue;
        }
        if (piece.poly.isEmpty) {
          b.face(piece.surface, const [], reversed: piece.reversed);
          continue;
        }
      }
      if (piece.poly.length < 2) continue;
      if (piece.poly.length < 3 &&
          piece.arc12 == null &&
          piece.outerArcs.isEmpty) {
        continue;
      }
      if (piece.uv != null &&
          (piece.surface is BrepCylinder ||
              piece.surface is BrepTorus ||
              piece.surface is BrepCone)) {
        final uv = piece.uv!;
        Vector3 at(double u, double v) {
          final s = piece.surface;
          if (s is BrepTorus) return s.evaluate(u, v);
          if (s is BrepCone) return s.evaluate(u, v);
          return (s as BrepCylinder).evaluate(u, v);
        }

        BrepCurve rim(double ua, double ub, double v) {
          final s = piece.surface;
          if (s is BrepTorus) {
            final c = s.center + s.axis * (s.minorRadius * math.sin(v));
            final rad = (s.majorRadius + s.minorRadius * math.cos(v)).abs();
            return BrepCircle(
              center: c,
              normal: s.axis,
              radius: math.max(rad, 1e-8),
              xDir: s.xDir,
              t0: ua,
              t1: ub,
            );
          }
          if (s is BrepCone) {
            return BrepCircle(
              center: s.origin + s.axis * v,
              normal: s.axis,
              radius: math.max(s.radiusAt(v), 1e-8),
              xDir: s.xDir,
              t0: ua,
              t1: ub,
            );
          }
          final cyl = s as BrepCylinder;
          return BrepCircle(
            center: cyl.origin + cyl.axis * v,
            normal: cyl.axis,
            radius: cyl.radius,
            xDir: cyl.xDir,
            t0: ua,
            t1: ub,
          );
        }

        final p00 = weld(at(uv.$1, uv.$3));
        final p10 = weld(at(uv.$2, uv.$3));
        final p11 = weld(at(uv.$2, uv.$4));
        final p01 = weld(at(uv.$1, uv.$4));
        BrepCurve wall(int ia, int ib, double u, double va, double vb) {
          final s = piece.surface;
          if (s is BrepTorus) {
            final cr = s.xDir * math.cos(u) + s.yDir * math.sin(u);
            var n = cr.cross(s.axis);
            if (n.length < 1e-12) n = s.xDir.clone();
            n.normalize();
            return BrepCircle(
              center: s.center + cr * s.majorRadius,
              normal: n,
              radius: s.minorRadius,
              xDir: cr,
              t0: va,
              t1: vb,
            );
          }
          return BrepLine(b.vertices[ia], b.vertices[ib]);
        }

        final e0 = sew(p00, p10, rim(uv.$1, uv.$2, uv.$3));
        final e1 = sew(p10, p11, wall(p10, p11, uv.$2, uv.$3, uv.$4));
        final e2 = sew(p11, p01, rim(uv.$2, uv.$1, uv.$4));
        final e3 = sew(p01, p00, wall(p01, p00, uv.$1, uv.$4, uv.$3));
        b.face(piece.surface, [
          BrepLoop([
            BrepCoedge(e0.$1, reversed: e0.$2),
            BrepCoedge(e1.$1, reversed: e1.$2),
            BrepCoedge(e2.$1, reversed: e2.$2),
            BrepCoedge(e3.$1, reversed: e3.$2),
          ]),
        ], reversed: piece.reversed);
        continue;
      }
      if (piece.surface is! BrepPlane) {
        // Project polygon onto a fitted plane.
      }
      var poly = piece.poly;
      Vector3? sphereKeepPole;
      if (piece.surface is BrepSphere &&
          piece.outerArcs.isNotEmpty &&
          poly.length >= 3) {
        final sph = piece.surface as BrepSphere;
        final maybe = poly.first;
        final onLoop = poly.skip(1).any((p) => (p - maybe).length < 1e-5);
        if (!onLoop && (maybe - sph.center).length > sph.radius * 0.4) {
          sphereKeepPole = maybe;
          poly = poly.sublist(1);
        }
      }
      if (poly.length < 2) continue;
      final ids = [for (final p in poly) weld(p)];
      var eids = <(int, bool)>[];
      for (var i = 0; i < ids.length; i++) {
        final j = (i + 1) % ids.length;
        final matchArcs = [
          if (piece.arc12 != null) piece.arc12!,
          ...piece.outerArcs,
        ];
        final BrepCurve curve;
        if (i < piece.outerCurves.length) {
          curve = piece.outerCurves[i];
        } else if (ids.length == 2 && matchArcs.length == 1) {
          curve = i == 0
              ? matchArcs.first
              : BrepLine(b.vertices[ids[i]], b.vertices[ids[j]]);
        } else {
          curve =
              _arcEdgeIfMatchAny(
                matchArcs,
                b.vertices[ids[i]],
                b.vertices[ids[j]],
              ) ??
              BrepLine(b.vertices[ids[i]], b.vertices[ids[j]]);
        }
        eids.add(sew(ids[i], ids[j], curve));
      }
      if (sphereKeepPole != null && eids.isNotEmpty) {
        final sph = piece.surface as BrepSphere;
        final small = sph.closestPoint(_centroid(poly));
        final large = sph.center * 2.0 - small;
        final wantSmall =
            (sphereKeepPole - small).length <= (sphereKeepPole - large).length;
        final e = b.edges[eids.first.$1];
        final picked = sph.loopFillPole(
          e.curve.evaluate(e.curve.t0),
          e.curve.tangent(e.curve.t0),
          coedgeReversed: eids.first.$2,
          faceReversed: piece.reversed,
          small: small,
          large: large,
        );
        final gotSmall = (picked - small).length <= (picked - large).length;
        if (gotSmall != wantSmall) {
          eids = [for (final x in eids.reversed) (x.$1, !x.$2)];
        }
      }
      Vector3 n;
      if (poly.length >= 3) {
        n = (poly[1] - poly[0]).cross(poly[2] - poly[0]);
      } else if (piece.arc12 != null) {
        final mid = piece.arc12!.evaluate(
          0.5 * (piece.arc12!.t0 + piece.arc12!.t1),
        );
        n = (poly[1] - poly[0]).cross(mid - poly[0]);
      } else if (piece.outerArcs.isNotEmpty) {
        final arc = piece.outerArcs.first;
        final mid = arc.evaluate(0.5 * (arc.t0 + arc.t1));
        n = (poly[1] - poly[0]).cross(mid - poly[0]);
      } else {
        continue;
      }
      if (n.length < 1e-12 && piece.surface is! BrepSphere) continue;
      final surf = piece.surface is BrepPlane || piece.surface is BrepSphere
          ? piece.surface
          : BrepPlane(origin: poly[0], normal: n.normalized());
      final loops = <BrepLoop>[
        BrepLoop([for (final e in eids) BrepCoedge(e.$1, reversed: e.$2)]),
      ];
      for (final circle in piece.holeCircles) {
        final seam = weld(circle.evaluate(0));
        final e = sew(seam, seam, circle);
        loops.add(BrepLoop([BrepCoedge(e.$1, reversed: !e.$2)]));
      }
      for (final curve in piece.innerClosed) {
        final seam = weld(curve.evaluate(curve.t0));
        final e = sew(seam, seam, curve);
        loops.add(BrepLoop([BrepCoedge(e.$1, reversed: !e.$2)]));
      }
      for (var hi = 0; hi < piece.holes.length; hi++) {
        final hole = piece.holes[hi];
        final curves = hi < piece.holeCurves.length
            ? piece.holeCurves[hi]
            : const <BrepCurve>[];
        final hasArc = curves.any((c) => c is BrepCircle);
        if (hole.length < 2) continue;
        if (hole.length < 3 && !hasArc) continue;
        final hid = [for (final p in hole) weld(p)];
        final heids = <(int, bool)>[];
        for (var i = 0; i < hid.length; i++) {
          final j = (i + 1) % hid.length;
          final curve = i < curves.length
              ? curves[i]
              : BrepLine(b.vertices[hid[i]], b.vertices[hid[j]]);
          heids.add(sew(hid[i], hid[j], curve));
        }
        loops.add(
          BrepLoop([
            for (final e in heids.reversed) BrepCoedge(e.$1, reversed: !e.$2),
          ]),
        );
      }
      if (surf is BrepPlane) mergeRoundHoles(loops);
      b.face(surf, loops, reversed: piece.reversed);
    }
    return b.build();
  }
}

class _Piece {
  _Piece(
    this.surface,
    this.poly,
    this.reversed, {
    this.uv,
    this.holes = const [],
    this.holeCircles = const [],
    this.holeCurves = const [],
    this.arc12,
    this.outerArcs = const [],
    this.outerCurves = const [],
    this.innerClosed = const [],
    this.sourceSolid,
    this.sourceFace,
  });
  final BrepSurface surface;
  final List<Vector3> poly;
  final bool reversed;
  final (double, double, double, double)? uv;
  final List<List<Vector3>> holes;
  final List<BrepCircle> holeCircles;
  final List<List<BrepCurve>> holeCurves;
  final BrepCircle? arc12;
  final List<BrepCircle> outerArcs;
  final List<BrepCurve> outerCurves;
  final List<BrepCurve> innerClosed;
  final BrepSolid? sourceSolid;
  final BrepFace? sourceFace;
}

class _ImprintSeg {
  _ImprintSeg(this.a, this.b, this.curve);
  final Vector3 a;
  final Vector3 b;
  final BrepCurve curve;
  Vector3 get mid {
    if (curve is BrepCircle) {
      final c = curve as BrepCircle;
      return c.evaluate(0.5 * (c.t0 + c.t1));
    }
    return (a + b) * 0.5;
  }
}

class _ImprintOuter {
  const _ImprintOuter({
    required this.pts,
    this.curves = const [],
    this.holes = const [],
    this.holeCurves = const [],
    // Read by the callers of _overlayLoops, but _overlayLoops cannot yet
    // produce analytic hole circles, so it is always the default.
    // ignore: unused_element_parameter
    this.holeCircles = const [],
    this.outerArcs = const [],
  });
  final List<Vector3> pts;
  final List<BrepCurve> curves;
  final List<List<Vector3>> holes;
  final List<List<BrepCurve>> holeCurves;
  final List<BrepCircle> holeCircles;
  final List<BrepCircle> outerArcs;
}

class _PrismImprint {
  const _PrismImprint({this.drop = false, this.outers = const []});
  final bool drop;
  final List<_ImprintOuter> outers;
}

/// Edges used by more than two faces (inscribed-cylinder pinches) are
/// duplicated so adjacent corners become separate solids.
BrepSolid _disconnectNonManifold(BrepSolid s, {bool splitPinches = true}) {
  final edgeFaces = <int, List<int>>{};
  for (var f = 0; f < s.faces.length; f++) {
    for (final loop in s.faces[f].loops) {
      for (final c in loop.coedges) {
        (edgeFaces[c.edge] ??= []).add(f);
      }
    }
  }
  final crowded = [
    for (final e in edgeFaces.entries)
      if (e.value.toSet().length > 2) e.key,
  ];
  final pinches = <int>[];
  for (final e in edgeFaces.entries) {
    if (!splitPinches) break;
    if (crowded.contains(e.key)) continue;
    final uses = e.value.toSet().toList();
    if (uses.length != 2) continue;
    if (s.faces[uses[0]].surface is BrepPlane &&
        s.faces[uses[1]].surface is BrepPlane) {
      continue;
    }
    if (_facesPinchAcrossEdge(s, e.key, uses[0], uses[1])) {
      pinches.add(e.key);
    }
  }
  if (crowded.isEmpty && pinches.isEmpty) return s;

  final b = BrepBuilder();
  for (final v in s.vertices) {
    b.vertex(v);
  }
  for (final e in s.edges) {
    b.edge(e.v0, e.v1, e.curve);
  }
  // First pair uses by connectivity away from contact edges. This also
  // handles entirely planar/faceted cutters: each separate corner shell has
  // two uses of its contact edge, even though four uses share its coordinates.
  final neighbors = List.generate(s.faces.length, (_) => <int>{});
  for (final entry in edgeFaces.entries) {
    if (crowded.contains(entry.key)) continue;
    final fs = entry.value.toSet().toList();
    if (fs.length != 2) continue;
    neighbors[fs[0]].add(fs[1]);
    neighbors[fs[1]].add(fs[0]);
  }
  final components = List<int>.filled(s.faces.length, -1);
  for (var i = 0; i < s.faces.length; i++) {
    if (components[i] != -1) continue;
    final pending = [i];
    components[i] = i;
    while (pending.isNotEmpty) {
      for (final next in neighbors[pending.removeLast()]) {
        if (components[next] != -1) continue;
        components[next] = i;
        pending.add(next);
      }
    }
  }
  final replace = <String, int>{};
  for (final eid in crowded) {
    final uses = edgeFaces[eid]!.toSet().toList();
    if (uses.length <= 2) continue;
    final e = s.edges[eid];
    final groups = <int, List<int>>{};
    for (final fi in uses) {
      (groups[components[fi]] ??= []).add(fi);
    }
    if (groups.length > 1 && groups.values.every((fs) => fs.length == 2)) {
      for (final pair in groups.values.skip(1)) {
        final id = b.edge(e.v0, e.v1, e.curve, origin: e.origin);
        for (final fi in pair) replace['$fi-$eid'] = id;
      }
      continue;
    }
    if (uses.every((fi) => s.faces[fi].surface is BrepPlane)) {
      final pairs = planarContactFacePairs(s, {for (final fi in uses) fi: eid});
      if (pairs != null) {
        for (final pair in pairs.skip(1)) {
          final id = b.edge(e.v0, e.v1, e.curve, origin: e.origin);
          for (final fi in pair) replace['$fi-$eid'] = id;
        }
      }
      continue;
    }
    final mid = (s.vertices[e.v0] + s.vertices[e.v1]) * 0.5;
    var axis = s.vertices[e.v1] - s.vertices[e.v0];
    if (axis.length2 < 1e-20) {
      if (e.curve is BrepLine) {
        axis = (e.curve as BrepLine).b - (e.curve as BrepLine).a;
      } else if (e.curve is BrepCircle) {
        axis = (e.curve as BrepCircle).normal;
      }
    }
    if (axis.length2 < 1e-20) continue;
    axis.normalize();
    var ref = axis.cross(Vector3(0, 0, 1));
    if (ref.length2 < 1e-12) ref = axis.cross(Vector3(0, 1, 0));
    ref.normalize();
    final bin = axis.cross(ref);
    final ranked = <(double, int)>[];
    for (final fi in uses) {
      final c = BrepClassify.faceCentroid(s, s.faces[fi]);
      final d = c - mid;
      final dp = d - axis * d.dot(axis);
      ranked.add((math.atan2(dp.dot(bin), dp.dot(ref)), fi));
    }
    ranked.sort((a, b) => a.$1.compareTo(b.$1));
    // Pair each plane with the nearest curved face so a corner stays
    // watertight; never leave two wall-halves on the same edge copy.
    final used = <int>{};
    var pairCount = 0;
    for (var i = 0; i < ranked.length; i++) {
      final fi = ranked[i].$2;
      if (used.contains(fi)) continue;
      if (s.faces[fi].surface is! BrepPlane) continue;
      int? best;
      var bestD = double.infinity;
      for (var j = 0; j < ranked.length; j++) {
        final fj = ranked[j].$2;
        if (fj == fi || used.contains(fj)) continue;
        if (s.faces[fj].surface is BrepPlane) continue;
        var d = (ranked[j].$1 - ranked[i].$1).abs();
        if (d > math.pi) d = math.pi * 2 - d;
        if (d < bestD) {
          bestD = d;
          best = fj;
        }
      }
      if (best == null) continue;
      used.add(fi);
      used.add(best);
      if (pairCount > 0) {
        final nid = b.edge(e.v0, e.v1, e.curve);
        replace['$fi-$eid'] = nid;
        replace['$best-$eid'] = nid;
      }
      pairCount++;
    }
    for (final fi in uses) {
      if (used.contains(fi)) continue;
      replace['$fi-$eid'] = b.edge(e.v0, e.v1, e.curve);
    }
  }

  for (final eid in pinches) {
    final uses = edgeFaces[eid]!.toSet().toList();
    if (uses.length != 2) continue;
    final e = s.edges[eid];
    final nid = b.edge(e.v0, e.v1, e.curve);
    replace['${uses[1]}-$eid'] = nid;
  }

  for (var fi = 0; fi < s.faces.length; fi++) {
    final f = s.faces[fi];
    b.face(
      f.surface,
      [
        for (final loop in f.loops)
          BrepLoop([
            for (final c in loop.coedges)
              BrepCoedge(
                replace['$fi-${c.edge}'] ?? c.edge,
                reversed: c.reversed,
              ),
          ]),
      ],
      reversed: f.reversed,
      origin: f.origin,
    );
  }
  return b.build();
}

/// True when two faces that share [eid] sit on opposite sides of a knife-edge
/// (inscribed cylinder tangent), not a normal wall junction.
bool _facesPinchAcrossEdge(BrepSolid s, int eid, int f0, int f1) {
  // Same-facing normals also occur at ordinary smooth seams. Their boundary
  // traversals are opposite; separating them would tear the solid open.
  final aUses = [
    for (final l in s.faces[f0].loops)
      for (final c in l.coedges)
        if (c.edge == eid) c,
  ];
  final bUses = [
    for (final l in s.faces[f1].loops)
      for (final c in l.coedges)
        if (c.edge == eid) c,
  ];
  if (aUses.length == 1 &&
      bUses.length == 1 &&
      aUses.single.reversed != bUses.single.reversed) {
    return false;
  }

  final e = s.edges[eid];
  final mid = (s.vertices[e.v0] + s.vertices[e.v1]) * 0.5;
  Vector3? outward(BrepFace f) {
    final uv = f.surface.parameters(mid);
    if (uv == null) return null;
    final n = f.outwardNormal(uv.$1, uv.$2);
    if (n.length2 < 1e-20) return null;
    return n.normalized();
  }

  final a = outward(s.faces[f0]);
  final b = outward(s.faces[f1]);
  if (a == null || b == null) return false;
  // Tangent pinch: both outwards point the same way. A 90°/concave
  // D-notch wall pair has near-orthogonal outwards.
  return a.dot(b) > 0.7;
}

/// Connected components by face adjacency (shared edges).
List<BrepSolid> brepSplitSolids(BrepSolid solid) {
  if (solid.faces.length <= 1) return [solid];
  final adj = <int, Set<int>>{};
  final edgeFaces = <int, List<int>>{};
  for (var f = 0; f < solid.faces.length; f++) {
    for (final loop in solid.faces[f].loops) {
      for (final c in loop.coedges) {
        (edgeFaces[c.edge] ??= []).add(f);
      }
    }
  }
  for (final fs in edgeFaces.values) {
    for (var i = 0; i < fs.length; i++) {
      for (var j = i + 1; j < fs.length; j++) {
        (adj[fs[i]] ??= {}).add(fs[j]);
        (adj[fs[j]] ??= {}).add(fs[i]);
      }
    }
  }
  final seen = <int>{};
  final groups = <List<int>>[];
  for (var i = 0; i < solid.faces.length; i++) {
    if (!seen.add(i)) continue;
    final stack = [i];
    final g = <int>[];
    while (stack.isNotEmpty) {
      final f = stack.removeLast();
      g.add(f);
      for (final n in adj[f] ?? const <int>{}) {
        if (seen.add(n)) stack.add(n);
      }
    }
    groups.add(g);
  }
  if (groups.length <= 1) return [solid];
  final shells = [for (final g in groups) _extractFaces(solid, g)];
  final volumes = [
    for (final shell in shells) MeshKit.signedVolume(BrepTessellate.run(shell)),
  ];
  if (volumes.any((volume) => !volume.isFinite || volume == 0)) {
    throw BrepUnsupported('Cannot classify a collapsed or nonfinite shell');
  }
  final exteriors = [
    for (var i = 0; i < groups.length; i++)
      if (volumes[i] > 0) i,
  ];
  if (exteriors.isEmpty)
    throw BrepUnsupported('Solid has no outward exterior shell');
  final combined = {
    for (final i in exteriors) i: [...groups[i]],
  };
  for (var i = 0; i < groups.length; i++) {
    if (volumes[i] >= 0) continue;
    final shell = shells[i];
    final probe = shell.vertices.isNotEmpty
        ? shell.vertices.first
        : shell.faces.first.surface.evaluate(0.7, 1.1);
    int? enclosing;
    for (final j in exteriors) {
      if ((enclosing == null || volumes[j] < volumes[enclosing]) &&
          BrepClassify.isInside(shells[j], probe))
        enclosing = j;
    }
    if (enclosing == null)
      throw BrepUnsupported('An inward cavity shell has no enclosing solid');
    combined[enclosing]!.addAll(groups[i]);
  }
  return [for (final indices in combined.values) _extractFaces(solid, indices)];
}

BrepSolid _extractFaces(BrepSolid s, List<int> faceIdx) {
  final b = BrepBuilder();
  final vMap = <int, int>{};
  int mapV(int i) => vMap.putIfAbsent(i, () => b.vertex(s.vertices[i]));
  final eMap = <int, int>{};
  int mapE(int i) {
    return eMap.putIfAbsent(i, () {
      final e = s.edges[i];
      return b.edge(mapV(e.v0), mapV(e.v1), e.curve, origin: e.origin);
    });
  }

  for (final fi in faceIdx) {
    final f = s.faces[fi];
    b.face(
      f.surface,
      [
        for (final loop in f.loops)
          BrepLoop([
            for (final c in loop.coedges)
              BrepCoedge(mapE(c.edge), reversed: c.reversed),
          ]),
      ],
      reversed: f.reversed,
      origin: f.origin,
    );
  }
  return b.build();
}
