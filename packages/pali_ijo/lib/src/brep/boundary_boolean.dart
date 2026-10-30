import 'dart:math' as math;
import 'package:vector_math/vector_math_64.dart';
import 'geom.dart';
import 'solid.dart';
import 'heal.dart';
import 'intersect.dart';
import 'classify.dart';
import 'tessellate.dart';
import 'provenance.dart';
import 'spline_convert.dart';
import 'spline_intersect.dart';
import 'curve_equivalence.dart';
import '../measure/bounds.dart';
import '../ijo_log.dart';

/// Boundary construction from shared, finite face/face intersections.
/// Kept separately from the rectangular-patch fast path while it is validated.
class BrepBoundaryBoolean {
  // Cache only during one synchronous operation, never across model edits.
  static Expando<Map<BrepFace, List<List<Vector2>>>>? _uvLoops;
  static Expando<BrepCurveLocator>? _curveLocators;
  static double _closestParam(BrepCurve c, Vector3 p) {
    final cache = _curveLocators;
    if (cache == null ||
        c is! BrepSplineCurve &&
            c is! BrepAffineCurve &&
            c is! BrepSurfaceIntersectionCurve)
      return c.closestParam(p);
    return (cache[c] ??= BrepCurveLocator(c)).closestParam(p);
  }

  static List<List<Vector2>> _faceLoops(BrepSolid s, BrepFace f) {
    List<List<Vector2>> build() => [
      for (final l in f.loops)
        _uv(f.surface, BrepTessellate.sampleLoop(s, l, 128)),
    ];
    final scope = _uvLoops;
    if (scope == null) return build();
    return (scope[s] ??= {}).putIfAbsent(f, build);
  }

  static BrepSolid subtract(BrepSolid a, BrepSolid b) => run(a, b, 1);

  static BrepSolid run(BrepSolid a, BrepSolid b, int op) {
    final previous = _uvLoops, previousLocators = _curveLocators;
    _uvLoops = Expando();
    _curveLocators = Expando();
    try {
      return BrepClassify.withTrimCache(() => _run(a, b, op));
    } finally {
      _uvLoops = previous;
      _curveLocators = previousLocators;
    }
  }

  static BrepSolid _run(BrepSolid a, BrepSolid b, int op) {
    a = BrepHeal.alignCircularTrims(a);
    b = BrepHeal.alignCircularTrims(b);
    final af = _faces(a, b), bf = _faces(b, a);
    final additions = <BrepFace, List<BrepCurve>>{};
    final vertices = <Vector3>[];
    final boundsA = {for (final f in af.faces) f: _faceBounds(af, f)},
        boundsB = {for (final f in bf.faces) f: _faceBounds(bf, f)};
    final opposed = <BrepFace, List<BrepFace>>{};
    for (final self in [af, bf]) {
      _imprintOpposedPlanarFaces(self, additions, vertices, opposed);
    }
    for (final fa in af.faces) {
      for (final fb in bf.faces) {
        final ba = boundsA[fa], bb = boundsB[fb];
        if (ba != null && bb != null && !ba.intersects(bb, tolerance: 1e-7))
          continue;
        if ((fa.surface is BrepSplineSurface ||
                fb.surface is BrepSplineSurface) &&
            ba != null &&
            bb != null &&
            _boundaryBoundsContact(ba, bb)) {
          continue;
        }
        if (BrepSplineIntersect.boundaryContactOnly(fa.surface, fb.surface)) {
          continue;
        }
        final curves = _intersections(fa.surface, fb.surface);
        for (final rawCurve in curves) {
          final curve = rawCurve is BrepLine
              ? _boundedLine(rawCurve, ba, bb)
              : rawCurve;
          if (curve == null) continue;
          final cuts = <double>[curve.t0, curve.t1];
          for (final entry in [(af, fa), (bf, fb)]) {
            for (final boundary in _constraints(entry.$1, entry.$2)) {
              cuts.addAll(_roots(curve, boundary));
            }
          }
          final ts = _unique(cuts);
          for (var i = 1; i < ts.length; i++) {
            final lo = ts[i - 1], hi = ts[i];
            if (hi - lo < 1e-10) continue;
            final p = curve.evaluate((lo + hi) / 2);
            if (!_contains(af, fa, p) || !_contains(bf, fb, p)) continue;
            final segment = _trim(curve, lo, hi);
            if ((segment.evaluate(segment.t0) - segment.evaluate(segment.t1))
                        .length <
                    1e-6 &&
                (segment.evaluate(segment.t0) -
                            segment.evaluate((segment.t0 + segment.t1) / 2))
                        .length <
                    1e-6)
              continue;
            (additions[fa] ??= []).add(segment);
            (additions[fb] ??= []).add(segment);
            vertices.addAll([
              segment.evaluate(segment.t0),
              segment.evaluate(segment.t1),
            ]);
          }
        }
      }
    }
    // Existing seam vertices must split their adjoining cap edges too.
    vertices.addAll(af.vertices);
    vertices.addAll(bf.vertices);
    final builder = BrepBuilder();
    final builtEdges = <(int, int), List<(int, int)>>{};
    final builtVertices = <Vector3>[];
    int vertex(Vector3 p) {
      for (var i = 0; i < builtVertices.length; i++) {
        if ((builtVertices[i] - p).length < 1e-6) return i;
      }
      builtVertices.add(p);
      return builder.vertex(p);
    }

    BrepCoedge edge(BrepCurve c, Vector3 start) {
      final p = c.evaluate(c.t0), q = c.evaluate(c.t1);
      final i = vertex(p), j = vertex(q);
      // A rounded midpoint key can put the same shared edge into adjacent
      // buckets at a half-grid boundary. Group by already-welded vertices,
      // then compare the curves geometrically before sharing their topology.
      final key = (math.min(i, j), math.max(i, j));
      final candidates = builtEdges.putIfAbsent(key, () => []);
      (int, int)? id;
      for (final candidate in candidates) {
        final old = builder.edges[candidate.$1].curve;
        final reverse = i == j
            ? old.tangent(old.t0).dot(c.tangent(c.t0)) < 0
            : candidate.$2 != i;
        if ([.25, .5, .75].every((f) {
          final t = reverse ? 1 - f : f;
          return (old.evaluate(old.t0 + (old.t1 - old.t0) * t) -
                      c.evaluate(c.t0 + (c.t1 - c.t0) * f))
                  .length <
              1e-6;
        })) {
          id = candidate;
          break;
        }
      }
      if (id == null) {
        id = (builder.edge(i, j, c), i);
        candidates.add(id);
      }
      final stored = builder.edges[id.$1].curve;
      final reverse = i == j
          ? stored.tangent(stored.t0).dot(c.tangent(c.t0)) < 0
          : (builtVertices[id.$2] - start).length > 1e-6;
      return BrepCoedge(id.$1, reversed: reverse);
    }

    for (final entry in [(af, bf, true), (bf, af, false)]) {
      final self = entry.$1, other = entry.$2, fromA = entry.$3;
      final otherBounds = conservativeSolidBounds(other);
      final selfBounds = fromA ? boundsA : boundsB;
      for (final face in self.faces) {
        final bounds = selfBounds[face];
        final disjoint =
            bounds != null &&
            otherBounds != null &&
            !bounds.intersects(otherBounds, tolerance: 1e-7);
        if ((op == 2 || op == 1 && !fromA) && disjoint) continue;
        // Untouched helical patches retain their exact trim loops, avoiding
        // expensive thread arrangements for a cut confined to the bolt head.
        // Other faces still need splitting at the axial chart seam vertices
        // introduced by _faces, even when they are far from the cutting tool.
        bool keepPoint(Vector3 p) =>
            !(opposed[face] ?? const <BrepFace>[]).any(
              (f) => _contains(self, f, p),
            ) &&
            _keep(face, p, other, fromA, op);
        final surface = face.surface;
        final untouchedThread =
            disjoint &&
            surface is BrepRuledSurface &&
            surface.lower is BrepHelix &&
            surface.upper is BrepHelix;
        if (untouchedThread ||
            surface is BrepSphere &&
                (additions[face] ?? []).isEmpty &&
                BrepTessellate.sphereFaceFrame(self, face) == null) {
          final p = face.loops.isEmpty
              ? face.surface.evaluate(.5, 1)
              : self.edges[face.loops.first.coedges.first.edge].curve.evaluate(
                  self.edges[face.loops.first.coedges.first.edge].curve.t0,
                );
          if (!disjoint && !keepPoint(p)) continue;
          final loops = <BrepLoop>[];
          for (final loop in face.loops) {
            final uses = <BrepCoedge>[];
            for (final co in loop.coedges) {
              final e = self.edges[co.edge];
              final use = edge(
                e.curve,
                self.vertices[co.reversed ? e.v1 : e.v0],
              );
              uses.add(
                e.v0 == e.v1
                    ? BrepCoedge(
                        use.edge,
                        reversed: co.reversed != use.reversed,
                      )
                    : use,
              );
            }
            loops.add(BrepLoop(uses));
          }
          builder.face(
            face.surface,
            loops,
            reversed: op == 1 && !fromA ? !face.reversed : face.reversed,
            origin: face.origin,
          );
          continue;
        }
        final segments = <BrepCurve>[];
        // Build complete face regions before deciding which side to keep.
        // Removing an original edge by its midpoint can open a loop when that
        // edge meets a coincident or tangent boundary in a later cut.
        for (final c in [
          for (final loop in face.loops)
            for (final co in loop.coedges) self.edges[co.edge].curve,
          ...additions[face] ?? <BrepCurve>[],
        ]) {
          final cuts = <double>[c.t0, c.t1];
          for (final p in vertices) {
            final t = _closestParam(c, p);
            if (t > math.min(c.t0, c.t1) + 1e-9 &&
                t < math.max(c.t0, c.t1) - 1e-9 &&
                (c.evaluate(t) - p).length < 1e-6)
              cuts.add(t);
          }
          final ts = _unique(cuts);
          for (var i = 1; i < ts.length; i++) {
            segments.add(_trim(c, ts[i - 1], ts[i]));
          }
        }
        final unique = <BrepCurve>[];
        for (final c in segments) {
          final p = c.evaluate(c.t0),
              q = c.evaluate(c.t1),
              m = c.evaluate((c.t0 + c.t1) / 2);
          if (unique.any(
            (e) =>
                ((e.evaluate(e.t0) - p).length < 1e-6 &&
                        (e.evaluate(e.t1) - q).length < 1e-6 ||
                    (e.evaluate(e.t1) - p).length < 1e-6 &&
                        (e.evaluate(e.t0) - q).length < 1e-6) &&
                (e.evaluate((e.t0 + e.t1) / 2) - m).length < 1e-5,
          ))
            continue;
          unique.add(c);
        }
        if (unique.isEmpty) {
          if (face.loops.isEmpty) {
            final p = face.surface.evaluate(.5, 1);
            if (keepPoint(p))
              builder.face(
                face.surface,
                [],
                reversed: op == 1 && !fromA ? !face.reversed : face.reversed,
                origin: face.origin,
              );
          }
          continue;
        }
        final loops = _chains(
          unique,
          face.surface,
          discardContact: (c) {
            final t = (c.t0 + c.t1) / 2, mid = c.evaluate(t);
            final uv = face.surface.parameters(mid)!;
            final across = face.surface
                .normalAt(uv.$1, uv.$2)
                .cross(c.tangent(t))
                .normalized();
            for (final distance in [1e-4, 1e-5]) {
              for (final sign in [-1, 1]) {
                final probe = face.surface.closestPoint(
                  mid + across * (sign * distance),
                );
                if (_contains(self, face, probe) && keepPoint(probe)) {
                  return false;
                }
              }
            }
            return true;
          },
        );
        final mapped =
            <
              ({
                List<(BrepCurve, Vector3)> chain,
                List<Vector2> uv,
                double area,
              })
            >[];
        for (final chain in loops) {
          final points = <Vector3>[];
          for (final (c, start) in chain) {
            final rev = (c.evaluate(c.t0) - start).length > 1e-6;
            final n = c is BrepLine ? 1 : 24;
            for (var k = 0; k < n; k++) {
              final t = rev
                  ? c.t1 + (c.t0 - c.t1) * k / n
                  : c.t0 + (c.t1 - c.t0) * k / n;
              points.add(c.evaluate(t));
            }
          }
          final frame = BrepTessellate.sphereFaceFrame(self, face);
          final uv = frame == null
              ? _uv(face.surface, points)
              : [
                  for (final p in points)
                    BrepTessellate.sphereProject(frame, p),
                ];
          if ((face.surface is BrepCylinder || face.surface is BrepCone) &&
              uv.isNotEmpty) {
            final domain = face.surface is BrepCylinder
                ? BrepTessellate.cylinderFaceDomain(self, face)
                : BrepTessellate.coneFaceDomain(self, face);
            final center = (domain.$1 + domain.$2) / 2;
            final mean = uv.fold(0.0, (double sum, p) => sum + p.x) / uv.length;
            final shift =
                ((center - mean) / (2 * math.pi)).round() * 2 * math.pi;
            for (final p in uv) {
              p.x += shift;
            }
          }
          var area = 0.0;
          for (var i = 0; i < uv.length; i++) {
            final p = uv[i], q = uv[(i + 1) % uv.length];
            area += p.x * q.y - q.x * p.y;
          }
          // A probe immediately beside a tangent intersection can still be
          // within the other surface's boundary tolerance. Try other edges
          // of this region before deciding from coincident-face orientation.
          bool? retain;
          bool? coincident;
          for (final pair in chain) {
            final c = pair.$1;
            final mid = c.evaluate((c.t0 + c.t1) / 2);
            var tangent = c.tangent((c.t0 + c.t1) / 2);
            if ((c.evaluate(c.t0) - pair.$2).length > 1e-6) tangent = -tangent;
            final params = face.surface.parameters(mid)!;
            final inward = face.surface
                .normalAt(params.$1, params.$2)
                .cross(tangent)
                .normalized();
            final probe = face.surface.closestPoint(mid + inward * 1e-4);
            if ((probe - mid).dot(inward) < 5e-5 ||
                !_contains(self, face, probe)) {
              continue;
            }
            final keep = keepPoint(probe);
            if (BrepClassify.side(other, probe) == BrepSide.on) {
              coincident ??= keep;
            } else {
              retain = keep;
              break;
            }
          }
          if (!(retain ?? coincident ?? false)) continue;
          mapped.add((chain: chain, uv: uv, area: area / 2));
        }
        mapped.sort((x, y) => y.area.abs().compareTo(x.area.abs()));
        final used = <int>{};
        for (var i = 0; i < mapped.length; i++) {
          if (used.contains(i) || mapped[i].area <= 1e-10) continue;
          final outer = mapped[i];
          final grouped = [outer];
          used.add(i);
          for (var j = i + 1; j < mapped.length; j++) {
            if (!used.contains(j) &&
                mapped[j].area < 0 &&
                mapped[j].area.abs() < outer.area - 1e-9 &&
                !mapped[j].chain.any(
                  (holeEdge) => outer.chain.any(
                    (outerEdge) => identical(holeEdge.$1, outerEdge.$1),
                  ),
                ) &&
                _inPoly(outer.uv, mapped[j].uv.first)) {
              grouped.add(mapped[j]);
              used.add(j);
            }
          }
          final reverse = op == 1 && !fromA ? !face.reversed : face.reversed;
          final outLoops = <BrepLoop>[];
          for (var k = 0; k < grouped.length; k++) {
            final item = grouped[k];
            final positive = k == 0 ? !reverse : reverse;
            var chain = item.chain;
            if ((item.area > 0) != positive) {
              chain = [
                for (final pair in chain.reversed)
                  (
                    pair.$1,
                    (pair.$1.evaluate(pair.$1.t0) - pair.$2).length < 1e-6
                        ? pair.$1.evaluate(pair.$1.t1)
                        : pair.$1.evaluate(pair.$1.t0),
                  ),
              ];
            }
            outLoops.add(
              BrepLoop([for (final pair in chain) edge(pair.$1, pair.$2)]),
            );
          }
          builder.face(
            face.surface,
            outLoops,
            reversed: reverse,
            origin: face.origin,
          );
        }
      }
    }
    return BrepProvenance.inherit(builder.build(), [
      a,
      b,
    ], op: 'boundary boolean $op');
  }

  /// Face extensions can create an internal overlap between opposite caps
  /// already in the same body. Imprint those finite planar boundaries before
  /// classifying regions, so a later cut removes the buried patches on both
  /// sides instead of leaving four uses of its new boundary edge.
  static void _imprintOpposedPlanarFaces(
    BrepSolid solid,
    Map<BrepFace, List<BrepCurve>> additions,
    List<Vector3> vertices,
    Map<BrepFace, List<BrepFace>> opposed,
  ) {
    final faces = solid.faces
        .where(
          (f) =>
              f.surface is BrepPlane &&
              f.loops.isNotEmpty &&
              f.loops.every(
                (l) => l.coedges.every(
                  (c) => solid.edges[c.edge].curve is BrepLine,
                ),
              ),
        )
        .toList();
    final bounds = {for (final f in faces) f: _faceBounds(solid, f)};
    for (var i = 0; i < faces.length; i++) {
      final a = faces[i], pa = a.surface as BrepPlane;
      for (final b in faces.skip(i + 1)) {
        final pb = b.surface as BrepPlane;
        if (a.outwardNormal(0, 0).dot(b.outwardNormal(0, 0)) > -1 + 1e-9 ||
            pa.height(pb.origin).abs() > 1e-8)
          continue;
        final ba = bounds[a], bb = bounds[b];
        if (ba != null && bb != null && !ba.intersects(bb, tolerance: 1e-7))
          continue;
        (opposed[a] ??= []).add(b);
        (opposed[b] ??= []).add(a);
        for (final pair in [(a, b), (b, a)]) {
          final constraints = _constraints(solid, pair.$2);
          for (final loop in pair.$1.loops) {
            for (final co in loop.coedges) {
              final curve = solid.edges[co.edge].curve;
              final cuts = _unique([
                curve.t0,
                curve.t1,
                for (final constraint in constraints)
                  ..._roots(curve, constraint),
              ]);
              for (var k = 1; k < cuts.length; k++) {
                final lo = cuts[k - 1], hi = cuts[k];
                if (hi - lo < 1e-10 ||
                    !_contains(solid, pair.$2, curve.evaluate((lo + hi) / 2)))
                  continue;
                final segment = _trim(curve, lo, hi);
                final start = segment.evaluate(segment.t0),
                    end = segment.evaluate(segment.t1);
                if ((end - start).length < 1e-6) continue;
                (additions[pair.$2] ??= []).add(segment);
                vertices.addAll([start, end]);
              }
            }
          }
        }
      }
    }
  }

  // Surface intersections can return long supporting lines. Restrict the
  // root search to the finite face pair before resolving its trim crossings.
  static BrepLine? _boundedLine(BrepLine line, IjoBounds? a, IjoBounds? b) {
    var lo = 0.0, hi = 1.0;
    final d = line.b - line.a;
    for (final box in [a, b]) {
      if (box == null) continue;
      for (var k = 0; k < 3; k++) {
        final min = box.min[k], max = box.max[k];
        if (d[k].abs() < 1e-12 * d.length) {
          if (line.a[k] < min - 1e-7 || line.a[k] > max + 1e-7) return null;
          continue;
        }
        final t0 = (min - line.a[k]) / d[k], t1 = (max - line.a[k]) / d[k];
        lo = math.max(lo, math.min(t0, t1));
        hi = math.min(hi, math.max(t0, t1));
        if (hi <= lo) return null;
      }
    }
    return BrepLine(line.evaluate(lo), line.evaluate(hi));
  }

  static bool _keep(
    BrepFace face,
    Vector3 p,
    BrepSolid other,
    bool fromA,
    int op,
  ) {
    // Coincident planar caps have an exact occupancy decision. A ray test
    // through a thread endpoint can land on several helical seams and miss
    // the solid just behind that cap.
    if (face.surface is BrepPlane) {
      final n = face.outwardNormal(0, 0);
      for (final candidate in other.faces) {
        final plane = candidate.surface;
        if (plane is! BrepPlane || plane.height(p).abs() > 1e-8) continue;
        final dot = candidate.outwardNormal(0, 0).dot(n);
        if (dot.abs() < 1 - 1e-9 || !_contains(other, candidate, p)) continue;
        if (op == 0) return dot > 0 && fromA;
        if (op == 2) return dot > 0 && fromA;
        return dot < 0 && fromA;
      }
    }
    final side = BrepClassify.side(other, p);
    if (side != BrepSide.on)
      return op == 0
          ? side == BrepSide.outside
          : op == 2
          ? side == BrepSide.inside
          : fromA
          ? side == BrepSide.outside
          : side == BrepSide.inside;
    final uv = face.surface.parameters(p)!;
    final n = face.outwardNormal(uv.$1, uv.$2);
    final back = BrepClassify.side(other, p - n * 1e-5) == BrepSide.inside;
    final front = BrepClassify.side(other, p + n * 1e-5) == BrepSide.inside;
    return op == 0
        ? !front && (fromA || !back)
        : op == 2
        ? back && (fromA || front)
        : fromA
        ? !back
        // Opposing coincident faces keep A's existing boundary. Emitting
        // the reversed tool face there too would duplicate the same area.
        : front && back;
  }

  // Coordinate-aligned edges can also bound a notched patch or a hole.
  // Replacing either with its UV bounding rectangle would fill missing area.
  static bool _hasRectangularAxialTrim(BrepSolid s, BrepFace f) {
    final cyl = f.surface;
    final axis = cyl is BrepCylinder ? cyl.axis : (cyl as BrepCone).axis;
    final origin = cyl is BrepCylinder ? cyl.origin : (cyl as BrepCone).origin;
    if (f.loops.length == 2 &&
        f.loops.every(
          (l) =>
              l.coedges.length == 1 &&
              s.edges[l.coedges.single.edge].curve is BrepCircle &&
              (s.edges[l.coedges.single.edge].curve as BrepCircle).isClosed,
        )) {
      return true;
    }
    if (f.loops.length != 1) return false;
    // An axial wall may encode its cap circles in one loop, linked
    // by a seam used in both directions. That seam is a chart cut, not trim.
    final uses = <int, List<BrepCoedge>>{};
    for (final co in f.loops.single.coedges) {
      (uses[co.edge] ??= []).add(co);
    }
    final boundary = uses.entries
        .where(
          (e) =>
              !(e.value.length == 2 &&
                  e.value[0].reversed != e.value[1].reversed),
        )
        .toList();
    if ((boundary.length == 2 ||
            boundary.length == 1 && cyl is BrepCone && cyl.topRadius == 0) &&
        boundary.every(
          (e) =>
              e.value.length == 1 &&
              s.edges[e.key].curve is BrepCircle &&
              (s.edges[e.key].curve as BrepCircle).isClosed,
        )) {
      return true;
    }
    final (u0, u1, z0, z1) = cyl is BrepCylinder
        ? BrepTessellate.cylinderFaceDomain(s, f)
        : BrepTessellate.coneFaceDomain(s, f);
    for (final co in f.loops.single.coedges) {
      final c = s.edges[co.edge].curve;
      if (c is BrepCircle) {
        final z = (c.center - origin).dot(axis);
        if ((z - z0).abs() > 1e-7 && (z - z1).abs() > 1e-7) return false;
      } else if (c is BrepLine) {
        for (final p in [c.a, c.b]) {
          final uv = cyl.parameters(p)!;
          final d0 = (p - cyl.evaluate(u0, uv.$2)).length;
          final d1 = (p - cyl.evaluate(u1, uv.$2)).length;
          if (math.min(d0, d1) > 1e-7) return false;
        }
      } else {
        return false;
      }
    }
    return true;
  }

  /// Great-circle chart seams bound each spherical patch away from parameter
  /// poles. Clip existing cap/band trims too, so later bores preserve them.
  static void _appendSpherePatches(BrepBuilder b, BrepSolid s, BrepFace face) {
    final sphere = face.surface as BrepSphere;
    final constraints = _constraints(s, face);
    for (final sx in [-1.0, 1.0]) {
      for (final sy in [-1.0, 1.0]) {
        for (final sz in [-1.0, 1.0]) {
          var directions = [
            Vector3(sx, 0, 0),
            Vector3(0, sy, 0),
            Vector3(0, 0, sz),
          ];
          if (sx * sy * sz < 0) directions = directions.reversed.toList();
          final frame = BrepPlane(
            origin: sphere.center,
            normal: Vector3(sx, sy, sz),
          );
          bool inOctant(Vector3 p) {
            final d = p - sphere.center;
            return d.x * sx >= -1e-8 && d.y * sy >= -1e-8 && d.z * sz >= -1e-8;
          }

          final candidates = <BrepCurve>[];
          for (final entry in <(BrepCurve, bool)>[
            for (final loop in face.loops)
              for (final co in loop.coedges) (s.edges[co.edge].curve, true),
            for (var i = 0; i < 3; i++)
              (
                BrepCircle(
                  center: sphere.center,
                  normal: directions[i].cross(directions[(i + 1) % 3]),
                  radius: sphere.radius,
                  xDir: directions[i],
                  t0: 0,
                  t1: math.pi / 2,
                ),
                false,
              ),
          ]) {
            final c = entry.$1;
            final cuts = <double>[c.t0, c.t1];
            if (entry.$2) {
              for (var axis = 0; axis < 3; axis++) {
                cuts.addAll(_roots(c, (p) => (p - sphere.center)[axis]));
              }
            } else {
              for (final constraint in constraints)
                cuts.addAll(_roots(c, constraint));
            }
            final ts = _unique(cuts);
            for (var i = 1; i < ts.length; i++) {
              final mid = c.evaluate((ts[i - 1] + ts[i]) / 2);
              if (!inOctant(mid) || !_contains(s, face, mid)) continue;
              final part = _trim(c, ts[i - 1], ts[i]);
              final p = part.evaluate(part.t0), q = part.evaluate(part.t1);
              if (candidates.any(
                (old) =>
                    ((old.evaluate(old.t0) - p).length < 1e-6 &&
                            (old.evaluate(old.t1) - q).length < 1e-6 ||
                        (old.evaluate(old.t0) - q).length < 1e-6 &&
                            (old.evaluate(old.t1) - p).length < 1e-6) &&
                    (old.evaluate((old.t0 + old.t1) / 2) - mid).length < 1e-6,
              ))
                continue;
              candidates.add(part);
            }
          }
          if (candidates.isEmpty) continue;
          final kept =
              <
                ({
                  List<(BrepCurve, Vector3)> chain,
                  List<Vector2> uv,
                  double area,
                })
              >[];
          for (final chain in _chains(
            candidates,
            sphere,
            discardContact: (c) {
              final mid = c.evaluate((c.t0 + c.t1) / 2),
                  across = (mid - sphere.center).normalized().cross(
                    c.tangent((c.t0 + c.t1) / 2),
                  );
              return [-1, 1].every((sign) {
                final probe = sphere.closestPoint(mid + across * (sign * 1e-4));
                return !inOctant(probe) || !_contains(s, face, probe);
              });
            },
          )) {
            final first = chain.first,
                c = first.$1,
                t = (c.t0 + c.t1) / 2,
                mid = c.evaluate(t);
            var tangent = c.tangent(t);
            if ((c.evaluate(c.t0) - first.$2).length > 1e-6) tangent = -tangent;
            final probe = sphere.closestPoint(
              mid + (mid - sphere.center).normalized().cross(tangent) * 1e-4,
            );
            if (!inOctant(probe) || !_contains(s, face, probe)) continue;
            final uv = <Vector2>[];
            for (final (curve, start) in chain) {
              final rev = (curve.evaluate(curve.t0) - start).length > 1e-6;
              for (var k = 0; k < 24; k++) {
                final t = rev
                    ? curve.t1 + (curve.t0 - curve.t1) * k / 24
                    : curve.t0 + (curve.t1 - curve.t0) * k / 24;
                uv.add(BrepTessellate.sphereProject(frame, curve.evaluate(t)));
              }
            }
            var area = 0.0;
            for (var i = 0; i < uv.length; i++)
              area +=
                  uv[i].x * uv[(i + 1) % uv.length].y -
                  uv[i].y * uv[(i + 1) % uv.length].x;
            kept.add((chain: chain, uv: uv, area: area / 2));
          }
          final outers = kept.where((p) => p.area > 1e-10).toList()
            ..sort((a, b) => a.area.compareTo(b.area));
          final holes = kept.where((p) => p.area < -1e-10).toList();
          for (final outer in outers) {
            final group = [outer];
            for (final hole in [...holes]) {
              if (_inPoly(outer.uv, hole.uv.first)) {
                group.add(hole);
                holes.remove(hole);
              }
            }
            final loops = <BrepLoop>[];
            for (final part in group) {
              final uses = <BrepCoedge>[];
              for (final (curve, start) in part.chain) {
                final a = b.vertex(curve.evaluate(curve.t0)),
                    z = b.vertex(curve.evaluate(curve.t1));
                uses.add(
                  BrepCoedge(
                    b.edge(a, z, curve),
                    reversed: (curve.evaluate(curve.t0) - start).length > 1e-6,
                  ),
                );
              }
              loops.add(
                BrepLoop(
                  face.reversed
                      ? [
                          for (final co in uses.reversed)
                            BrepCoedge(co.edge, reversed: !co.reversed),
                        ]
                      : uses,
                ),
              );
            }
            b.face(sphere, loops, reversed: face.reversed, origin: face.origin);
          }
          if (holes.isNotEmpty)
            throw BrepUnsupported('Spherical trim has an unbounded chart hole');
        }
      }
    }
  }

  // Finite rational charts retain the exact torus while avoiding ambiguous
  // periodic branches in the face arrangement after a general plane cut.
  static bool _appendTorusPatches(
    BrepBuilder b,
    BrepSolid solid,
    BrepFace face,
  ) {
    final tor = face.surface as BrepTorus;
    for (final loop in face.loops)
      for (final co in loop.coedges) {
        final curve = solid.edges[co.edge].curve;
        if (curve is! BrepCircle) return false;
        final uv = [
          for (final t in [0.0, .25, .5, .75, 1.0])
            tor.parameters(
              curve.evaluate(curve.t0 + (curve.t1 - curve.t0) * t),
            )!,
        ];
        bool constant(bool u) => uv.every(
          (q) =>
              math
                  .sin(
                    ((u ? q.$1 : q.$2) - (u ? uv.first.$1 : uv.first.$2)) / 2,
                  )
                  .abs() <
              1e-7,
        );
        if (!constant(true) && !constant(false)) return false;
      }
    final domain =
        BrepTessellate.torusFaceDomain(solid, tor, face) ??
        (0.0, 2 * math.pi, 0.0, 2 * math.pi);
    final nu = math.max(1, ((domain.$2 - domain.$1) / (math.pi / 2)).ceil());
    final nv = math.max(1, ((domain.$4 - domain.$3) / (math.pi / 2)).ceil());
    BrepCircle unit(double a, double z) => BrepCircle(
      center: Vector3.zero(),
      normal: Vector3(0, 0, 1),
      xDir: Vector3(1, 0, 0),
      radius: 1,
      t0: a,
      t1: z,
    );
    BrepCircle latitude(double v, double a, double z) => BrepCircle(
      center: tor.center + tor.axis * (tor.minorRadius * math.sin(v)),
      normal: tor.axis,
      xDir: tor.xDir,
      radius: tor.majorRadius + tor.minorRadius * math.cos(v),
      t0: a,
      t1: z,
    );
    BrepCircle meridian(double u, double a, double z) {
      final radial = tor.xDir * math.cos(u) + tor.yDir * math.sin(u);
      return BrepCircle(
        center: tor.center + radial * tor.majorRadius,
        normal: radial.cross(tor.axis),
        xDir: radial,
        radius: tor.minorRadius,
        t0: a,
        t1: z,
      );
    }

    for (var i = 0; i < nu; i++)
      for (var j = 0; j < nv; j++) {
        final u0 = domain.$1 + (domain.$2 - domain.$1) * i / nu,
            u1 = domain.$1 + (domain.$2 - domain.$1) * (i + 1) / nu;
        final v0 = domain.$3 + (domain.$4 - domain.$3) * j / nv,
            v1 = domain.$3 + (domain.$4 - domain.$3) * (j + 1) / nv;
        final u = BrepSplineConvert.curve(unit(u0, u1)),
            v = BrepSplineConvert.curve(unit(v0, v1));
        final patch = BrepSplineSurface(
          degreeU: 2,
          degreeV: 2,
          knotsU: u.knots,
          knotsV: v.knots,
          points: [
            for (final a in u.points)
              [
                for (final c in v.points)
                  tor.center +
                      (tor.xDir * a.x + tor.yDir * a.y) *
                          (tor.majorRadius + tor.minorRadius * c.x) +
                      tor.axis * (tor.minorRadius * c.y),
              ],
          ],
          weights: [
            for (final a in u.weights) [for (final c in v.weights) a * c],
          ],
        );
        final curves = [
          latitude(v0, u0, u1),
          meridian(u1, v0, v1),
          latitude(v1, u0, u1),
          meridian(u0, v0, v1),
        ];
        final uses = <BrepCoedge>[];
        for (var k = 0; k < 4; k++) {
          final c = curves[k];
          uses.add(
            BrepCoedge(
              b.edge(b.vertex(c.evaluate(c.t0)), b.vertex(c.evaluate(c.t1)), c),
              reversed: k >= 2,
            ),
          );
        }
        b.face(
          patch,
          [BrepLoop(uses)],
          reversed: face.reversed,
          origin: face.origin,
        );
      }
    return true;
  }

  static bool _appendClosedRuledPatches(
    BrepBuilder b,
    BrepSolid s,
    BrepFace f,
  ) {
    final surface = f.surface;
    if (surface is! BrepRuledSurface) return false;
    final uses = f.loops.expand((l) => l.coedges).toList();
    final seams = uses
        .where((co) => s.edges[co.edge].curve is BrepLine)
        .toList();
    final separateRims =
        f.loops.length == 2 && f.loops.every((l) => l.coedges.length == 1);
    final sewnSeam =
        f.loops.length == 1 &&
        uses.length == 4 &&
        seams.length == 2 &&
        seams[0].edge == seams[1].edge &&
        seams[0].reversed != seams[1].reversed;
    if (!separateRims && !sewnSeam) return false;
    final lower = surface.lower, upper = surface.upper;
    if (lower.evaluate(lower.t0).distanceTo(lower.evaluate(lower.t1)) > 1e-8 ||
        upper.evaluate(upper.t0).distanceTo(upper.evaluate(upper.t1)) > 1e-8)
      return false;
    bool matches(BrepCurve a, BrepCurve b) =>
        sameRationalCurve(a, b) || sameRationalCurve(a, b, reversed: true);
    final rims = [
      for (final co in uses)
        if (s.edges[co.edge].curve is! BrepLine) s.edges[co.edge].curve,
    ];
    if (rims.length != 2) return false;
    if (!(matches(rims[0], lower) && matches(rims[1], upper) ||
        matches(rims[1], lower) && matches(rims[0], upper)))
      return false;
    // Two closed rims have zero area in a periodic UV chart. Give the
    // arrangement finite charts just as for circular cylinders, preserving
    // exact rational rails and shared seam points instead of dropping the wall.
    for (var i = 0; i < 4; i++) {
      BrepCurve part(BrepCurve c) => _trim(
        c,
        c.t0 + (c.t1 - c.t0) * i / 4,
        c.t0 + (c.t1 - c.t0) * (i + 1) / 4,
      );
      final a = part(lower), z = part(upper);
      final p = [
        a.evaluate(a.t0),
        a.evaluate(a.t1),
        z.evaluate(z.t1),
        z.evaluate(z.t0),
      ];
      final ids = p.map(b.vertex).toList();
      b.face(
        BrepRuledSurface(a, z),
        [
          BrepLoop([
            BrepCoedge(b.edge(ids[0], ids[1], a)),
            BrepCoedge(b.edge(ids[1], ids[2], BrepLine(p[1], p[2]))),
            BrepCoedge(b.edge(ids[3], ids[2], z), reversed: true),
            BrepCoedge(b.edge(ids[3], ids[0], BrepLine(p[3], p[0]))),
          ]),
        ],
        reversed: f.reversed,
        origin: f.origin,
      );
    }
    return true;
  }

  static BrepSolid _faces(BrepSolid s, BrepSolid other) {
    final b = BrepBuilder();
    for (final v in s.vertices) {
      b.vertex(v);
    }
    for (final e in s.edges) {
      b.edge(e.v0, e.v1, e.curve, origin: e.origin);
    }
    for (final f in s.faces) {
      final surface = f.surface;
      if (_appendClosedRuledPatches(b, s, f)) continue;
      if (surface is BrepTorus && _appendTorusPatches(b, s, f)) continue;
      if (surface is BrepSphere &&
          BrepTessellate.sphereFaceFrame(s, f) == null &&
          other.faces.any((of) => facesIntersect(s, f, other, of))) {
        if (f.loops.any(
          (l) =>
              l.coedges.length != 1 ||
              s.edges[l.coedges.single.edge].curve is! BrepCircle ||
              !(s.edges[l.coedges.single.edge].curve as BrepCircle).isClosed,
        )) {
          throw BrepUnsupported(
            'Spherical trim crosses an unresolved chart boundary',
          );
        }
        _appendSpherePatches(b, s, f);
        continue;
      }
      final rectangular =
          (surface is BrepCylinder || surface is BrepCone) &&
          _hasRectangularAxialTrim(s, f) &&
          f.loops.every(
            (l) => l.coedges.every(
              (c) =>
                  s.edges[c.edge].curve is BrepLine ||
                  s.edges[c.edge].curve is BrepCircle &&
                      (s.edges[c.edge].curve as BrepCircle).normal
                              .cross(
                                surface is BrepCylinder
                                    ? surface.axis
                                    : (surface as BrepCone).axis,
                              )
                              .length <
                          1e-6,
            ),
          );
      if (!rectangular) {
        b.face(surface, f.loops, reversed: f.reversed, origin: f.origin);
        continue;
      }
      final cyl = surface;
      final axis = cyl is BrepCylinder ? cyl.axis : (cyl as BrepCone).axis;
      final origin = cyl is BrepCylinder
          ? cyl.origin
          : (cyl as BrepCone).origin;
      final xDir = cyl is BrepCylinder ? cyl.xDir : (cyl as BrepCone).xDir;
      double radiusAt(double z) =>
          cyl is BrepCylinder ? cyl.radius : (cyl as BrepCone).radiusAt(z);
      final (u0, u1, z0, z1) = cyl is BrepCylinder
          ? BrepTessellate.cylinderFaceDomain(s, f)
          : BrepTessellate.coneFaceDomain(s, f);
      // Rebuilding also visits faces untouched by the other solid. A complete
      // axial wall bounded by two circles has zero-area loops in its UV chart;
      // give it finite chart patches too, or the arrangement drops the wall.
      final n = math.max(1, ((u1 - u0) / math.pi).ceil());
      for (var i = 0; i < n; i++) {
        final lo = u0 + (u1 - u0) * i / n, hi = u0 + (u1 - u0) * (i + 1) / n;
        final p = [
          cyl.evaluate(lo, z0),
          cyl.evaluate(hi, z0),
          cyl.evaluate(hi, z1),
          cyl.evaluate(lo, z1),
        ];
        final ids = p.map(b.vertex).toList();
        final curves = <BrepCurve>[
          BrepCircle(
            center: origin + axis * z0,
            normal: axis,
            radius: radiusAt(z0),
            xDir: xDir,
            t0: lo,
            t1: hi,
          ),
          BrepLine(p[1], p[2]),
          BrepCircle(
            center: origin + axis * z1,
            normal: axis,
            radius: radiusAt(z1),
            xDir: xDir,
            t0: lo,
            t1: hi,
          ),
          BrepLine(p[3], p[0]),
        ];
        final co = <BrepCoedge>[];
        if ((p[1] - p[0]).length > 1e-8) {
          co.add(BrepCoedge(b.edge(ids[0], ids[1], curves[0])));
        }
        co.add(BrepCoedge(b.edge(ids[1], ids[2], curves[1])));
        if ((p[3] - p[2]).length > 1e-8) {
          co.add(BrepCoedge(b.edge(ids[3], ids[2], curves[2]), reversed: true));
        }
        co.add(BrepCoedge(b.edge(ids[3], ids[0], curves[3])));
        b.face(cyl, [BrepLoop(co)], reversed: f.reversed, origin: f.origin);
      }
    }
    return b.build();
  }

  // Opposite support bounds meet only on face boundaries. Existing boundary
  // edges are already split and classified by the collector below.
  static bool _boundaryBoundsContact(IjoBounds a, IjoBounds b) {
    for (var i = 0; i < 3; i++) {
      if ((a.max[i] - b.min[i]).abs() < 1e-8 ||
          (b.max[i] - a.min[i]).abs() < 1e-8) {
        return true;
      }
    }
    return false;
  }

  static IjoBounds? _faceBounds(BrepSolid solid, BrepFace face) {
    if (face.surface is BrepRuledSurface) {
      final s = face.surface as BrepRuledSurface;
      // Positive rational weights enclose conics in their control hull.
      if ([s.lower, s.upper].every(
        (c) =>
            c is BrepLine ||
            c is BrepCircle ||
            c is BrepEllipse ||
            c is BrepSplineCurve,
      )) {
        return IjoBounds.fromPoints([
          ...BrepSplineConvert.curve(s.lower).points,
          ...BrepSplineConvert.curve(s.upper).points,
        ]);
      }
    }
    return conservativeSolidBounds(
      BrepSolid(vertices: solid.vertices, edges: solid.edges, faces: [face]),
    );
  }

  static bool facesIntersect(
    BrepSolid a,
    BrepFace fa,
    BrepSolid b,
    BrepFace fb,
  ) {
    final ba = _faceBounds(a, fa), bb = _faceBounds(b, fb);
    if (ba != null && bb != null && !ba.intersects(bb, tolerance: 1e-7))
      return false;
    try {
      for (final curve in _intersections(fa.surface, fb.surface)) {
        final cuts = <double>[curve.t0, curve.t1];
        for (final boundary in [
          ..._constraints(a, fa),
          ..._constraints(b, fb),
        ]) {
          cuts.addAll(_roots(curve, boundary));
        }
        final ts = _unique(cuts);
        for (var i = 1; i < ts.length; i++) {
          final p = curve.evaluate((ts[i - 1] + ts[i]) / 2);
          if (_contains(a, fa, p) && _contains(b, fb, p)) return true;
        }
      }
    } on BrepUnsupported {
      return false;
    }
    return false;
  }

  static List<BrepCurve> _intersections(BrepSurface a, BrepSurface b) {
    if (a is BrepRuledSurface && b is BrepRuledSurface) {
      final helices = _helicalPair(a, b);
      if (helices != null) return helices;
      Vector3? extrusionAxis(BrepRuledSurface s) {
        final shift =
            s.upper.evaluate(s.upper.t0) - s.lower.evaluate(s.lower.t0);
        if (shift.length < 1e-9 ||
            !sameRationalCurve(s.lower, s.upper, translation: shift)) {
          return null;
        }
        return shift.normalized();
      }

      final ax = extrusionAxis(a), bx = extrusionAxis(b);
      if (ax != null &&
          bx != null &&
          ax.cross(bx).length < 1e-10 &&
          (sameRationalCurve(a.lower, b.lower, ignoreAxis: ax) ||
              sameRationalCurve(
                a.lower,
                b.lower,
                ignoreAxis: ax,
                reversed: true,
              ))) {
        // Coincident swept walls have no isolated intersection curve. Their
        // finite cap boundaries split them before material-side classification.
        return [];
      }
    }
    if (a is BrepRuledSurface && b is BrepPlane) return _ruledPlane(a, b);
    if (b is BrepRuledSurface && a is BrepPlane) return _ruledPlane(b, a);
    if (a is BrepRuledSurface && BrepQuadric.supports(b))
      return BrepIntersect.ruledQuadricBranches(a, b);
    if (b is BrepRuledSurface && BrepQuadric.supports(a))
      return BrepIntersect.ruledQuadricBranches(b, a);
    if (a is BrepSplineSurface ||
        b is BrepSplineSurface ||
        a is BrepRuledSurface && b is BrepRuledSurface)
      return BrepSplineIntersect.intersect(a, b);
    if (a is BrepRuledSurface || b is BrepRuledSurface)
      throw BrepUnsupported(
        'This ruled-surface pair needs a curved intersection solver',
      );
    if (a is BrepCylinder &&
        b is BrepCylinder &&
        a.axis.cross(b.axis).length > 1e-6) {
      return a.radius <= b.radius
          ? BrepIntersect.quadricBranches(a, b)
          : BrepIntersect.quadricBranches(b, a);
    }
    if (a is BrepCylinder && b is BrepSphere)
      return BrepIntersect.quadricBranches(a, b);
    if (b is BrepCylinder && a is BrepSphere)
      return BrepIntersect.quadricBranches(b, a);
    return BrepIntersect.surfaces(a, b);
  }

  /// Equal-lead coaxial screw surfaces reduce to two straight profile
  /// segments in (radius, z-lead*angle/2pi). Their crossings sweep helices;
  /// disjoint profiles cannot acquire an intersection by subdivision.
  static List<BrepCurve>? _helicalPair(BrepRuledSurface a, BrepRuledSurface b) {
    final rails = [a.lower, a.upper, b.lower, b.upper];
    if (rails.any((c) => c is! BrepHelix)) return null;
    final p = rails.cast<BrepHelix>(),
        al = p[0],
        au = p[1],
        bl = p[2],
        bu = p[3];
    if (p.any((c) => c.radiusSlope != 0)) return null;
    final n = al.axis, h = al.lead / (2 * math.pi);
    if (p.any(
      (c) =>
          (c.axis - n).length > 1e-10 ||
          (c.lead - al.lead).abs() > 1e-10 ||
          (c.origin - al.origin).cross(n).length > 1e-8,
    ))
      return null;
    for (final pair in [(al, au), (bl, bu)]) {
      if ((pair.$1.xDir - pair.$2.xDir).length > 1e-10 ||
          (pair.$1.t0 - pair.$2.t0).abs() > 1e-10 ||
          (pair.$1.t1 - pair.$2.t1).abs() > 1e-10)
        return null;
    }
    final delta = math.atan2(bl.xDir.dot(al.yDir), bl.xDir.dot(al.xDir));
    final amin = math.min(al.t0, al.t1), amax = math.max(al.t0, al.t1);
    final bmin = math.min(bl.t0, bl.t1), bmax = math.max(bl.t0, bl.t1);
    final first = ((bmin + delta - amax) / (2 * math.pi) - 1e-10).ceil();
    final last = ((bmax + delta - amin) / (2 * math.pi) + 1e-10).floor();
    final ar = au.radius - al.radius, az = (au.origin - al.origin).dot(n);
    final br = bu.radius - bl.radius, bz = (bu.origin - bl.origin).dot(n);
    final det = br * az - ar * bz;
    // Parallel profiles either miss or share a surface. The latter is resolved
    // by the same material-side classification as coincident planar faces.
    if (det.abs() <
        1e-12 *
            math.max(1, math.sqrt((ar * ar + az * az) * (br * br + bz * bz))))
      return [];
    final out = <BrepCurve>[];
    for (var k = first; k <= last; k++) {
      final lo = math.max(amin, bmin + delta - 2 * math.pi * k),
          hi = math.min(amax, bmax + delta - 2 * math.pi * k);
      if (hi - lo < 1e-10) continue;
      final dr = bl.radius - al.radius,
          dz = (bl.origin - al.origin).dot(n) - h * delta + al.lead * k;
      final v = (-dr * bz + br * dz) / det, w = (ar * dz - az * dr) / det;
      if (v < -1e-8 || v > 1 + 1e-8 || w < -1e-8 || w > 1 + 1e-8) continue;
      out.add(
        BrepHelix(
          origin: al.origin + n * (az * v),
          axis: n,
          xDir: al.xDir,
          radius: al.radius + ar * v,
          lead: al.lead,
          t0: lo,
          t1: hi,
        ),
      );
    }
    return out;
  }

  static List<BrepCurve> _ruledPlane(BrepRuledSurface s, BrepPlane p) {
    final hull = [
      ...BrepSplineConvert.curve(s.lower).points,
      ...BrepSplineConvert.curve(s.upper).points,
    ];
    if (hull.every((q) => p.height(q) > 1e-8) ||
        hull.every((q) => p.height(q) < -1e-8))
      return [];
    final curve = BrepPlaneSectionCurve(surface: s, plane: p);
    final cuts = <double>[0, 1];
    // Rulings parallel to the plane either miss it or lie in it. An isolated
    // coincident ruling is itself a straight section; a coincident patch has
    // no isolated curve and is handled by neighborhood classification.
    final den = [for (var i = 0; i <= 16; i++) s.dv(i / 16, 0).dot(p.normal)];
    if (den.every((v) => v.abs() < 1e-10)) {
      final lower = s.lower, roots = _roots(lower, p.height);
      return [
        for (final t in roots)
          () {
            final u = (t - lower.t0) / (lower.t1 - lower.t0);
            return BrepLine(s.evaluate(u, 0), s.evaluate(u, 1));
          }(),
      ];
    }
    for (final c in [s.lower, s.upper]) {
      for (final t in _roots(c, p.height)) cuts.add((t - c.t0) / (c.t1 - c.t0));
    }
    // A denominator zero partitions the rational section into separate branches.
    for (var i = 1; i < den.length; i++) {
      if (den[i - 1] * den[i] < 0) {
        var a = (i - 1) / 16, b = i / 16, fa = den[i - 1];
        for (var k = 0; k < 45; k++) {
          final m = (a + b) / 2, fm = s.dv(m, 0).dot(p.normal);
          if (fa * fm <= 0) {
            b = m;
          } else {
            a = m;
            fa = fm;
          }
        }
        cuts.add((a + b) / 2);
      }
    }
    final sorted = _unique(cuts), out = <BrepCurve>[];
    for (var i = 1; i < sorted.length; i++) {
      final a = sorted[i - 1], b = sorted[i], v = curve.ruling((a + b) / 2);
      if (v >= -1e-8 && v <= 1 + 1e-8) {
        if (s.dv(a, 0).dot(p.normal).abs() < 1e-10 ||
            s.dv(b, 0).dot(p.normal).abs() < 1e-10)
          throw BrepUnsupported('Section reaches a singular coincident ruling');
        out.add(curve.trimmed(a, b));
      }
    }
    return out;
  }

  static List<double Function(Vector3)> _constraints(BrepSolid s, BrepFace f) {
    final out = <double Function(Vector3)>[];
    for (final loop in f.loops) {
      for (final co in loop.coedges) {
        final c = s.edges[co.edge].curve;
        if (c is BrepSurfaceIntersectionCurve) {
          final mid = c.evaluate((c.t0 + c.t1) / 2),
              uv = f.surface.parameters(mid)!;
          final n = f.surface.normalAt(uv.$1, uv.$2),
              p = c.parametersAt((c.t0 + c.t1) / 2);
          final a = BrepSurfaceIntersectionCurve.normalizedJet(
                c.first,
                p.x,
                p.y,
              ),
              b = BrepSurfaceIntersectionCurve.normalizedJet(
                c.second,
                p.z,
                p.w,
              );
          final target =
              n.dot(a.$2.cross(a.$3).normalized()).abs() >
                  n.dot(b.$2.cross(b.$3).normalized()).abs()
              ? c.second
              : c.first;
          out.add((p) {
            final uv = target.parameters(p), q = target.evaluate(uv.$1, uv.$2);
            return (p - q).dot(target.normalAt(uv.$1, uv.$2));
          });
        } else if (c is BrepRuledQuadricCurve) {
          if (BrepProvenance.surfaceKey(f.surface) ==
              BrepProvenance.surfaceKey(c.surface)) {
            out.add((p) => BrepQuadric.value(c.target, p));
          } else {
            out.add((p) {
              final q = c.surface.closestPoint(p), uv = c.surface.parameters(q);
              return (p - q).dot(c.surface.normalAt(uv.$1, uv.$2));
            });
          }
        } else if (c is BrepPlaneSectionCurve) {
          if (f.surface is BrepPlane) {
            out.add((p) {
              final q = c.surface.closestPoint(p), uv = c.surface.parameters(q);
              return (p - q).dot(c.surface.normalAt(uv.$1, uv.$2));
            });
          } else {
            out.add(c.plane.height);
          }
        } else if (c is BrepSplineCurve ||
            c is BrepAffineCurve ||
            c is BrepHelix) {
          // Imported trims need not be coordinate lines. Their signed local
          // side is determined by the actual boundary tangent on this face.
          out.add((p) {
            final t = _closestParam(c, p),
                q = c.evaluate(t),
                uv = f.surface.parameters(q)!;
            final inward = f.surface.normalAt(uv.$1, uv.$2).cross(c.tangent(t));
            return (p - q).dot(inward);
          });
        } else if (c is BrepLine) {
          final mid = (c.a + c.b) / 2;
          final uv = f.surface.parameters(mid)!;
          final n = (c.b - c.a)
              .cross(f.surface.normalAt(uv.$1, uv.$2))
              .normalized();
          out.add((p) => (p - c.a).dot(n));
        } else if (c is BrepCircle) {
          if (f.surface is BrepPlane) {
            out.add((p) => (p - c.center).length2 - c.radius * c.radius);
          } else {
            out.add((p) => (p - c.center).dot(c.normal));
          }
        } else if (c is BrepEllipse) {
          if (f.surface is BrepPlane) {
            out.add((p) {
              final d = p - c.center;
              return math.pow(d.dot(c.xDir) / c.rMajor, 2) +
                  math.pow(d.dot(c.yDir) / c.rMinor, 2) -
                  1;
            });
          } else {
            out.add((p) => (p - c.center).dot(c.normal));
          }
        } else {
          final intersection = c as BrepIntersectionCurve;
          final isCarrier =
              BrepProvenance.surfaceKey(f.surface) ==
              BrepProvenance.surfaceKey(intersection.carrier);
          final boundary = isCarrier
              ? intersection.target
              : intersection.carrier;
          if (boundary is BrepCylinder) {
            out.add((p) {
              final d = p - boundary.origin;
              final radial = d - boundary.axis * d.dot(boundary.axis);
              return radial.length2 - boundary.radius * boundary.radius;
            });
          } else {
            final sphere = boundary as BrepSphere;
            out.add(
              (p) =>
                  (p - sphere.center).length2 - sphere.radius * sphere.radius,
            );
          }
        }
      }
    }
    return out;
  }

  static List<double> _roots(BrepCurve c, double Function(Vector3) f) {
    final roots = <double>[];
    const n = 128;
    var lo = c.t0, flo = f(c.evaluate(lo));
    // A coincident boundary contributes no isolated crossing parameters.
    if ([
      0.0,
      .17,
      .39,
      .61,
      .83,
      1.0,
    ].every((t) => f(c.evaluate(c.t0 + (c.t1 - c.t0) * t)).abs() < 1e-9))
      return roots;
    for (var i = 1; i <= n; i++) {
      final hi = c.t0 + (c.t1 - c.t0) * i / n, fhi = f(c.evaluate(hi));
      if (flo.abs() < 1e-9) roots.add(lo);
      // A sample already accepted as zero must not also create a second
      // crossing from floating-point sign noise at a tangent contact.
      if (flo * fhi < 0 && flo.abs() >= 1e-9 && fhi.abs() >= 1e-9) {
        var a = lo, b = hi, fa = flo;
        for (var k = 0; k < 45; k++) {
          final m = (a + b) / 2, fm = f(c.evaluate(m));
          if (fa * fm <= 0) {
            b = m;
          } else {
            a = m;
            fa = fm;
          }
        }
        roots.add((a + b) / 2);
      }
      lo = hi;
      flo = fhi;
    }
    if (flo.abs() < 1e-9) roots.add(lo);
    return roots;
  }

  static List<double> _unique(List<double> xs) {
    xs.sort();
    final out = <double>[];
    for (final x in xs) {
      if (out.isEmpty || x - out.last > 1e-9) out.add(x);
    }
    return out;
  }

  static BrepCurve _trim(BrepCurve c, double lo, double hi) {
    if (c is BrepHelix) return c.trimmed(lo, hi);
    if (c is BrepSurfaceIntersectionCurve) return c.trimmed(lo, hi);
    if (c is BrepRuledQuadricCurve) return c.trimmed(lo, hi);
    if (c is BrepPlaneSectionCurve) return c.trimmed(lo, hi);
    if (c is BrepSplineCurve) return c.trimmed(lo, hi);
    if (c is BrepAffineCurve)
      return BrepAffineCurve(_trim(c.base, lo, hi), c.transform);
    if (c is BrepLine) return BrepLine(c.evaluate(lo), c.evaluate(hi));
    if (c is BrepCircle)
      return BrepCircle(
        center: c.center,
        normal: c.normal,
        radius: c.radius,
        xDir: c.xDir,
        t0: lo,
        t1: hi,
      );
    if (c is BrepEllipse)
      return BrepEllipse(
        center: c.center,
        normal: c.normal,
        xDir: c.xDir,
        rMajor: c.rMajor,
        rMinor: c.rMinor,
        t0: lo,
        t1: hi,
      );
    return (c as BrepIntersectionCurve).trimmed(lo, hi);
  }

  static List<Vector2> _uv(BrepSurface s, List<Vector3> points) => [
    for (final p in BrepTessellate.axialLoopParameters(s, points))
      Vector2(p.x, p.y),
  ];

  static bool _contains(BrepSolid s, BrepFace f, Vector3 p) {
    final q = f.surface.parameters(p)!;
    if ((f.surface.evaluate(q.$1, q.$2) - p).length > 1e-5) return false;
    if (f.surface is BrepSphere)
      return (BrepClassify.closestOnFace(s, f, p) - p).length < 1e-5;
    if (f.loops.isEmpty) return true;
    if (f.surface is BrepCone && _hasRectangularAxialTrim(s, f)) {
      final d = BrepTessellate.coneFaceDomain(s, f);
      var u = q.$1;
      u += (((d.$1 + d.$2) / 2 - u) / (2 * math.pi)).round() * 2 * math.pi;
      return u >= d.$1 - 1e-8 &&
          u <= d.$2 + 1e-8 &&
          q.$2 >= d.$3 - 1e-8 &&
          q.$2 <= d.$4 + 1e-8;
    }
    final planar = BrepClassify.containsPlaneFace(s, f, p);
    if (planar != null) return planar;
    // Polygonal UV loops are a broad classification aid. Near a curved trim,
    // their chord error can exceed the side probe and accept both orientations
    // of the same loop. Use the exact edge and its oriented local side there.
    double best = 1e-3;
    bool? nearSide;
    for (final loop in f.loops)
      for (final co in loop.coedges) {
        final c = s.edges[co.edge].curve;
        if (c is! BrepHelix &&
            c is! BrepSplineCurve &&
            c is! BrepAffineCurve &&
            c is! BrepSurfaceIntersectionCurve &&
            c is! BrepRuledQuadricCurve &&
            c is! BrepPlaneSectionCurve)
          continue;
        final t = _closestParam(c, p),
            point = c.evaluate(t),
            distance = (p - point).length;
        if (distance >= best) continue;
        if (distance < 1e-8) return true;
        final lo = math.min(c.t0, c.t1), hi = math.max(c.t0, c.t1);
        if (t <= lo + (hi - lo) * 1e-7 || t >= hi - (hi - lo) * 1e-7) continue;
        final uv = f.surface.parameters(point)!;
        var tangent = c.tangent(t);
        if (co.reversed) tangent = -tangent;
        final inward = f.outwardNormal(uv.$1, uv.$2).cross(tangent);
        best = distance;
        nearSide = (p - point).dot(inward) >= 0;
      }
    if (nearSide != null) return nearSide;

    if (f.surface is BrepCylinder &&
        f.loops.length == 2 &&
        f.loops.every(
          (loop) =>
              loop.coedges.length == 1 &&
              s.edges[loop.coedges.single.edge].curve is BrepCircle &&
              (s.edges[loop.coedges.single.edge].curve as BrepCircle).isClosed,
        )) {
      final surface = f.surface as BrepCylinder;
      final z = surface.parameters(p)!.$2;
      final limits = [
        for (final loop in f.loops)
          surface
              .parameters(
                (s.edges[loop.coedges.single.edge].curve as BrepCircle).center,
              )!
              .$2,
      ]..sort();
      return z >= limits.first - 1e-7 && z <= limits.last + 1e-7;
    }
    var inside = false;
    for (final poly in _faceLoops(s, f)) {
      var u = q.$1;
      if (f.surface is BrepCylinder || f.surface is BrepCone) {
        final center =
            poly.map((v) => v.x).reduce((a, b) => a + b) / poly.length;
        while (u - center > math.pi) {
          u -= 2 * math.pi;
        }
        while (u - center < -math.pi) {
          u += 2 * math.pi;
        }
      }
      if (_inPoly(poly, Vector2(u, q.$2))) inside = !inside;
    }
    return inside;
  }

  static bool _inPoly(List<Vector2> poly, Vector2 p) {
    var inside = false;
    for (var i = 0; i < poly.length; i++) {
      final a = poly[i], b = poly[(i + 1) % poly.length];
      final d = b - a;
      final t = d.length2 < 1e-20
          ? 0.0
          : ((p - a).dot(d) / d.length2).clamp(0.0, 1.0);
      if ((p - a - d * t).length < 1e-7) return true;
      if ((a.y > p.y) != (b.y > p.y) &&
          p.x < a.x + (p.y - a.y) * (b.x - a.x) / (b.y - a.y))
        inside = !inside;
    }
    return inside;
  }

  static List<List<(BrepCurve, Vector3)>> _chains(
    List<BrepCurve> input,
    BrepSurface surface, {
    bool Function(BrepCurve)? discardContact,
  }) {
    final curves = <BrepCurve>[];
    final closed = <List<(BrepCurve, Vector3)>>[];
    for (final c in input) {
      final p = c.evaluate(c.t0),
          q = c.evaluate(c.t1),
          m = c.evaluate((c.t0 + c.t1) / 2);
      if ((p - q).length < 1e-6) {
        if ((p - m).length < 1e-6) continue;
        BrepCurve? reverse;
        if (c is BrepCircle)
          reverse = BrepCircle(
            center: c.center,
            normal: -c.normal,
            radius: c.radius,
            xDir: c.xDir,
            t0: -c.t1,
            t1: -c.t0,
          );
        if (c is BrepEllipse)
          reverse = BrepEllipse(
            center: c.center,
            normal: -c.normal,
            xDir: c.xDir,
            rMajor: c.rMajor,
            rMinor: c.rMinor,
            t0: -c.t1,
            t1: -c.t0,
          );
        if (reverse != null) {
          closed.add([(c, p)]);
          closed.add([(reverse, reverse.evaluate(reverse.t0))]);
        } else {
          for (var k = 0; k < 4; k++) {
            curves.add(
              _trim(
                c,
                c.t0 + (c.t1 - c.t0) * k / 4,
                c.t0 + (c.t1 - c.t0) * (k + 1) / 4,
              ),
            );
          }
        }
      } else {
        curves.add(c);
      }
    }
    final vertices = <Vector3>[];
    int vertex(Vector3 p) {
      final i = vertices.indexWhere((v) => (v - p).length < 1e-6);
      if (i >= 0) return i;
      vertices.add(p);
      return vertices.length - 1;
    }

    final endpoints = <List<int>>[];
    final outgoing = <int, List<int>>{};
    for (var i = 0; i < curves.length; i++) {
      final c = curves[i],
          a = vertex(c.evaluate(c.t0)),
          b = vertex(c.evaluate(c.t1));
      endpoints.add([a, b]);
      (outgoing[a] ??= []).add(i * 2);
      (outgoing[b] ??= []).add(i * 2 + 1);
    }
    // A coincident cut can leave a wire with no retained area on either
    // side. Remove only dangling edges whose neighbouring face regions were
    // both rejected; an unexplained open boundary must still fail.
    if (discardContact != null) {
      final removed = <int>{};
      final pending = outgoing.keys
          .where((v) => outgoing[v]!.length == 1)
          .toList();
      while (pending.isNotEmpty) {
        final v = pending.removeLast(), uses = outgoing[v]!;
        if (uses.length != 1) continue;
        final e = uses.single ~/ 2;
        if (!discardContact(curves[e])) continue;
        removed.add(e);
        for (final end in endpoints[e]) {
          outgoing[end]!.removeWhere((h) => h ~/ 2 == e);
          if (outgoing[end]!.length == 1) pending.add(end);
        }
      }
      if (removed.isNotEmpty) {
        return [
          ...closed,
          ..._chains([
            for (var i = 0; i < curves.length; i++)
              if (!removed.contains(i)) curves[i],
          ], surface),
        ];
      }
    }
    if (outgoing.values.any((v) => v.length < 2)) {
      if (IjoLog.verbose) {
        IjoLog.line('Open arrangement on ${surface.toJson()}');
        for (final e in outgoing.entries.where((e) => e.value.length < 2)) {
          final c = curves[e.value.single ~/ 2];
          IjoLog.line(
            'dangling ${vertices[e.key]} ${c.runtimeType} ${c.evaluate(c.t0)} -> ${c.evaluate(c.t1)}',
          );
        }
      }
      throw BrepUnsupported(
        'Open intersection arrangement on ${surface.runtimeType}',
      );
    }
    int origin(int h) => endpoints[h ~/ 2][h % 2];
    int end(int h) => endpoints[h ~/ 2][1 - h % 2];
    double angle(int h) {
      final c = curves[h ~/ 2], t = h.isEven ? c.t0 : c.t1;
      final p = c.evaluate(t);
      final sample = c.evaluate(t + (h.isEven ? 1 : -1) * (c.t1 - c.t0) * 1e-4);
      if (surface is BrepSphere) {
        final normal = (p - surface.center).normalized(),
            x = brepFrameU(normal);
        final delta = sample - p;
        return math.atan2(delta.dot(normal.cross(x)), delta.dot(x));
      }
      final a = surface.parameters(p)!, b = surface.parameters(sample)!;
      var du = b.$1 - a.$1;
      if (surface is BrepCylinder || surface is BrepCone) {
        while (du > math.pi) {
          du -= 2 * math.pi;
        }
        while (du < -math.pi) {
          du += 2 * math.pi;
        }
      }
      return math.atan2(b.$2 - a.$2, du);
    }

    final angles = [for (var h = 0; h < curves.length * 2; h++) angle(h)];
    for (final edges in outgoing.values) {
      edges.sort((a, b) => angles[a].compareTo(angles[b]));
    }
    final visited = <int>{}, loops = <List<(BrepCurve, Vector3)>>[];
    for (var seed = 0; seed < curves.length * 2; seed++) {
      if (visited.contains(seed)) continue;
      final chain = <(BrepCurve, Vector3)>[];
      var h = seed;
      do {
        if (!visited.add(h))
          throw BrepUnsupported('Ambiguous intersection arrangement');
        chain.add((curves[h ~/ 2], vertices[origin(h)]));
        final back = h ^ 1, candidates = outgoing[end(h)]!;
        final index = candidates.indexOf(back);
        h = candidates[(index + candidates.length - 1) % candidates.length];
      } while (h != seed);
      loops.add(chain);
    }
    return [...loops, ...closed];
  }
}
