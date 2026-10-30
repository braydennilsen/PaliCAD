import 'dart:math' as math;

import 'package:vector_math/vector_math_64.dart';

import '../measure/bounds.dart';
import 'boolean.dart';
import 'contact.dart';
import 'curve_equivalence.dart';
import 'geom.dart';
import 'heal.dart';
import 'primitives.dart';
import 'solid.dart';
import 'xform.dart';
import 'tessellate.dart';

class BrepExtrude {
  BrepExtrude._();

  /// Prism from a planar [ring] swept along [direction].
  static BrepSolid polygon(List<Vector3> ring, Vector3 direction) =>
      polygonWithHoles(ring, const [], direction);

  /// Extrude sampled curves while preserving explicitly straight source edges.
  /// Points alone cannot distinguish a regular polygon from a sampled circle.
  static BrepSolid sampledProfile(
    List<Vector3> ring,
    Vector3 direction, {
    List<List<Vector3>> holes = const [],
    List<(Vector3, Vector3)> linearEdges = const [],
    List<BrepSplineCurve> sourceSplines = const [],
  }) {
    if (sourceSplines.isNotEmpty) {
      final rings = [
        prepareRing(ring, direction),
        for (final h in holes) prepareRing(h, direction),
      ];
      final preserved = [
        for (final r in rings)
          _preserveSplines(r, sourceSplines, linearEdges, direction),
      ];
      if (preserved.any((c) => c != null)) {
        final b = BrepBuilder(), loops = <BrepLoop>[];
        final normal = _profileNormal(rings.first, direction);
        for (var j = 0; j < rings.length; j++) {
          final r = rings[j];
          final cylinder = _hasLinearEdge(r, linearEdges)
              ? null
              : BrepPrimitives.tryCylinderFromRing(r, normal);
          final circle = cylinder?.edges
              .map((e) => e.curve)
              .whereType<BrepCircle>()
              .where((c) => (c.center - r.first).dot(normal).abs() < 1e-7)
              .firstOrNull;
          final curves =
              preserved[j] ??
              (circle == null ? null : <BrepCurve>[circle]) ??
              fitCircularRuns(r, normal, linearEdges: linearEdges)?.curves ??
              [
                for (var i = 0; i < r.length; i++)
                  BrepLine(r[i], r[(i + 1) % r.length]),
              ];
          final vertices = [for (final c in curves) b.vertex(c.evaluate(c.t0))];
          loops.add(
            BrepLoop([
              for (var i = 0; i < curves.length; i++)
                BrepCoedge(
                  b.edge(
                    vertices[i],
                    vertices[(i + 1) % curves.length],
                    curves[i],
                  ),
                ),
            ]),
          );
        }
        b.face(BrepPlane(origin: rings.first.first, normal: normal), loops);
        final face = b.build();
        return fromFace(face, face.faces.single, direction);
      }
    }
    if (holes.isNotEmpty) {
      return profileWithHoles(ring, holes, direction, linearEdges: linearEdges);
    }
    if (!_hasLinearEdge(ring, linearEdges)) {
      final cylinder = BrepPrimitives.tryCylinderFromRing(ring, direction);
      if (cylinder != null) return cylinder;
    }
    return profile(ring, direction, linearEdges: linearEdges);
  }

  /// Retain whole source cubic spans only when the arrangement follows their
  /// samples in parameter order. Split/intersected regions retain the existing
  /// sampled-boundary path; this never replaces a different selected boundary.
  static List<BrepCurve>? _preserveSplines(
    List<Vector3> ring,
    List<BrepSplineCurve> sources,
    List<(Vector3, Vector3)> linearEdges,
    Vector3 direction,
  ) {
    if (ring.length < 4) return null;
    const eps = 1e-6;
    final matches = <int, (int, BrepCurve)>{};
    final protectedEdges = [...linearEdges];
    var unsafeCurve = false;
    final normal = _profileNormal(ring, direction);
    final x = normal
        .cross(normal.x.abs() < .8 ? Vector3(1, 0, 0) : Vector3(0, 1, 0))
        .normalized();
    final y = normal.cross(x);
    Vector2 project(Vector3 p) => Vector2(p.dot(x), p.dot(y));
    final polygon = ring.map(project).toList();
    bool crosses(Vector2 a, Vector2 b, Vector2 c, Vector2 d) {
      if (math.max(a.x, b.x) < math.min(c.x, d.x) ||
          math.max(c.x, d.x) < math.min(a.x, b.x) ||
          math.max(a.y, b.y) < math.min(c.y, d.y) ||
          math.max(c.y, d.y) < math.min(a.y, b.y))
        return false;
      final ab = b - a, cd = d - c, ac = c - a;
      final den = ab.x * cd.y - ab.y * cd.x;
      if (den.abs() < 1e-12) return false;
      final t = (ac.x * cd.y - ac.y * cd.x) / den;
      final u = (ac.x * ab.y - ac.y * ab.x) / den;
      return t > 1e-8 && t < 1 - 1e-8 && u > 1e-8 && u < 1 - 1e-8;
    }

    bool preservesBoundary(BrepCurve c, int first, int count) {
      // A span can pass through all its coarse samples yet loop across a
      // neighbouring branch between them. Restoring it would change the
      // selected arrangement region. Retain that run's native straight trims.
      final samples = [
        for (var k = 0; k <= 64; k++)
          project(c.evaluate(c.t0 + (c.t1 - c.t0) * k / 64)),
      ];
      for (var k = 0; k < samples.length - 1; k++) {
        final a = samples[k], b = samples[k + 1];
        for (var j = k + 2; j < samples.length - 1; j++) {
          if (crosses(a, b, samples[j], samples[j + 1])) return false;
        }
        for (var j = 0; j < ring.length; j++) {
          if ((j - first + ring.length) % ring.length < count) continue;
          if (crosses(a, b, polygon[j], polygon[(j + 1) % ring.length]))
            return false;
        }
      }
      return true;
    }

    for (var i = 0; i < ring.length; i++) {
      for (final c in sources) {
        for (final reverse in [false, true]) {
          final start = reverse ? c.t1 : c.t0, end = reverse ? c.t0 : c.t1;
          if ((c.evaluate(start) - ring[i]).length > eps) continue;
          var previous = 0.0;
          for (var count = 1; count <= ring.length; count++) {
            final p = ring[(i + count) % ring.length];
            if ((c.evaluate(end) - p).length <= eps) {
              if (count >= 3) {
                if (preservesBoundary(c, i, count)) {
                  matches[i] = (count, c.trimmed(start, end));
                } else {
                  unsafeCurve = true;
                  for (var j = 0; j < count; j++) {
                    protectedEdges.add((
                      ring[(i + j) % ring.length],
                      ring[(i + j + 1) % ring.length],
                    ));
                  }
                }
              }
              break;
            }
            final t = c.closestParam(p), u = (t - start) / (end - start);
            if (u <= previous || u >= 1 || (c.evaluate(t) - p).length > eps) {
              break;
            }
            previous = u;
          }
        }
      }
    }
    if (matches.isEmpty)
      return unsafeCurve
          ? [
              for (var i = 0; i < ring.length; i++)
                BrepLine(ring[i], ring[(i + 1) % ring.length]),
            ]
          : null;
    final start = matches.keys.first, result = <BrepCurve>[];
    var used = 0;
    while (used < ring.length) {
      final index = (start + used) % ring.length, match = matches[index];
      if (match != null) {
        if (used + match.$1 > ring.length) return null;
        result.add(match.$2);
        used += match.$1;
      } else {
        final chain = <Vector3>[ring[index]];
        do {
          used++;
          chain.add(ring[(start + used) % ring.length]);
        } while (used < ring.length &&
            !matches.containsKey((start + used) % ring.length));
        // Reuse circular fitting for the rest of a mixed arc/spline outline.
        final normal = _profileNormal(ring, direction);
        final fitted = chain.length >= 6
            ? fitCircularRuns(
                chain,
                normal,
                linearEdges: [...protectedEdges, (chain.last, chain.first)],
              )
            : null;
        if (fitted == null) {
          result.addAll([
            for (var i = 1; i < chain.length; i++)
              BrepLine(chain[i - 1], chain[i]),
          ]);
        } else {
          // The synthetic closing chord is not part of the selected boundary.
          final kept = fitted.curves
              .where(
                (c) =>
                    !((c.evaluate(c.t0) - chain.last).length <= eps &&
                        (c.evaluate(c.t1) - chain.first).length <= eps),
              )
              .toList();
          // Rotate the fitted open chain back to its original start.
          final first = kept.indexWhere(
            (c) => (c.evaluate(c.t0) - chain.first).length <= eps,
          );
          if (first < 0) return null;
          result.addAll([...kept.skip(first), ...kept.take(first)]);
        }
      }
    }
    var cursor = start;
    final safe = <BrepCurve>[];
    for (final c in result) {
      var count = 1;
      while (count <= ring.length &&
          (ring[(cursor + count) % ring.length] - c.evaluate(c.t1)).length >
              eps) {
        count++;
      }
      if (count > ring.length) return null;
      if (c is BrepCircle && !preservesBoundary(c, cursor, count)) {
        safe.addAll([
          for (var j = 0; j < count; j++)
            BrepLine(
              ring[(cursor + j) % ring.length],
              ring[(cursor + j + 1) % ring.length],
            ),
        ]);
      } else {
        safe.add(c);
      }
      cursor = (cursor + count) % ring.length;
    }
    result
      ..clear()
      ..addAll(safe);
    for (var i = 0; i < result.length; i++) {
      final next = result[(i + 1) % result.length];
      if ((result[i].evaluate(result[i].t1) - next.evaluate(next.t0)).length >
          eps) {
        return null;
      }
    }
    return result;
  }

  static bool _isLinearEdge(
    Vector3 a,
    Vector3 b,
    List<(Vector3, Vector3)> linearEdges,
  ) {
    final d = b - a;
    final length = d.length;
    if (length < 1e-12) return false;
    final axis = d / length;
    for (final (p, q) in linearEdges) {
      // Collinear segments may have been split or merged by the arrangement.
      if ((p - a).cross(axis).length > 1e-6 ||
          (q - a).cross(axis).length > 1e-6)
        continue;
      final u = (p - a).dot(axis), v = (q - a).dot(axis);
      if (math.min(length, math.max(u, v)) - math.max(0, math.min(u, v)) > 1e-8)
        return true;
    }
    return false;
  }

  static bool _hasLinearEdge(
    List<Vector3> ring,
    List<(Vector3, Vector3)> linearEdges,
  ) {
    for (var i = 0; i < ring.length; i++) {
      if (_isLinearEdge(ring[i], ring[(i + 1) % ring.length], linearEdges)) {
        return true;
      }
    }
    return false;
  }

  /// [polygon] plus circular-run fitting. Full sampled circles still go
  /// through [BrepPrimitives.tryCylinderFromRing] at the call site.
  static BrepSolid profile(
    List<Vector3> ring,
    Vector3 direction, {
    List<(Vector3, Vector3)> linearEdges = const [],
  }) {
    if (ring.length < 3 || direction.length < 1e-12) {
      throw BrepUnsupported('Extrude needs a closed profile and a direction');
    }
    final prepared = prepareRing(ring, direction);
    final fitted = fitCircularRuns(
      prepared,
      _profileNormal(prepared, direction),
      linearEdges: linearEdges,
    );
    if (fitted != null) {
      return fromCurves(fitted.pts, fitted.curves, direction);
    }
    return polygon(ring, direction);
  }

  /// Sweep fitted outer and inner boundaries together, retaining analytic
  /// cavity walls rather than replacing sampled arcs with polygon facets.
  static BrepSolid profileWithHoles(
    List<Vector3> ring,
    List<List<Vector3>> holes,
    Vector3 direction, {
    List<(Vector3, Vector3)> linearEdges = const [],
  }) {
    final outer = prepareRing(ring, direction);
    final normal = _profileNormal(outer, direction);
    final b = BrepBuilder();
    final loops = <BrepLoop>[];
    for (final input in [outer, ...holes]) {
      final points = prepareRing(input, direction);
      if (points.length < 3 ||
          points.any((p) => (p - outer.first).dot(normal).abs() > 1e-7)) {
        throw BrepUnsupported(
          'Extrude boundaries must be coplanar closed loops',
        );
      }
      final cylinder = _hasLinearEdge(points, linearEdges)
          ? null
          : BrepPrimitives.tryCylinderFromRing(points, normal);
      final circular = cylinder?.edges
          .map((e) => e.curve)
          .whereType<BrepCircle>()
          .where((c) => (c.center - outer.first).dot(normal).abs() < 1e-7)
          .firstOrNull;
      final fitted = fitCircularRuns(points, normal, linearEdges: linearEdges);
      final curves = circular != null
          ? <BrepCurve>[circular]
          : fitted?.curves ??
                <BrepCurve>[
                  for (var i = 0; i < points.length; i++)
                    BrepLine(points[i], points[(i + 1) % points.length]),
                ];
      final vertices = [for (final c in curves) b.vertex(c.evaluate(c.t0))];
      loops.add(
        BrepLoop([
          for (var i = 0; i < curves.length; i++)
            BrepCoedge(
              b.edge(vertices[i], vertices[(i + 1) % curves.length], curves[i]),
            ),
        ]),
      );
    }
    b.face(BrepPlane(origin: outer.first, normal: normal), loops);
    final source = b.build();
    return fromFace(source, source.faces.single, direction);
  }

  /// Prism with inner loops. Hole walls face the cavity (outward = into the hole).
  static BrepSolid polygonWithHoles(
    List<Vector3> ring,
    List<List<Vector3>> holes,
    Vector3 direction,
  ) {
    if (ring.length < 3 || direction.length < 1e-12) {
      throw BrepUnsupported('Extrude needs a closed profile and a direction');
    }
    final prepared = prepareRing(ring, direction);
    if (prepared.length < 3) {
      throw BrepUnsupported('Extrude profile is degenerate');
    }
    for (final hole in holes) {
      _profileNormal(hole, direction);
      if (hole.any(
        (p) =>
            (p - prepared.first)
                .dot(_profileNormal(prepared, direction))
                .abs() >
            1e-7,
      )) {
        throw BrepUnsupported('Extrude holes must lie in the profile plane');
      }
    }
    final holeRings = <List<Vector3>>[
      for (final hole in holes)
        if (hole.length >= 3) prepareRing(hole, direction),
    ];
    final n = _profileNormal(prepared, direction);
    var u = n.cross(Vector3(0, 0, 1));
    if (u.length < 1e-8) u = n.cross(Vector3(0, 1, 0));
    u.normalize();
    final b = BrepBuilder();
    final nBot = -n;

    ({List<int> eBot, List<int> eTop, List<int> eVert}) belt(
      List<Vector3> loop,
    ) {
      final bot = [for (final p in loop) b.vertex(p)];
      final top = [for (final p in loop) b.vertex(p + direction)];
      final eBot = <int>[];
      final eTop = <int>[];
      final eVert = <int>[];
      for (var i = 0; i < loop.length; i++) {
        final j = (i + 1) % loop.length;
        eBot.add(b.edge(bot[i], bot[j], BrepLine(loop[i], loop[j])));
        eTop.add(
          b.edge(
            top[i],
            top[j],
            BrepLine(loop[i] + direction, loop[j] + direction),
          ),
        );
        eVert.add(
          b.edge(bot[i], top[i], BrepLine(loop[i], loop[i] + direction)),
        );
      }
      return (eBot: eBot, eTop: eTop, eVert: eVert);
    }

    final outer = belt(prepared);
    final inners = [for (final h in holeRings) belt(h)];

    b.face(BrepPlane(origin: prepared[0], normal: nBot, xDir: u), [
      BrepLoop([
        for (var i = prepared.length - 1; i >= 0; i--)
          BrepCoedge(outer.eBot[i], reversed: true),
      ]),
      for (final h in inners)
        BrepLoop([
          for (var i = 0; i < h.eBot.length; i++) BrepCoedge(h.eBot[i]),
        ]),
    ]);
    b.face(BrepPlane(origin: prepared[0] + direction, normal: n, xDir: u), [
      BrepLoop([for (final e in outer.eTop) BrepCoedge(e)]),
      for (final h in inners)
        BrepLoop([
          for (var i = h.eTop.length - 1; i >= 0; i--)
            BrepCoedge(h.eTop[i], reversed: true),
        ]),
    ]);
    for (var i = 0; i < prepared.length; i++) {
      final j = (i + 1) % prepared.length;
      final fn = (prepared[j] - prepared[i]).cross(direction);
      if (fn.length < 1e-12) continue;
      b.face(BrepPlane(origin: prepared[i], normal: fn.normalized()), [
        BrepLoop([
          BrepCoedge(outer.eBot[i]),
          BrepCoedge(outer.eVert[j]),
          BrepCoedge(outer.eTop[i], reversed: true),
          BrepCoedge(outer.eVert[i], reversed: true),
        ]),
      ]);
    }
    for (var hi = 0; hi < holeRings.length; hi++) {
      final loop = holeRings[hi];
      final h = inners[hi];
      for (var i = 0; i < loop.length; i++) {
        final j = (i + 1) % loop.length;
        var fn = (loop[j] - loop[i]).cross(direction);
        if (fn.length < 1e-12) continue;
        fn.normalize();
        b.face(BrepPlane(origin: loop[i], normal: -fn), [
          BrepLoop([
            BrepCoedge(h.eBot[i], reversed: true),
            BrepCoedge(h.eVert[i]),
            BrepCoedge(h.eTop[i]),
            BrepCoedge(h.eVert[j], reversed: true),
          ]),
        ]);
      }
    }
    return pairPlanarLineContacts(BrepHeal.run(b.build()));
  }

  /// Slide a planar cap along its outward normal when every vertex has a
  /// generator edge parallel to that normal (box, corner, prism).
  static BrepSolid? tryMovePlanarFace(
    BrepSolid solid,
    int faceIndex,
    double distance,
  ) {
    solid = BrepXform.materializeForEditing(solid);
    final rails = _railsForPlanarFace(solid, faceIndex);
    if (rails == null) return null;
    if (distance.abs() < 1e-12) return solid;
    if (distance < 0 && -distance >= rails.minRail - 1e-6) {
      return BrepSolid(vertices: const [], edges: const [], faces: const []);
    }
    return _shiftVertices(solid, rails.faceVerts, rails.n * distance);
  }

  /// Inward union along a planar cap: no-op while the prism stays inside the
  /// solid; if it pokes out the far side, extend the opposite cap.
  static BrepSolid? tryUnionInwardFace(
    BrepSolid solid,
    int faceIndex,
    double distance,
  ) {
    if (distance >= -1e-12) return null;
    solid = BrepXform.materializeForEditing(solid);
    final rails = _railsForPlanarFace(solid, faceIndex);
    if (rails == null) return null;
    final depth = -distance;
    if (depth <= rails.minRail + 1e-6) return solid;
    final overshoot = depth - rails.minRail;
    final moved = <int>{for (final o in rails.otherOf.values) o};
    return _shiftVertices(solid, moved, rails.n * -overshoot);
  }

  static ({
    BrepFace face,
    Vector3 n,
    Set<int> faceVerts,
    Map<int, int> otherOf,
    double minRail,
  })?
  _railsForPlanarFace(BrepSolid solid, int faceIndex) {
    if (faceIndex < 0 || faceIndex >= solid.faces.length) return null;
    final face = solid.faces[faceIndex];
    if (face.surface is! BrepPlane || face.loops.isEmpty) return null;
    final pl = face.surface as BrepPlane;
    final n = (face.reversed ? -pl.normal : pl.normal).normalized();
    final faceEdges = <int>{};
    final faceVerts = <int>{};
    for (final loop in face.loops) {
      for (final c in loop.coedges) {
        faceEdges.add(c.edge);
        final e = solid.edges[c.edge];
        faceVerts.add(e.v0);
        faceVerts.add(e.v1);
      }
    }
    if (faceVerts.isEmpty) return null;

    final otherOf = <int, int>{};
    for (final vi in faceVerts) {
      int? other;
      for (var ei = 0; ei < solid.edges.length; ei++) {
        if (faceEdges.contains(ei)) continue;
        final e = solid.edges[ei];
        if (e.v0 != vi && e.v1 != vi) continue;
        if (e.curve is! BrepLine) continue;
        final o = e.v0 == vi ? e.v1 : e.v0;
        final dir = solid.vertices[o] - solid.vertices[vi];
        if (dir.length < 1e-12) continue;
        if (dir.normalized().cross(n).length > 1e-9) continue;
        other = o;
        break;
      }
      if (other == null) return null;
      otherOf[vi] = other;
    }

    // Sliding vertices preserves topology only for a complete straight prism.
    // Rails at the selected rim alone say nothing about internal shoulders,
    // cavities or other shells: a U's front rim has long rails, but its bridge
    // ends well before those rails do. Such cuts must use the boolean path.
    final opposite = otherOf.values.toSet();
    if (opposite.length != faceVerts.length ||
        {...faceVerts, ...opposite}.length != solid.vertices.length) {
      return null;
    }
    for (final f in solid.faces) {
      if (f.loops.isEmpty || f.loops.any((loop) => loop.coedges.isEmpty)) {
        return null;
      }
      final surface = f.surface;
      if (surface is BrepPlane) {
        final dot = surface.normal.dot(n).abs();
        if (dot > 1e-9 && (dot - 1).abs() > 1e-9) return null;
      } else if (surface is BrepCylinder) {
        if (surface.axis.cross(n).length > 1e-9) return null;
      } else {
        final ruled = surface is BrepRuledSurface
            ? surface
            : surface is BrepSplineSurface
            ? surface.ruled
            : null;
        if (ruled == null) return null;
        final shift =
            ruled.upper.evaluate(ruled.upper.t0) -
            ruled.lower.evaluate(ruled.lower.t0);
        if (shift.length < 1e-9 ||
            shift.normalized().cross(n).length > 1e-9 ||
            !sameRationalCurve(ruled.lower, ruled.upper, translation: shift)) {
          return null;
        }
      }
    }

    var minRail = double.infinity;
    for (final vi in faceVerts) {
      final len = (solid.vertices[otherOf[vi]!] - solid.vertices[vi]).length;
      minRail = math.min(minRail, len);
    }
    for (final vi in faceVerts) {
      final rail = solid.vertices[otherOf[vi]!] - solid.vertices[vi];
      if ((rail.dot(n) + minRail).abs() > 1e-6) return null;
    }
    return (
      face: face,
      n: n,
      faceVerts: faceVerts,
      otherOf: otherOf,
      minRail: minRail,
    );
  }

  static BrepSolid _shiftVertices(
    BrepSolid solid,
    Set<int> movedVerts,
    Vector3 delta,
  ) {
    // The rail check establishes two parallel caps of a complete prism.
    // Stretch the entire BREP between them: moving only endpoint vertices
    // leaves spline trim curves and their ruled wall surfaces at the old cap.
    final normal = delta.normalized();
    final fixed =
        solid.vertices[List.generate(
          solid.vertices.length,
          (i) => i,
        ).firstWhere((i) => !movedVerts.contains(i))];
    final span = (solid.vertices[movedVerts.first] - fixed).dot(normal);
    final transform = BrepAffineTransform.fromFunction(
      (p) => p + delta * ((p - fixed).dot(normal) / span),
    );
    return BrepXform.materializeForEditing(BrepXform.affine(solid, transform));
  }

  /// Glue a prism onto [faceIndex] along [direction], dropping the original
  /// face. The direct attachment is safe only outside the source solid.
  /// Otherwise union the swept volume so other walls are trimmed at contact.
  static BrepSolid attachFromFace(
    BrepSolid solid,
    int faceIndex,
    Vector3 direction,
  ) {
    if (faceIndex < 0 || faceIndex >= solid.faces.length) {
      throw BrepUnsupported('Attach extrude face is missing');
    }
    final face = solid.faces[faceIndex];
    if (face.loops.isEmpty || direction.length < 1e-12) {
      throw BrepUnsupported(
        'Attach extrude needs a closed face and a direction',
      );
    }
    if (face.loops.length != 1) {
      throw BrepUnsupported(
        'Attach extrude of faces with holes is not implemented yet',
      );
    }
    final surface = face.surface;
    if (surface is BrepPlane &&
        direction.dot(face.outwardNormal(0, 0)) > 0 &&
        !_outsideSupportingPlane(solid, face)) {
      return BrepHeal.run(
        BrepBoolean.run(
          solid,
          fromFace(solid, face, direction),
          BrepBoolean.unionOp,
        ),
      );
    }
    final n = direction.normalized();
    var u = n.cross(Vector3(0, 0, 1));
    if (u.length < 1e-8) u = n.cross(Vector3(0, 1, 0));
    u.normalize();

    final b = BrepBuilder();
    for (final v in solid.vertices) {
      b.vertex(v);
    }
    for (final e in solid.edges) {
      b.edge(e.v0, e.v1, e.curve);
    }
    for (var i = 0; i < solid.faces.length; i++) {
      if (i == faceIndex) continue;
      final f = solid.faces[i];
      b.face(f.surface, f.loops, reversed: f.reversed);
    }

    final loop = face.loops.first;
    final startVis = <int>[];
    for (final c in loop.coedges) {
      final e = solid.edges[c.edge];
      startVis.add(c.reversed ? e.v1 : e.v0);
    }
    final farOf = <int, int>{};
    final vertOf = <int, int>{};
    for (final vi in startVis.toSet()) {
      final far = solid.vertices[vi] + direction;
      farOf[vi] = b.vertex(far);
      vertOf[vi] = b.edge(vi, farOf[vi]!, BrepLine(solid.vertices[vi], far));
    }

    final eFar = <int>[];
    for (var i = 0; i < loop.coedges.length; i++) {
      final c = loop.coedges[i];
      final e = solid.edges[c.edge];
      final vs = c.reversed ? e.v1 : e.v0;
      final ve = c.reversed ? e.v0 : e.v1;
      final a = solid.vertices[vs] + direction;
      final d = solid.vertices[ve] + direction;
      final curve = _curveAlong(e.curve, c.reversed, a, d, direction);
      eFar.add(b.edge(farOf[vs]!, farOf[ve]!, curve));
      if (curve is BrepCircle && curve.normal.cross(n).length < 1e-6) {
        final axis = curve.normal.dot(n) >= 0 ? curve.normal : -curve.normal;
        // Concave arcs face into their cylinder, unlike convex roundovers.
        // A closed edge graph does not detect a wrongly oriented surface.
        final mid = (curve.t0 + curve.t1) / 2;
        final inward =
            curve
                .tangent(mid)
                .cross(direction)
                .dot(curve.evaluate(mid) - curve.center) <
            0;
        b.face(
          BrepCylinder(
            origin: curve.center,
            axis: axis,
            radius: curve.radius,
            xDir: curve.xDir,
          ),
          [
            BrepLoop([
              BrepCoedge(c.edge, reversed: c.reversed),
              BrepCoedge(vertOf[ve]!),
              BrepCoedge(eFar.last, reversed: true),
              BrepCoedge(vertOf[vs]!, reversed: true),
            ]),
          ],
          reversed: inward,
        );
        continue;
      }
      final fn = (solid.vertices[ve] - solid.vertices[vs]).cross(direction);
      if (fn.length < 1e-12) continue;
      b.face(BrepPlane(origin: solid.vertices[vs], normal: fn.normalized()), [
        BrepLoop([
          BrepCoedge(c.edge, reversed: c.reversed),
          BrepCoedge(vertOf[ve]!),
          BrepCoedge(eFar.last, reversed: true),
          BrepCoedge(vertOf[vs]!, reversed: true),
        ]),
      ]);
    }
    b.face(
      BrepPlane(
        origin: solid.vertices[startVis.first] + direction,
        normal: n,
        xDir: u,
      ),
      [
        BrepLoop([for (final e in eFar) BrepCoedge(e)]),
      ],
    );
    return BrepHeal.run(b.build());
  }

  // Prove that an outward sweep cannot meet any existing material. A closed
  // edge graph alone permits self-overlapping shells. Work in the face frame
  // so rotated parts have the same test, and use conservative curved bounds
  // rather than vertices (which can miss the bulge of a curved return wall).
  static bool _outsideSupportingPlane(BrepSolid solid, BrepFace face) {
    final plane = face.surface as BrepPlane;
    final n = face.outwardNormal(0, 0).normalized();
    final x = plane.xDir;
    final y = n.cross(x).normalized();
    final local = BrepXform.map(solid, (p) {
      final d = p - plane.origin;
      return Vector3(d.dot(x), d.dot(y), d.dot(n));
    });
    final bounds = conservativeSolidBounds(local);
    if (bounds == null) return false;
    return bounds.max.z <= math.max(1e-8, bounds.diagonal * 1e-10);
  }

  /// Sweep a B-Rep face along [direction], keeping circular edges as cylinders.
  /// [originShift] moves the profile before the sweep (face-extrude slab).
  static BrepSolid fromFace(
    BrepSolid solid,
    BrepFace face,
    Vector3 direction, {
    Vector3? originShift,
  }) {
    final transformed = BrepXform.unwrapAffine(solid);
    if (transformed != null) {
      final i = solid.faces.indexOf(face), m = transformed.transform;
      if (i < 0) throw BrepUnsupported('Extrude face is missing');
      return BrepXform.affine(
        fromFace(
          transformed.solid,
          transformed.solid.faces[i],
          m.inverseVector(direction),
          originShift: originShift == null
              ? null
              : m.inverseVector(originShift),
        ),
        m,
      );
    }
    final surface = face.surface;
    if (surface is! BrepPlane ||
        face.loops.isEmpty ||
        !direction.x.isFinite ||
        !direction.y.isFinite ||
        !direction.z.isFinite ||
        direction.length < 1e-12) {
      throw BrepUnsupported(
        'Extrude needs a planar face and a finite direction',
      );
    }
    final shift = originShift ?? Vector3.zero();
    var n = surface.normal;
    if (n.dot(direction) < 0) n = -n;
    final height = n.dot(direction);
    if (height < 1e-10 * direction.length)
      throw BrepUnsupported('Extrude direction lies in the face plane');
    final lateral = direction - n * height;
    if (lateral.length > 1e-9 * direction.length) {
      final origin = surface.origin + shift;
      final m = BrepAffineTransform.fromFunction(
        (p) => p + lateral * ((p - origin).dot(n) / height),
      );
      return BrepXform.affine(
        fromFace(solid, face, n * height, originShift: shift),
        m,
      );
    }
    // Saved fillet/hollow trims can have reversed or stale arc parameters.
    // Align them before sampling the profile or constructing swept walls;
    // repairing the Boolean operands later cannot fix a malformed cutter.
    // Face loops and edge indices remain unchanged, so the selection is valid.
    solid = BrepHeal.alignCircularTrims(solid);
    final b = BrepBuilder(), lower = <BrepLoop>[], upper = <BrepLoop>[];
    for (var li = 0; li < face.loops.length; li++) {
      final original = face.loops[li];
      final points = BrepTessellate.sampleLoop(solid, original, 32);
      var area = Vector3.zero();
      for (var i = 0; i < points.length; i++) {
        area += (points[i] - surface.origin).cross(
          points[(i + 1) % points.length] - surface.origin,
        );
      }
      if (area.length < 1e-12)
        throw BrepUnsupported('Extrude boundary has zero area');
      final reverse = (area.dot(n) < 0) != (li > 0);
      final uses = reverse
          ? [
              for (final co in original.coedges.reversed)
                BrepCoedge(co.edge, reversed: !co.reversed),
            ]
          : original.coedges;
      final curves = <BrepCurve>[];
      for (final co in uses) {
        final edge = solid.edges[co.edge];
        curves.add(
          _curveAlong(
            edge.curve,
            co.reversed,
            solid.vertices[co.reversed ? edge.v1 : edge.v0] + shift,
            solid.vertices[co.reversed ? edge.v0 : edge.v1] + shift,
            shift,
          ),
        );
      }
      final count = curves.length;
      final v0 = [for (final c in curves) b.vertex(c.evaluate(c.t0))];
      final v1 = [
        for (final c in curves) b.vertex(c.evaluate(c.t0) + direction),
      ];
      final e0 = [
        for (var i = 0; i < count; i++)
          b.edge(v0[i], v0[(i + 1) % count], curves[i]),
      ];
      final e1 = [
        for (var i = 0; i < count; i++)
          b.edge(
            v1[i],
            v1[(i + 1) % count],
            _translatedCurve(curves[i], direction),
          ),
      ];
      final rails = [
        for (var i = 0; i < count; i++)
          b.edge(v0[i], v1[i], BrepLine(b.vertices[v0[i]], b.vertices[v1[i]])),
      ];
      lower.add(
        BrepLoop([for (final e in e0.reversed) BrepCoedge(e, reversed: true)]),
      );
      upper.add(BrepLoop([for (final e in e1) BrepCoedge(e)]));
      final ruled = curves.any((c) => c is! BrepLine && c is! BrepCircle);
      for (var i = 0; i < count; i++) {
        final curve = curves[i], topCurve = b.edges[e1[i]].curve;
        final loop = BrepLoop([
          BrepCoedge(e0[i]),
          BrepCoedge(rails[(i + 1) % count]),
          BrepCoedge(e1[i], reversed: true),
          BrepCoedge(rails[i], reversed: true),
        ]);
        if (ruled && curve is! BrepLine) {
          b.face(BrepRuledSurface(curve, topCurve), [loop]);
        } else if (curve is BrepCircle) {
          final mid = (curve.t0 + curve.t1) / 2;
          final inward =
              curve
                  .tangent(mid)
                  .cross(direction)
                  .dot(curve.evaluate(mid) - curve.center) <
              0;
          b.face(
            BrepCylinder(
              origin: curve.center,
              axis: n,
              radius: curve.radius,
              xDir: curve.xDir,
            ),
            curve.isClosed
                ? [
                    BrepLoop([BrepCoedge(e0[i])]),
                    BrepLoop([BrepCoedge(e1[i], reversed: true)]),
                  ]
                : [loop],
            reversed: inward,
          );
        } else {
          b.face(
            BrepPlane(
              origin: curve.evaluate(curve.t0),
              normal: curve.tangent(curve.t0).cross(direction),
            ),
            [loop],
          );
        }
      }
    }
    b.face(
      BrepPlane(origin: surface.origin + shift, normal: -n, xDir: surface.xDir),
      lower,
    );
    b.face(
      BrepPlane(
        origin: surface.origin + shift + direction,
        normal: n,
        xDir: surface.xDir,
      ),
      upper,
    );
    return pairPlanarLineContacts(BrepHeal.sewClosedRims(b.build()));
  }

  /// Sweep mixed line/arc loops. Concave arcs (a D-notch) reverse the
  /// cylinder so the wall faces out of the solid.
  static BrepSolid fromCurves(
    List<Vector3> pts,
    List<BrepCurve> curves,
    Vector3 direction,
  ) {
    if (pts.isEmpty || pts.length != curves.length)
      throw BrepUnsupported('Extrude needs a closed curve loop');
    final b = BrepBuilder();
    for (final p in pts) {
      b.vertex(p);
    }
    final edges = [
      for (var i = 0; i < pts.length; i++)
        b.edge(i, (i + 1) % pts.length, curves[i]),
    ];
    final samples = [
      for (final c in curves)
        for (var j = 0; j < 4; j++) c.evaluate(c.t0 + (c.t1 - c.t0) * j / 4),
    ];
    final n = _profileNormal(samples, direction);
    b.face(BrepPlane(origin: pts.first, normal: n), [
      BrepLoop([for (final e in edges) BrepCoedge(e)]),
    ]);
    final source = b.build();
    return fromFace(source, source.faces.single, direction);
  }

  /// Circular runs of 6+ samples become one [BrepCircle]; the rest stay lines.
  ///
  /// [linearEdges] protects polygon sides and line segments from arc fitting.
  /// Without source metadata, six or more cocircular points are ambiguous.
  static ({List<Vector3> pts, List<BrepCurve> curves})? fitCircularRuns(
    List<Vector3> ring,
    Vector3 axis, {
    List<(Vector3, Vector3)> linearEdges = const [],
  }) {
    if (ring.length < 6) return null;
    final n = ring.length;
    final ax = axis.normalized();
    var start = 0;
    var bestTurn = -1.0;
    for (var i = 0; i < n; i++) {
      final prev = ring[(i - 1 + n) % n];
      final cur = ring[i];
      final next = ring[(i + 1) % n];
      final e1 = cur - prev;
      final e2 = next - cur;
      final l1 = e1.length;
      final l2 = e2.length;
      if (l1 < 1e-12 || l2 < 1e-12) continue;
      final turn = e1.cross(e2).length / (l1 * l2);
      if (turn > bestTurn) {
        bestTurn = turn;
        start = i;
      }
    }
    final pts = [for (var i = 0; i < n; i++) ring[(start + i) % n]];
    final closedPts = [...pts, pts.first];

    final straight = [
      for (var i = 0; i < n; i++)
        _isLinearEdge(pts[i], pts[(i + 1) % n], linearEdges),
    ];
    final outPts = <Vector3>[];
    final outCurves = <BrepCurve>[];
    var foundArc = false;
    var i = 0;
    while (i < n) {
      ({
        Vector3 center,
        double radius,
        Vector3 xDir,
        double span,
        Vector3 normal,
      })?
      best;
      var runLen = 0;
      // The closing edge belongs to the ring too. An arc ending at the first
      // vertex must include that endpoint instead of leaving a short chord.
      var maxRun = n - i + 1;
      for (var k = i; k < n; k++) {
        if (straight[k]) {
          maxRun = k - i + 1;
          break;
        }
      }
      for (var len = maxRun; len >= 6; len--) {
        // All vertices of a regular polygon also lie on a circle. Only fit
        // a run spanning the ring when source metadata establishes its chord.
        if (len >= n && !straight.any((edge) => edge)) continue;
        final fit = _fitArc(closedPts.sublist(i, i + len), ax);
        if (fit != null) {
          best = fit;
          runLen = len;
          break;
        }
      }
      outPts.add(pts[i]);
      if (best != null) {
        foundArc = true;
        outCurves.add(
          BrepCircle(
            center: best.center,
            normal: best.normal,
            radius: best.radius,
            xDir: best.xDir,
            t0: 0,
            t1: best.span,
          ),
        );
        i += runLen - 1;
      } else {
        outCurves.add(BrepLine(pts[i], pts[(i + 1) % n]));
        i += 1;
      }
    }
    if (!foundArc || outPts.length < 2 || outCurves.length != outPts.length) {
      return null;
    }
    return (pts: outPts, curves: outCurves);
  }

  static ({
    Vector3 center,
    double radius,
    Vector3 xDir,
    double span,
    Vector3 normal,
  })?
  _fitArc(List<Vector3> run, Vector3 axis) {
    if (run.length < 6) return null;
    final a = run.first;
    final b = run[run.length ~/ 2];
    final c = run.last;
    var xAxis = brepFrameU(axis);
    final yAxis = axis.cross(xAxis);
    Vector2 to2(Vector3 p) {
      final d = p - a;
      return Vector2(d.dot(xAxis), d.dot(yAxis));
    }

    for (final p in run) {
      if ((p - a).dot(axis).abs() > 1e-4) return null;
    }
    final circ = _circumcircle(to2(a), to2(b), to2(c));
    if (circ == null) return null;
    final r = circ.$2;
    if (r < 1e-6) return null;
    final tol = math.max(1e-4, r * 1e-3);
    final center = a + xAxis * circ.$1.x + yAxis * circ.$1.y;
    var maxSag = 0.0;
    for (final p in run) {
      final radial = p - center;
      final along = radial.dot(axis);
      if (along.abs() > 1e-4) return null;
      if (((radial - axis * along).length - r).abs() > tol) return null;
    }
    final chord = c - a;
    final chordLen = chord.length;
    if (chordLen < 1e-9) return null;
    final chordN = chord / chordLen;
    for (var k = 1; k < run.length - 1; k++) {
      final d = run[k] - a;
      final sag = (d - chordN * d.dot(chordN)).length;
      if (sag > maxSag) maxSag = sag;
    }
    if (maxSag < math.max(1e-4, r * 0.02)) return null;

    var radial0 = a - center;
    radial0 -= axis * radial0.dot(axis);
    if (radial0.length < 1e-12) return null;
    radial0.normalize();
    final y = axis.cross(radial0);
    var span = 0.0;
    var prev = 0.0;
    for (var k = 1; k < run.length; k++) {
      var rad = run[k] - center;
      rad -= axis * rad.dot(axis);
      if (rad.length < 1e-12) return null;
      var t = math.atan2(rad.dot(y), rad.dot(radial0));
      var d = t - prev;
      while (d > math.pi) {
        d -= math.pi * 2;
      }
      while (d < -math.pi) {
        d += math.pi * 2;
      }
      if (k > 1 && span.abs() > 1e-8 && d.sign != span.sign && d.abs() > 1e-4) {
        return null;
      }
      span += d;
      prev = t;
    }
    if (span.abs() < 15 * math.pi / 180) return null;
    if (span.abs() > math.pi * 2 - 1e-3) return null;

    var normal = axis.clone();
    var spanOut = span;
    if (spanOut < 0) {
      normal = -normal;
      spanOut = -spanOut;
    }
    var xDir = a - center;
    xDir -= normal * xDir.dot(normal);
    if (xDir.length < 1e-12) return null;
    xDir.normalize();
    return (
      center: center,
      radius: r,
      xDir: xDir,
      span: spanOut,
      normal: normal,
    );
  }

  static (Vector2, double)? _circumcircle(Vector2 a, Vector2 b, Vector2 c) {
    final d = 2 * (a.x * (b.y - c.y) + b.x * (c.y - a.y) + c.x * (a.y - b.y));
    if (d.abs() < 1e-18) return null;
    final a2 = a.dot(a);
    final b2 = b.dot(b);
    final c2 = c.dot(c);
    final ux = (a2 * (b.y - c.y) + b2 * (c.y - a.y) + c2 * (a.y - b.y)) / d;
    final uy = (a2 * (c.x - b.x) + b2 * (a.x - c.x) + c2 * (b.x - a.x)) / d;
    final center = Vector2(ux, uy);
    final r = (a - center).length;
    if (r < 1e-12) return null;
    return (center, r);
  }

  static BrepCurve _curveAlong(
    BrepCurve curve,
    bool reversed,
    Vector3 start,
    Vector3 end,
    Vector3 shift,
  ) {
    if (curve is BrepAffineCurve) {
      // Orient in the source curve's coordinates, then apply the exact map.
      // start/end already include the extrusion's origin shift; remove it
      // before mapping back so nested translations are not applied twice.
      final oriented = _curveAlong(
        curve.base,
        reversed,
        curve.transform.inversePoint(start - shift),
        curve.transform.inversePoint(end - shift),
        Vector3.zero(),
      );
      if (oriented is BrepLine) return BrepLine(start, end);
      final m = curve.transform;
      return BrepAffineCurve(
        oriented,
        BrepAffineTransform(m.x, m.y, m.z, m.offset + shift),
      );
    }
    BrepCurve? oriented;
    if (curve is BrepSurfaceIntersectionCurve)
      oriented = curve.trimmed(
        reversed ? curve.t1 : curve.t0,
        reversed ? curve.t0 : curve.t1,
      );
    if (curve is BrepRuledQuadricCurve)
      oriented = curve.trimmed(
        reversed ? curve.t1 : curve.t0,
        reversed ? curve.t0 : curve.t1,
      );
    if (curve is BrepPlaneSectionCurve)
      oriented = curve.trimmed(
        reversed ? curve.t1 : curve.t0,
        reversed ? curve.t0 : curve.t1,
      );
    if (curve is BrepSplineCurve)
      oriented = curve.trimmed(
        reversed ? curve.t1 : curve.t0,
        reversed ? curve.t0 : curve.t1,
      );
    if (oriented != null) return _translatedCurve(oriented, shift);
    if (curve is BrepLine) return BrepLine(start, end);
    if (curve is BrepEllipse) {
      return BrepEllipse(
        center: curve.center + shift,
        normal: reversed ? -curve.normal : curve.normal,
        xDir: curve.xDir,
        rMajor: curve.rMajor,
        rMinor: curve.rMinor,
        t0: reversed ? -curve.t1 : curve.t0,
        t1: reversed ? -curve.t0 : curve.t1,
      );
    }
    if (curve is BrepCircle) {
      // Preserve the signed parameter sweep rather than inferring an arc
      // from its endpoints. A semicircle has two equally short solutions;
      // atan2 rounding used to select the opposite half after reversal.
      final startT = reversed ? curve.t1 : curve.t0;
      final sweep = (reversed ? -1 : 1) * (curve.t1 - curve.t0);
      return BrepCircle(
        center: curve.center + shift,
        normal: sweep >= 0 ? curve.normal : -curve.normal,
        radius: curve.radius,
        xDir: curve.evaluate(startT) - curve.center,
        t0: 0,
        t1: sweep.abs(),
      );
    }
    throw BrepUnsupported(
      'Extrude of ${curve.runtimeType} edges is not implemented yet',
    );
  }

  static BrepCurve _translatedCurve(BrepCurve curve, Vector3 d) {
    if (curve is BrepSurfaceIntersectionCurve ||
        curve is BrepRuledQuadricCurve ||
        curve is BrepPlaneSectionCurve ||
        curve is BrepSplineCurve ||
        curve is BrepAffineCurve ||
        curve is BrepIntersectionCurve) {
      return BrepAffineCurve(
        curve,
        BrepAffineTransform(
          Vector3(1, 0, 0),
          Vector3(0, 1, 0),
          Vector3(0, 0, 1),
          d,
        ),
      );
    }
    if (curve is BrepEllipse)
      return BrepEllipse(
        center: curve.center + d,
        normal: curve.normal,
        xDir: curve.xDir,
        rMajor: curve.rMajor,
        rMinor: curve.rMinor,
        t0: curve.t0,
        t1: curve.t1,
      );
    if (curve is BrepLine) return BrepLine(curve.a + d, curve.b + d);
    if (curve is BrepCircle) {
      return BrepCircle(
        center: curve.center + d,
        normal: curve.normal,
        radius: curve.radius,
        xDir: curve.xDir,
        t0: curve.t0,
        t1: curve.t1,
      );
    }
    throw BrepUnsupported(
      'Extrude of ${curve.runtimeType} edges is not implemented yet',
    );
  }

  static Vector3 _profileNormal(List<Vector3> ring, Vector3 direction) {
    if (ring.length < 3 ||
        [
          direction,
          ...ring,
        ].any((p) => !p.x.isFinite || !p.y.isFinite || !p.z.isFinite) ||
        direction.length < 1e-12) {
      throw BrepUnsupported('Extrude needs a finite profile and direction');
    }
    var n = Vector3.zero();
    final origin = ring.first;
    for (var i = 0; i < ring.length; i++) {
      n += (ring[i] - origin).cross(ring[(i + 1) % ring.length] - origin);
    }
    if (n.length < 1e-12)
      throw BrepUnsupported('Extrude profile has zero area');
    n.normalize();
    final extent = ring.map((p) => (p - origin).length).reduce(math.max);
    if (ring.any(
      (p) => (p - origin).dot(n).abs() > math.max(1e-8, extent * 1e-9),
    )) {
      throw BrepUnsupported('Extrude profile must be planar');
    }
    if (n.dot(direction).abs() < 1e-10 * direction.length)
      throw BrepUnsupported('Extrude direction lies in the profile plane');
    return n.dot(direction) < 0 ? -n : n;
  }

  /// Drop collinear vertices and wind [ring] CCW about [direction].
  static List<Vector3> prepareRing(List<Vector3> ring, Vector3 direction) {
    _profileNormal(ring, direction);
    var pts = collapseCollinear(ring);
    if (pts.length < 3) pts = [for (final p in ring) p.clone()];
    return _ensureWinding(pts, direction);
  }

  /// Remove vertices that sit on a straight run of edges.
  static List<Vector3> collapseCollinear(List<Vector3> ring) {
    if (ring.length < 4) {
      return [for (final p in ring) p.clone()];
    }
    var pts = [for (final p in ring) p.clone()];
    for (var pass = 0; pass < 8; pass++) {
      if (pts.length < 4) break;
      final out = <Vector3>[];
      for (var i = 0; i < pts.length; i++) {
        final prev = pts[(i - 1 + pts.length) % pts.length];
        final cur = pts[i];
        final next = pts[(i + 1) % pts.length];
        if (_keepVertex(prev, cur, next)) out.add(cur);
      }
      if (out.length == pts.length || out.length < 3) break;
      pts = out;
    }
    return pts.length >= 3 ? pts : [for (final p in ring) p.clone()];
  }

  static bool _keepVertex(Vector3 prev, Vector3 cur, Vector3 next) {
    final e1 = cur - prev;
    final e2 = next - cur;
    final l1 = e1.length;
    final l2 = e2.length;
    if (l1 < 1e-9 || l2 < 1e-9) return false;
    final cr = e1.cross(e2).length / (l1 * l2);
    final dot = e1.dot(e2) / (l1 * l2);
    return cr > 1e-6 || dot < 0.999999;
  }

  static List<Vector3> _ensureWinding(List<Vector3> ring, Vector3 direction) {
    var area = Vector3.zero();
    for (var i = 0; i < ring.length; i++) {
      area += ring[i].cross(ring[(i + 1) % ring.length]);
    }
    if (area.dot(direction) < 0) {
      return ring.reversed.map((p) => p.clone()).toList();
    }
    return ring;
  }

  /// Outer boundary of a triangulated planar profile, as a closed ring.
  static List<Vector3> outerRing(
    List<Vector3> vertices,
    List<List<int>> triangles,
  ) {
    String ek(int a, int b) => a < b ? '$a-$b' : '$b-$a';
    final count = <String, int>{};
    final directed = <(int, int)>[];
    for (final t in triangles) {
      if (t.length < 3) continue;
      for (var i = 1; i < t.length - 1; i++) {
        final tri = [t[0], t[i], t[i + 1]];
        for (var e = 0; e < 3; e++) {
          final a = tri[e];
          final b = tri[(e + 1) % 3];
          final k = ek(a, b);
          count[k] = (count[k] ?? 0) + 1;
          directed.add((a, b));
        }
      }
    }
    final nexts = <int, List<int>>{};
    for (final (a, b) in directed) {
      if ((count[ek(a, b)] ?? 0) != 1) continue;
      (nexts[a] ??= []).add(b);
    }
    if (nexts.isEmpty) {
      throw BrepUnsupported('Extrude profile has no outer boundary');
    }
    final start = nexts.keys.first;
    final ringIdx = <int>[start];
    final used = <String>{};
    var cur = start;
    while (ringIdx.length <= vertices.length) {
      final opts = nexts[cur];
      if (opts == null || opts.isEmpty) break;
      int? nxt;
      for (final cand in opts) {
        final k = '$cur-$cand';
        if (used.add(k)) {
          nxt = cand;
          break;
        }
      }
      if (nxt == null) break;
      if (nxt == start) break;
      ringIdx.add(nxt);
      cur = nxt;
    }
    if (ringIdx.length < 3) {
      throw BrepUnsupported('Extrude profile boundary is degenerate');
    }
    return [for (final i in ringIdx) vertices[i].clone()];
  }
}
