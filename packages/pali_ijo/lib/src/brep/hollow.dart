import 'dart:math' as math;

import 'package:vector_math/vector_math_64.dart';

import '../mesh/mesh_data.dart';
import 'fillet.dart';
import 'convex_planes.dart';
import 'curve_equivalence.dart';
import 'geom.dart';
import 'heal.dart';
import 'intersect.dart';
import 'primitives.dart';
import 'provenance.dart';
import 'solid.dart';
import 'spline_convert.dart';
import 'slice.dart';
import 'tessellate.dart';
import 'xform.dart';

part 'hollow_ruled_offset.dart';

/// Face-offset shell: inset every face except the opening, drop the opening,
/// and stitch the outer rim to the inner rim.
class BrepHollow {
  BrepHollow._();

  /// Shells [solid] to [wallThicknessMm], leaving [openingFaces] open.
  ///
  /// The outer faces keep their tags. The inner wall does not inherit: it is a
  /// face offset inward, so it lies on a *different* plane or a smaller-radius
  /// cylinder than the face it was derived from, and comes back tagged as
  /// created by the hollow. That is the honest answer — the inside of a shell
  /// is new surface, not a re-trim of the outside.
  static BrepSolid run(
    BrepSolid solid,
    List<int> openingFaces,
    Vector3 openingNormal, {
    required double wallThicknessMm,
  }) => BrepProvenance.inheritFrom(
    _run(solid, openingFaces, openingNormal, wallThicknessMm: wallThicknessMm),
    solid,
    op: 'hollow',
  );

  /// Shell a picked display region while retaining analytic geometry.
  /// A whole sphere has one face, so removing that face would remove the
  /// entire body. Give the picked pole a circular opening into the cavity.
  static BrepSolid fromMeshOpening(
    BrepSolid solid,
    MeshData mesh,
    List<int> meshFaces,
    Vector3 normal, {
    required double wallThicknessMm,
  }) {
    if (wallThicknessMm.isFinite && wallThicknessMm.abs() < 1e-9) {
      return solid.clone();
    }
    if (solid.faces.length == 1 &&
        solid.faces.single.surface is BrepSphere &&
        solid.faces.single.loops.isEmpty &&
        !solid.faces.single.reversed) {
      final sphere = solid.faces.single.surface as BrepSphere;
      final wall = wallThicknessMm.abs();
      if (!wall.isFinite || wall >= sphere.radius || normal.length2 < 1e-20) {
        throw BrepUnsupported(
          'Hollow thickness must be smaller than the sphere radius',
        );
      }
      final n = normal.normalized();
      var height =
          sphere.radius - wall - math.max(wall * .25, sphere.radius * .005);
      for (final fi in meshFaces) {
        if (fi < 0 || fi >= mesh.faces.length) continue;
        for (final vi in mesh.faces[fi]) {
          height = math.min(height, (mesh.vertices[vi] - sphere.center).dot(n));
        }
      }
      final parts = BrepSlice.parts(solid, sphere.center + n * height, n);
      for (final part in parts) {
        for (var fi = 0; fi < part.faces.length; fi++) {
          final face = part.faces[fi];
          if (face.surface is BrepPlane &&
              face.outwardNormal(0, 0).dot(n) > .999) {
            return run(part, [fi], n, wallThicknessMm: wallThicknessMm);
          }
        }
      }
      throw BrepUnsupported('Hollow could not create the spherical opening');
    }
    return run(
      solid,
      facesFromMeshOpening(solid, mesh, meshFaces, normal),
      normal,
      wallThicknessMm: wallThicknessMm,
    );
  }

  /// Smallest convex cylinder / sphere / torus-minor radius that an inward
  /// offset would shrink. Null when the solid has no such rounds.
  ///
  /// Wall thickness must stay below this to retain those rounded inner walls.
  /// Offsets that collapse a round reject instead of substituting a mesh pocket.
  static double? minConvexOffsetRadius(BrepSolid solid) {
    double? minR;
    void consider(double r) {
      if (r <= 1e-6) return;
      minR = minR == null ? r : math.min(minR!, r);
    }

    for (final face in solid.faces) {
      final s = face.surface;
      if (s is BrepCylinder) {
        final sn = s.normalAt(0, 0);
        final out = face.outwardNormal(0, 0);
        if (out.dot(sn) > 0) consider(s.radius);
      } else if (s is BrepSphere) {
        if (!face.reversed) consider(s.radius);
      } else if (s is BrepTorus) {
        if (!face.reversed) consider(s.minorRadius);
      }
    }
    return minR;
  }

  static BrepSolid _run(
    BrepSolid solid,
    List<int> openingFaces,
    Vector3 openingNormal, {
    required double wallThicknessMm,
  }) {
    if (!wallThicknessMm.isFinite) {
      throw BrepUnsupported('Hollow thickness must be finite');
    }
    final t = wallThicknessMm.abs();
    if (t < 1e-9) return solid.clone();
    for (final fi in openingFaces) {
      if (fi < 0 || fi >= solid.faces.length) {
        throw BrepUnsupported('Hollow opening face is out of range');
      }
    }
    // Placed bodies may wrap even their planes in affine surfaces. Recover
    // those exact planes before intersecting offsets; never modify the source.
    BrepPlane? planeOf(BrepSurface surface) {
      if (surface is BrepPlane) return surface;
      if (surface is BrepAffineSurface) {
        final base = planeOf(surface.base);
        if (base != null) {
          return BrepPlane(
            origin: surface.transform.point(base.origin),
            normal: surface.transform.normal(base.normal),
            xDir: surface.transform.vector(base.xDir).normalized(),
          );
        }
      }
      return null;
    }

    solid = BrepXform.materializeForEditing(
      solid,
      preserveAffineSweeps: true,
    ).clone();
    bool lineCurve(BrepCurve c) =>
        c is BrepLine || c is BrepAffineCurve && lineCurve(c.base);
    for (final edge in solid.edges) {
      final c = edge.curve;
      if (c is BrepAffineCurve && lineCurve(c)) {
        edge.curve = BrepLine(c.evaluate(c.t0), c.evaluate(c.t1));
      }
    }
    for (final face in solid.faces) {
      face.surface = planeOf(face.surface) ?? face.surface;
    }
    final opening = openingFaces.toSet();
    if (BrepConvexPlanes.supports(solid)) {
      final clipped = BrepConvexPlanes.clip(solid, [
        for (var fi = 0; fi < solid.faces.length; fi++)
          BrepPlane(
            origin:
                (solid.faces[fi].surface as BrepPlane).origin -
                solid.faces[fi].outwardNormal(0, 0) *
                    (opening.contains(fi) ? 0 : t),
            normal: solid.faces[fi].outwardNormal(0, 0),
          ),
      ], sources: List.generate(solid.faces.length, (i) => i));
      if (opening.any((fi) => !clipped.faceMap.containsKey(fi))) {
        throw BrepUnsupported('Hollow thickness closes the selected opening');
      }
      return BrepHeal.sewClosedRims(
        _stitch(
          solid,
          clipped.solid,
          opening,
          innerFaceMap: clipped.faceMap,
          planarRims: true,
        ),
      );
    }
    final inner = _offsetInward(solid, t, opening);
    return BrepHeal.sewClosedRims(_stitch(solid, inner, opening));
  }

  static BrepSolid _stitch(
    BrepSolid outer,
    BrepSolid inner,
    Set<int> opening, {
    Map<int, int>? innerFaceMap,
    bool planarRims = false,
  }) {
    final b = BrepBuilder();
    int weld(Vector3 p) {
      for (var i = 0; i < b.vertices.length; i++) {
        if ((b.vertices[i] - p).length < 1e-6) return i;
      }
      return b.vertex(p);
    }

    final edgeOf = <String, int>{};
    (int, bool) sew(int i, int j, BrepCurve curve) {
      // Legacy planar booleans can store a line opposite to its vertex order.
      // New rim topology (and STEP edge directions) must agree with the wire.
      if (planarRims && curve is BrepLine) {
        curve = BrepLine(b.vertices[i], b.vertices[j]);
      }
      final mid = curve.evaluate(0.5 * (curve.t0 + curve.t1));
      final q =
          '${(mid.x * 1e4).round()}_${(mid.y * 1e4).round()}_${(mid.z * 1e4).round()}';
      final lo = i <= j ? i : j;
      final hi = i <= j ? j : i;
      final k = '$lo-$hi-$q';
      final existing = edgeOf[k];
      if (existing != null) {
        final e = b.edges[existing];
        if (i == j) {
          final c = e.curve;
          final oldTangent = c.tangent((c.t0 + c.t1) * .5) * (c.t1 - c.t0).sign;
          final newTangent =
              curve.tangent((curve.t0 + curve.t1) * .5) *
              (curve.t1 - curve.t0).sign;
          return (existing, oldTangent.dot(newTangent) < 0);
        }
        return (existing, e.v0 != i);
      }
      final id = b.edge(i, j, curve);
      edgeOf[k] = id;
      return (id, false);
    }

    void addSolid(BrepSolid s, {required bool reverse}) {
      final omitted = reverse && innerFaceMap != null
          ? opening.map((fi) => innerFaceMap[fi]).toSet()
          : opening;
      for (var fi = 0; fi < s.faces.length; fi++) {
        if (omitted.contains(fi)) continue;
        final f = s.faces[fi];
        // Closed-circle sphere wires select the retained cap independently
        // of the face normal (also used by classification and saved files).
        // Reversing the inner shell must not select the opposite cap.
        final reverseWire =
            reverse &&
            !(f.surface is BrepSphere &&
                f.loops.every(
                  (loop) =>
                      loop.coedges.length == 1 &&
                      s.edges[loop.coedges.single.edge].curve is BrepCircle &&
                      (s.edges[loop.coedges.single.edge].curve as BrepCircle)
                          .isClosed,
                ));
        final loops = <BrepLoop>[];
        for (final loop in f.loops) {
          final cs = <BrepCoedge>[];
          for (final c in loop.coedges) {
            final e = s.edges[c.edge];
            final v0 = weld(s.vertices[e.v0]);
            final v1 = weld(s.vertices[e.v1]);
            final se = sew(v0, v1, e.curve);
            cs.add(
              BrepCoedge(se.$1, reversed: se.$2 ^ c.reversed ^ reverseWire),
            );
          }
          loops.add(BrepLoop(reverseWire ? cs.reversed.toList() : cs));
        }
        b.face(f.surface, loops, reversed: reverse ? !f.reversed : f.reversed);
      }
    }

    addSolid(outer, reverse: false);
    addSolid(inner, reverse: true);

    for (final fi in opening) {
      final innerFi = innerFaceMap?[fi] ?? fi;
      if (fi < 0 || fi >= outer.faces.length || innerFi >= inner.faces.length) {
        continue;
      }
      final of = outer.faces[fi];
      final inf = inner.faces[innerFi];
      if (planarRims) {
        // Offset corners may split or disappear, so the inner opening need
        // not have the same number of edges. One planar face with an inner
        // trim joins both wires without inventing a correspondence.
        BrepLoop copyLoop(BrepSolid s, BrepLoop loop, bool reverse) {
          final coedges = <BrepCoedge>[];
          for (final c in loop.coedges) {
            final e = s.edges[c.edge];
            final se = sew(
              weld(s.vertices[e.v0]),
              weld(s.vertices[e.v1]),
              e.curve,
            );
            coedges.add(
              BrepCoedge(se.$1, reversed: se.$2 ^ c.reversed ^ reverse),
            );
          }
          return BrepLoop(reverse ? coedges.reversed.toList() : coedges);
        }

        b.face(of.surface, [
          copyLoop(outer, of.loops.single, false),
          copyLoop(inner, inf.loops.single, true),
        ], reversed: of.reversed);
        continue;
      }
      final loopCount = math.min(of.loops.length, inf.loops.length);
      if (loopCount == 0) continue;
      // Both loop conventions occur in saved and Boolean-created faces:
      // some follow the surface normal, others the outward face normal.
      // Determine the planar rim orientation from the actual outer wire.
      var reverseRim = of.reversed;
      if (of.surface is BrepPlane) {
        final points = BrepTessellate.sampleLoop(outer, of.loops.first, 32);
        final origin = (of.surface as BrepPlane).origin;
        var area = Vector3.zero();
        for (var i = 0; i < points.length; i++) {
          area += (points[i] - origin).cross(
            points[(i + 1) % points.length] - origin,
          );
        }
        reverseRim = area.dot(of.outwardNormal(0, 0)) < 0;
      }
      for (var li = 0; li < loopCount; li++) {
        _stitchOpeningLoop(
          b,
          weld,
          sew,
          outer,
          inner,
          of.loops[li],
          inf.loops[li],
          of.outwardNormal(0, 0),
          reverseRim,
        );
      }
    }
    return b.build();
  }

  static void _stitchOpeningLoop(
    BrepBuilder b,
    int Function(Vector3) weld,
    (int, bool) Function(int, int, BrepCurve) sew,
    BrepSolid outer,
    BrepSolid inner,
    BrepLoop oLoop,
    BrepLoop iLoop,
    Vector3 openingNormal,
    bool openingReversed,
  ) {
    final oSegs = _loopSegments(outer, oLoop);
    var iSegs = _loopSegments(inner, iLoop);
    if (oSegs.isEmpty || iSegs.isEmpty) return;
    if (_isClosedRing(oSegs) && _isClosedRing(iSegs)) {
      _stitchAnnulus(b, weld, sew, oSegs.first, iSegs.first, openingNormal);
      return;
    }
    if (oSegs.length != iSegs.length) {
      throw BrepUnsupported('Hollow opening loops do not match after offset');
    }
    iSegs = _alignLoop(oSegs, iSegs);

    for (var k = 0; k < oSegs.length; k++) {
      final o = oSegs[k];
      final i = iSegs[k];
      final va = weld(o.$1);
      final vb = weld(o.$2);
      final vib = weld(i.$2);
      final via = weld(i.$1);
      if ({va, vb, vib, via}.length < 3) continue;
      final e0 = sew(va, vb, o.$3);
      final e1 = sew(vb, vib, BrepLine(o.$2, i.$2));
      final e2 = sew(vib, via, _flipped(i.$3));
      final e3 = sew(via, va, BrepLine(i.$1, o.$1));
      var n = (o.$2 - o.$1).cross(i.$1 - o.$1);
      if (n.length < 1e-12) n = (o.$2 - o.$1).cross(i.$2 - o.$1);
      if (n.length < 1e-12) continue;
      // The strip inherits the opening loop's winding; reverse it when
      // that wire runs against the outward orientation of the opening.
      final plane = BrepPlane(origin: o.$1, normal: n.normalized());
      bool onPlane(BrepCurve curve) {
        bool rational(BrepCurve c) =>
            c is BrepLine ||
            c is BrepCircle ||
            c is BrepEllipse ||
            c is BrepSplineCurve ||
            c is BrepAffineCurve && rational(c.base);
        if (!rational(curve)) return false;
        // Control-net containment proves planarity between samples as well.
        return BrepSplineConvert.curve(
          curve,
        ).points.every((p) => plane.height(p).abs() < 1e-7);
      }

      final planar = onPlane(o.$3) && onPlane(i.$3);
      b.face(planar ? plane : BrepRuledSurface(o.$3, i.$3), [
        BrepLoop([
          BrepCoedge(e0.$1, reversed: e0.$2),
          BrepCoedge(e1.$1, reversed: e1.$2),
          BrepCoedge(e2.$1, reversed: e2.$2),
          BrepCoedge(e3.$1, reversed: e3.$2),
        ]),
      ], reversed: openingReversed);
    }
  }

  static bool _isClosedRing(List<(Vector3, Vector3, BrepCurve)> segs) {
    if (segs.length != 1) return false;
    return (segs.first.$1 - segs.first.$2).length < 1e-8;
  }

  static void _stitchAnnulus(
    BrepBuilder b,
    int Function(Vector3) weld,
    (int, bool) Function(int, int, BrepCurve) sew,
    (Vector3, Vector3, BrepCurve) outer,
    (Vector3, Vector3, BrepCurve) inner,
    Vector3 openingNormal,
  ) {
    final oc = outer.$3, ic = inner.$3;
    if (oc is BrepCircle && ic is BrepCircle) {
      final axis = openingNormal.normalized();
      final delta = oc.center - ic.center;
      final height = delta.dot(axis);
      if (height.abs() > 1e-7) {
        // A cap tangent to a fillet moves inward with the offset fillet.
        // Its equal-radius rims need a cylindrical lip, not a flat annulus.
        if ((delta - axis * height).length > 1e-6 ||
            (oc.radius - ic.radius).abs() > 1e-6 ||
            oc.normal.dot(axis).abs() < 1 - 1e-8 ||
            ic.normal.dot(axis).abs() < 1 - 1e-8) {
          throw BrepUnsupported('Hollow needs a nonplanar opening transition');
        }
        final vo = weld(outer.$1), vi = weld(inner.$1);
        final eo = sew(vo, vo, oc), ei = sew(vi, vi, ic);
        b.face(
          BrepCylinder(
            origin: ic.center,
            axis: axis,
            radius: ic.radius,
            xDir: ic.xDir,
          ),
          [
            BrepLoop([BrepCoedge(eo.$1, reversed: eo.$2)]),
            BrepLoop([BrepCoedge(ei.$1, reversed: !ei.$2)]),
          ],
          reversed: height > 0,
        );
        return;
      }
    }
    // Offsetting a through-hole grows its radius. Its inner shell rim is
    // therefore the outer boundary of this annulus, not the hole. Planar
    // tessellation requires the enclosing loop first regardless of which
    // shell supplied it. The rim faces out through the selected opening.
    double ringSize(BrepCurve curve) => switch (curve) {
      BrepCircle c => c.radius * c.radius,
      BrepEllipse e => e.rMajor * e.rMinor,
      _ => _ringArea(curve),
    };
    if (ringSize(inner.$3) > ringSize(outer.$3)) {
      final swap = outer;
      outer = inner;
      inner = swap;
    }
    final va = weld(outer.$1);
    final via = weld(inner.$1);
    if (va == via) return;
    final e0 = sew(va, va, outer.$3);
    final e1 = sew(via, via, inner.$3);
    b.face(BrepPlane(origin: outer.$1, normal: openingNormal.normalized()), [
      BrepLoop([BrepCoedge(e0.$1, reversed: e0.$2)]),
      BrepLoop([BrepCoedge(e1.$1, reversed: !e1.$2)]),
    ]);
  }

  static double _ringArea(BrepCurve curve) {
    final origin = curve.evaluate(curve.t0);
    var area = Vector3.zero(), previous = Vector3.zero();
    for (var i = 1; i <= 128; i++) {
      final p =
          curve.evaluate(curve.t0 + (curve.t1 - curve.t0) * i / 128) - origin;
      area += previous.cross(p);
      previous = p;
    }
    return area.length / (2 * math.pi);
  }

  static List<(Vector3, Vector3, BrepCurve)> _loopSegments(
    BrepSolid s,
    BrepLoop loop,
  ) {
    final out = <(Vector3, Vector3, BrepCurve)>[];
    for (final c in loop.coedges) {
      final e = s.edges[c.edge];
      var a = c.reversed ? s.vertices[e.v1] : s.vertices[e.v0];
      var b = c.reversed ? s.vertices[e.v0] : s.vertices[e.v1];
      final curve = c.reversed ? _flipped(e.curve) : e.curve;
      if (curve is BrepCircle && (a - b).length < 1e-8) {
        final on = curve.evaluate(curve.t0);
        if ((a - curve.center).length < curve.radius * 0.5) {
          a = on;
          b = on;
        }
      }
      out.add((a.clone(), b.clone(), curve));
    }
    return out;
  }

  static List<(Vector3, Vector3, BrepCurve)> _alignLoop(
    List<(Vector3, Vector3, BrepCurve)> outer,
    List<(Vector3, Vector3, BrepCurve)> inner,
  ) {
    var best = inner;
    var bestD = _loopPairCost(outer, inner);
    for (var r = 1; r < inner.length; r++) {
      final rot = [
        for (var i = 0; i < inner.length; i++) inner[(i + r) % inner.length],
      ];
      final d = _loopPairCost(outer, rot);
      if (d < bestD) {
        bestD = d;
        best = rot;
      }
    }
    final reversed = [
      for (var i = inner.length - 1; i >= 0; i--)
        (inner[i].$2, inner[i].$1, _flipped(inner[i].$3)),
    ];
    for (var r = 0; r < reversed.length; r++) {
      final rot = [
        for (var i = 0; i < reversed.length; i++)
          reversed[(i + r) % reversed.length],
      ];
      final d = _loopPairCost(outer, rot);
      if (d < bestD) {
        bestD = d;
        best = rot;
      }
    }
    return best;
  }

  static double _loopPairCost(
    List<(Vector3, Vector3, BrepCurve)> a,
    List<(Vector3, Vector3, BrepCurve)> b,
  ) {
    var d = 0.0;
    for (var i = 0; i < a.length; i++) {
      final am = (a[i].$1 + a[i].$2) * 0.5;
      final bm = (b[i].$1 + b[i].$2) * 0.5;
      d += (am - bm).length;
      d += (a[i].$1 - b[i].$1).length + (a[i].$2 - b[i].$2).length;
    }
    return d;
  }

  static BrepCurve _flipped(BrepCurve c) {
    if (c is BrepSplineCurve) return c.trimmed(c.t1, c.t0);
    if (c is BrepAffineCurve) {
      return BrepAffineCurve(_flipped(c.base), c.transform);
    }
    if (c is BrepPlaneSectionCurve) return c.trimmed(c.t1, c.t0);
    if (c is BrepSurfaceIntersectionCurve) return c.trimmed(c.t1, c.t0);
    if (c is BrepLine) return BrepLine(c.b, c.a);
    if (c is BrepCircle) {
      return BrepCircle(
        center: c.center,
        normal: c.normal,
        radius: c.radius,
        xDir: c.xDir,
        t0: c.t1,
        t1: c.t0,
      );
    }
    if (c is BrepEllipse) {
      return BrepEllipse(
        center: c.center,
        normal: c.normal,
        xDir: c.xDir,
        rMajor: c.rMajor,
        rMinor: c.rMinor,
        t0: c.t1,
        t1: c.t0,
      );
    }
    return c;
  }

  static BrepSolid _offsetInward(BrepSolid s, double t, Set<int> opening) {
    if (s.faces.length == 1 && s.faces.first.surface is BrepSphere) {
      final sp = s.faces.first.surface as BrepSphere;
      final r = sp.radius - t;
      if (r <= 1e-8) {
        throw BrepUnsupported('Hollow thickness removes the entire interior');
      }
      return BrepPrimitives.sphere(center: sp.center, radius: r);
    }

    // An opening tangent to a round must move with the offset round while
    // constructing the inner rim. Keeping its plane fixed makes the shrunk
    // round and that plane disjoint, so no vertex solver can satisfy both.
    bool offsetOpening(int fi) {
      if (!opening.contains(fi)) return true;
      final face = s.faces[fi];
      for (final loop in face.loops) {
        for (final co in loop.coedges) {
          final curve = s.edges[co.edge].curve;
          final p = curve.evaluate((curve.t0 + curve.t1) / 2);
          final uv = face.surface.parameters(p)!;
          final n = face.outwardNormal(uv.$1, uv.$2);
          for (final adj in BrepFillet.facesForEdge(s, co.edge)) {
            if (opening.contains(adj) || s.faces[adj].surface is BrepPlane) {
              continue;
            }
            final other = s.faces[adj];
            final ouv = other.surface.parameters(p)!;
            if (n.dot(other.outwardNormal(ouv.$1, ouv.$2)) > 1 - 1e-8) {
              return true;
            }
          }
        }
      }
      return false;
    }

    final offsetSurfs = <BrepSurface>[
      for (var i = 0; i < s.faces.length; i++)
        _offsetSurface(s.faces[i], offsetOpening(i) ? t : 0.0),
    ];
    final vFaces = <int, Set<int>>{};
    for (var fi = 0; fi < s.faces.length; fi++) {
      for (final loop in s.faces[fi].loops) {
        for (final c in loop.coedges) {
          final e = s.edges[c.edge];
          (vFaces[e.v0] ??= {}).add(fi);
          (vFaces[e.v1] ??= {}).add(fi);
        }
      }
    }

    final inner = s.clone();
    for (var vi = 0; vi < inner.vertices.length; vi++) {
      final fis = (vFaces[vi] ?? const <int>{}).toList();
      if (fis.isEmpty) continue;
      inner.vertices[vi] = _projectToSurfaces(
        s.vertices[vi],
        fis,
        offsetSurfs,
        s,
      );
    }
    for (var fi = 0; fi < inner.faces.length; fi++) {
      inner.faces[fi].surface = offsetSurfs[fi];
    }
    for (var ei = 0; ei < inner.edges.length; ei++) {
      inner.edges[ei].curve = _offsetCurve(s, inner, ei, offsetSurfs);
    }
    _snapClosedCircles(inner);
    for (var i = 0; i < s.edges.length; i++) {
      final old = s.edges[i], edge = inner.edges[i];
      if (old.curve is! BrepLine || old.v0 == old.v1) continue;
      final a = s.vertices[old.v1] - s.vertices[old.v0];
      final b = inner.vertices[edge.v1] - inner.vertices[edge.v0];
      if (a.dot(b) <= 1e-10 * a.length) {
        throw BrepUnsupported('Hollow thickness collapses an inner edge');
      }
    }
    return BrepHeal.alignCircularTrims(inner);
  }

  /// Cylinder caps store a full circle on an axis vertex. After offset, put
  /// that vertex on the circle so the seam is a wall generator.
  static void _snapClosedCircles(BrepSolid s) {
    for (final e in s.edges) {
      final c = e.curve;
      if (c is! BrepCircle || e.v0 != e.v1) continue;
      final on = c.evaluate(c.t0);
      if ((s.vertices[e.v0] - c.center).length < c.radius * 0.5) {
        s.vertices[e.v0] = on;
      }
    }
  }

  static BrepSurface _offsetSurface(BrepFace face, double t) {
    final s = face.surface;
    if (t == 0) return s;
    if (s is BrepPlane) {
      final n = face.reversed ? -s.normal : s.normal;
      return BrepPlane(
        origin: s.origin - n * t,
        normal: s.normal,
        xDir: s.xDir,
      );
    }
    if (s is BrepCylinder) {
      final sn = s.normalAt(0, 0);
      final out = face.outwardNormal(0, 0);
      final convex = out.dot(sn) > 0;
      if (convex && s.radius <= t + 1e-3) {
        throw BrepUnsupported(
          'Hollow thickness ${t.toStringAsFixed(2)} mm removes a '
          '${s.radius.toStringAsFixed(2)} mm fillet. Use a thinner wall '
          'to keep rounded inner corners.',
          code: BrepUnsupported.filletOffset,
        );
      }
      final r = convex ? s.radius - t : s.radius + t;
      return BrepCylinder(
        origin: s.origin,
        axis: s.axis,
        radius: math.max(r, 1e-4),
        xDir: s.xDir,
      );
    }
    if (s is BrepSphere) {
      final convex = !face.reversed;
      if (convex && s.radius <= t + 1e-3) {
        throw BrepUnsupported(
          'Hollow thickness ${t.toStringAsFixed(2)} mm removes a '
          '${s.radius.toStringAsFixed(2)} mm fillet. Use a thinner wall '
          'to keep rounded inner corners.',
          code: BrepUnsupported.filletOffset,
        );
      }
      final r = (convex ? s.radius - t : s.radius + t).abs();
      return BrepSphere(center: s.center, radius: math.max(r, 1e-4));
    }
    if (s is BrepTorus) {
      final r = face.reversed ? s.minorRadius + t : s.minorRadius - t;
      if (r <= 1e-8 || r >= s.majorRadius - 1e-8) {
        throw BrepUnsupported(
          'Hollow thickness collapses the toroidal wall',
          code: BrepUnsupported.filletOffset,
        );
      }
      return BrepTorus(
        center: s.center,
        axis: s.axis,
        majorRadius: s.majorRadius,
        minorRadius: r,
        xDir: s.xDir,
      );
    }
    if (s is BrepCone) {
      final sl = math.atan2(s.radius - s.topRadius, s.height);
      final out = face.outwardNormal(0, s.height * 0.5);
      final radial = s.xDir * math.cos(0) + s.yDir * math.sin(0);
      final convex = out.dot(radial) > 0;
      // At fixed axial height, the radial offset is t / cos(slope).
      final dr = t / math.cos(sl);
      final r0 = convex ? s.radius - dr : s.radius + dr;
      final r1 = convex ? s.topRadius - dr : s.topRadius + dr;
      // A pointed cone keeps an apex after offset, but at a different axial
      // position. Shorten its parameter span instead of treating the negative
      // radius at the OLD tip as a collapsed cavity. Frusta still require both
      // end radii: losing a cap would need a separate topology change.
      if (convex && s.topRadius == 0 && r0 > 1e-8) {
        return BrepCone(
          origin: s.origin,
          axis: s.axis,
          radius: r0,
          height: s.height * r0 / s.radius,
          xDir: s.xDir,
        );
      }
      if (convex && s.radius == 0 && r1 > 1e-8) {
        final shift = s.height * dr / s.topRadius;
        return BrepCone(
          origin: s.origin + s.axis * shift,
          axis: s.axis,
          radius: 0,
          topRadius: r1,
          height: s.height - shift,
          xDir: s.xDir,
        );
      }
      if (r0 <= 1e-8 || r1 <= 1e-8) {
        throw BrepUnsupported('Hollow thickness crosses the cone apex');
      }
      return BrepCone(
        origin: s.origin,
        axis: s.axis,
        radius: r0,
        topRadius: r1,
        height: s.height,
        xDir: s.xDir,
      );
    }
    final swept = _offsetExtrusion(face, t);
    if (swept != null) return swept;
    throw BrepUnsupported(
      'Hollow offset of ${s.runtimeType} is not implemented yet',
    );
  }

  static Vector3 _projectToSurfaces(
    Vector3 p,
    List<int> fis,
    List<BrepSurface> surfs,
    BrepSolid orig,
  ) {
    final planar = <int>[
      for (final fi in fis)
        if (surfs[fi] is BrepPlane) fi,
    ];
    var use = fis;
    var onClosedCircle = false;
    for (final e in orig.edges) {
      if (e.curve is BrepCircle && (e.v0 == e.v1)) {
        // Degenerate cylinder-cap encoding: one vertex at the axis.
        if ((orig.vertices[e.v0] - p).length < 1e-6) {
          onClosedCircle = true;
          break;
        }
      }
    }
    if (!onClosedCircle) {
      for (final fi in fis) {
        final s = orig.faces[fi].surface;
        if (s is BrepCylinder) {
          final d = p - s.origin;
          final rad = (d - s.axis * d.dot(s.axis)).length;
          if (rad < s.radius * 0.5 && planar.isNotEmpty) {
            use = planar;
          }
          break;
        }
      }
    }
    var q = p.clone();
    double residual(Vector3 at) {
      var worst = 0.0;
      for (final fi in use) {
        worst = math.max(worst, (surfs[fi].closestPoint(at) - at).length);
      }
      return worst;
    }

    // At an apex all generators meet. A single sampled cone normal cannot
    // constrain that vertex: a Newton step can slide it onto one generator.
    // Preserve the singular vertex explicitly on the offset cone's axis.
    for (final fi in use) {
      final before = orig.faces[fi].surface;
      final after = surfs[fi];
      if (before is! BrepCone || after is! BrepCone) continue;
      final slope = (before.topRadius - before.radius) / before.height;
      final newSlope = (after.topRadius - after.radius) / after.height;
      if (slope.abs() < 1e-15 || newSlope.abs() < 1e-15) continue;
      final oldApex = before.origin - before.axis * (before.radius / slope);
      if ((p - oldApex).length > 1e-7) continue;
      final apex = after.origin - after.axis * (after.radius / newSlope);
      if (residual(apex) > 1e-6) {
        throw BrepUnsupported(
          'Hollow offset surfaces do not meet at the cone apex',
        );
      }
      return apex;
    }

    for (var iter = 0; iter < 48; iter++) {
      final error = residual(q);
      if (error < 1e-9) return q;
      final matrix = Matrix3.zero();
      final rhs = Vector3.zero();
      final span = <Vector3>[];
      for (final fi in use) {
        final surface = surfs[fi];
        final closest = surface.closestPoint(q);
        final uv = surface.parameters(closest)!;
        final n = surface.normalAt(uv.$1, uv.$2);
        final distance = (q - closest).dot(n);
        final independent = n.clone();
        for (final axis in span) {
          independent.sub(axis * independent.dot(axis));
        }
        if (independent.length > 1e-7) span.add(independent.normalized());
        for (var i = 0; i < 3; i++) {
          rhs[i] -= n[i] * distance;
          for (var j = 0; j < 3; j++) {
            matrix.setEntry(i, j, matrix.entry(i, j) + n[i] * n[j]);
          }
        }
      }
      // Pin only unconstrained tangent directions. A tiny diagonal in an
      // otherwise singular normal matrix amplifies roundoff into seam motion.
      for (var i = 0; i < 3; i++) {
        for (var j = 0; j < 3; j++) {
          var nullProjection = i == j ? 1.0 : 0.0;
          for (final axis in span) {
            nullProjection -= axis[i] * axis[j];
          }
          matrix.setEntry(i, j, matrix.entry(i, j) + nullProjection);
        }
      }
      if (matrix.invert() == 0) break;
      final step = matrix.transform(rhs);
      var fraction = 1.0;
      while (fraction > 1.0 / 1024 && residual(q + step * fraction) >= error) {
        fraction *= .5;
      }
      q += step * fraction;
    }
    if (residual(q) > 1e-6) {
      throw BrepUnsupported(
        'Hollow offset surfaces do not meet at a stable vertex',
      );
    }
    return q;
  }

  static BrepCurve _offsetCurve(
    BrepSolid orig,
    BrepSolid inner,
    int edgeId,
    List<BrepSurface> offsetSurfs,
  ) {
    final e = inner.edges[edgeId];
    final a = inner.vertices[e.v0];
    final b = inner.vertices[e.v1];
    final oc = orig.edges[edgeId].curve;
    final neighbors = BrepFillet.facesForEdge(orig, edgeId);
    if (oc is! BrepLine) {
      final planes = [
        for (final fi in neighbors)
          if (offsetSurfs[fi] is BrepPlane) offsetSurfs[fi] as BrepPlane,
      ];
      final sweeps = [
        for (final fi in neighbors)
          if (offsetSurfs[fi] is BrepRuledSurface)
            offsetSurfs[fi] as BrepRuledSurface,
      ];
      if (planes.length == 1 && sweeps.length == 1) {
        final wall = sweeps.single;
        var ua = wall.parameters(a).$1, ub = wall.parameters(b).$1;
        if (e.v0 == e.v1 &&
            (oc.evaluate(oc.t0) - oc.evaluate(oc.t1)).length < 1e-8) {
          ua = 0;
          ub = 1;
          if (wall.du(0, 0).dot(oc.tangent(oc.t0) * (oc.t1 - oc.t0).sign) < 0) {
            ua = 1;
            ub = 0;
          }
        }
        if (e.v0 != e.v1 &&
            (ub - ua) *
                    wall
                        .du(ua, 0)
                        .dot(oc.tangent(oc.t0) * (oc.t1 - oc.t0).sign) <=
                1e-12) {
          throw BrepUnsupported(
            'Hollow thickness collapses a curved inner edge. Use a thinner wall.',
          );
        }
        final curve = BrepPlaneSectionCurve(
          surface: wall,
          plane: planes.single,
          t0: ua,
          t1: ub,
        );
        if ((curve.evaluate(ua) - a).length > 1e-6 ||
            (curve.evaluate(ub) - b).length > 1e-6) {
          throw BrepUnsupported(
            'Hollow curved boundary misses its offset vertex',
          );
        }
        final lower = wall.lower;
        final delta = wall.evaluate(0, 1) - wall.evaluate(0, 0);
        final plane = planes.single,
            denominator = delta.dot(planes.single.normal);
        if (lower is BrepSplineCurve &&
            denominator.abs() > 1e-12 * delta.length) {
          // A plane section of a translational spline sweep is an affine
          // projection of its rail: preserve the exact rational control net.
          return BrepSplineCurve(
            degree: lower.degree,
            points: [
              for (final p in lower.points)
                p - delta * (plane.height(p) / denominator),
            ],
            knots: lower.knots,
            weights: lower.weights,
            start: lower.t0 + (lower.t1 - lower.t0) * ua,
            end: lower.t0 + (lower.t1 - lower.t0) * ub,
          );
        }
        return curve;
      }
      if (planes.length == 2) {
        final direction = b - a;
        final straight =
            direction.length > 1e-10 &&
            [0.0, .25, .5, .75, 1.0].every((u) {
              final p =
                  oc.evaluate(oc.t0 + (oc.t1 - oc.t0) * u) - oc.evaluate(oc.t0);
              return p.cross(direction).length < 1e-7 * direction.length;
            });
        if (straight) return BrepLine(a, b);
      }
    }
    if (oc is BrepLine) {
      return BrepLine(a, b);
    }
    if (oc is BrepCircle) {
      final adj = BrepFillet.facesForEdge(orig, edgeId);
      BrepPlane? plane;
      BrepCylinder? cyl;
      BrepTorus? torus;
      for (final fi in adj) {
        final s = offsetSurfs[fi];
        if (s is BrepPlane) plane = s;
        if (s is BrepCylinder) cyl = s;
        if (s is BrepTorus) torus = s;
      }
      if (plane != null && cyl != null) {
        final denom = plane.normal.dot(cyl.axis);
        Vector3 center;
        if (denom.abs() > 1e-8) {
          final along = -plane.height(cyl.origin) / denom;
          center = cyl.origin + cyl.axis * along;
        } else {
          center = plane.closestPoint(oc.center);
        }
        return BrepCircle(
          center: center,
          normal: oc.normal,
          radius: cyl.radius,
          xDir: oc.xDir,
          t0: oc.t0,
          t1: oc.t1,
        );
      }
      if (torus != null) {
        // A torus has both latitude circles and tube-section circles. Their
        // centres are not generally the torus centre (e.g. a tangent cap).
        final delta = a - torus.center;
        final axial = torus.axis * delta.dot(torus.axis);
        final radial = delta - axial;
        final parallel = oc.normal.dot(torus.axis).abs();
        final Vector3 center;
        if (parallel > 1 - 1e-8) {
          center = torus.center + axial;
        } else if (parallel < 1e-8 && radial.length > 1e-8) {
          center = torus.center + radial.normalized() * torus.majorRadius;
        } else {
          throw BrepUnsupported(
            'Hollow needs a noncircular torus intersection',
          );
        }
        return BrepCircle(
          center: center,
          normal: oc.normal,
          radius: (a - center).length,
          xDir: oc.xDir,
          t0: oc.t0,
          t1: oc.t1,
        );
      }
      if (plane != null) {
        final center = plane.closestPoint(oc.center);
        return BrepCircle(
          center: center,
          normal: oc.normal,
          radius: (a - center).length,
          xDir: oc.xDir,
          t0: oc.t0,
          t1: oc.t1,
        );
      }
      return BrepCircle(
        center: oc.center,
        normal: oc.normal,
        radius: (a - oc.center).length,
        xDir: oc.xDir,
        t0: oc.t0,
        t1: oc.t1,
      );
    }
    if (oc is BrepEllipse) {
      final adj = BrepFillet.facesForEdge(orig, edgeId);
      if (adj.length == 2) {
        final intersections = BrepIntersect.surfaces(
          offsetSurfs[adj[0]],
          offsetSurfs[adj[1]],
        );
        for (final candidate in intersections.whereType<BrepEllipse>()) {
          if ((candidate.closestPoint(a) - a).length > 1e-6 ||
              (candidate.closestPoint(b) - b).length > 1e-6) {
            continue;
          }
          final oriented = BrepEllipse(
            center: candidate.center,
            normal: oc.normal,
            xDir: candidate.xDir,
            rMajor: candidate.rMajor,
            rMinor: candidate.rMinor,
          );
          final start = oriented.closestParam(a);
          var end = oriented.closestParam(b);
          if (oc.isClosed) {
            end = start + 2 * math.pi;
          } else {
            while (end <= start) {
              end += 2 * math.pi;
            }
          }
          return BrepEllipse(
            center: oriented.center,
            normal: oriented.normal,
            xDir: oriented.xDir,
            rMajor: oriented.rMajor,
            rMinor: oriented.rMinor,
            t0: start,
            t1: end,
          );
        }
      }
      throw BrepUnsupported(
        'Hollow could not construct the offset elliptical boundary',
      );
    }
    throw BrepUnsupported(
      'Hollow needs an offset intersection for ${oc.runtimeType}',
    );
  }

  /// Map picked mesh triangles to B-Rep faces for an opening.
  ///
  /// [MeshData.faceIds] are a hint only — host box meshes number faces
  /// differently from [BrepPrimitives.box]. Each candidate must face
  /// [openingNormal] and sit on the picked opening, or we fall back to
  /// a geometric match.
  static List<int> facesFromMeshOpening(
    BrepSolid solid,
    MeshData mesh,
    List<int> meshFaces,
    Vector3 openingNormal,
  ) {
    final n = openingNormal.length2 < 1e-20
        ? Vector3(0, 0, 1)
        : openingNormal.normalized();
    var centroid = Vector3.zero();
    var counted = 0;
    for (final f in meshFaces) {
      if (f < 0 || f >= mesh.faces.length) continue;
      for (final i in mesh.faces[f]) {
        centroid += mesh.vertices[i];
        counted++;
      }
    }
    if (counted > 0) centroid /= counted.toDouble();

    List<int> pick(Iterable<int> candidates) {
      final scored = <(int, double)>[];
      for (final fi in candidates) {
        if (fi < 0 || fi >= solid.faces.length) continue;
        final face = solid.faces[fi];
        final fn = face.outwardNormal(0, 0);
        if (fn.length2 < 1e-20 || fn.dot(n) < 0.85) continue;
        final sample = _faceSample(solid, fi);
        final rel = sample - centroid;
        scored.add((fi, rel.dot(n).abs() + (rel - n * rel.dot(n)).length));
      }
      if (scored.isEmpty) return const [];
      scored.sort((a, b) => a.$2.compareTo(b.$2));
      final best = scored.first.$2;
      return [
        for (final s in scored)
          if (s.$2 <= best + 1e-4) s.$1,
      ];
    }

    final ids = mesh.faceIds;
    if (ids != null) {
      final hinted = <int>{};
      for (final f in meshFaces) {
        if (f < 0 || f >= ids.length) continue;
        final id = ids[f];
        if (id >= 0 && id < solid.faces.length) hinted.add(id);
      }
      final keep = pick(hinted);
      if (keep.isNotEmpty) return keep.toSet().toList()..sort();
    }
    final keep = pick([for (var i = 0; i < solid.faces.length; i++) i]);
    if (keep.isEmpty) {
      throw BrepUnsupported('Hollow could not match the opening face');
    }
    return keep.toSet().toList()..sort();
  }

  static Vector3 _faceSample(BrepSolid s, int fi) {
    final face = s.faces[fi];
    var sum = Vector3.zero();
    var n = 0;
    for (final loop in face.loops) {
      for (final c in loop.coedges) {
        final e = s.edges[c.edge];
        sum += s.vertices[e.v0];
        n++;
      }
    }
    if (n == 0) return face.surface.evaluate(0, 0);
    return sum / n.toDouble();
  }
}
