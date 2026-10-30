import 'dart:math' as math;

import 'package:vector_math/vector_math_64.dart';

import '../mesh/mesh_data.dart';
import '../mesh/mesh_polygon.dart';
import '../tessellation.dart';
import '../ijo_log.dart';
import 'geom.dart';
import 'solid.dart';
import 'xform.dart';

class TessellatedBrep {
  TessellatedBrep(this.mesh);
  final MeshData mesh;
}

/// Circular display chords must stay inside a narrow annular face. Apply the
/// density to the entire body so its adjoining faces share the same samples.
class _BoundaryTessellation extends MeshTessellation {
  _BoundaryTessellation(MeshTessellation base, this.minimum)
    : super(
        base.quality,
        surfaceChordToleranceMm: base.surfaceChordToleranceMm,
        maxSurfaceRefinement: base.maxSurfaceRefinement,
      );
  final int minimum;
  @override
  int get cylinderSegments => math.max(super.cylinderSegments, minimum);
  @override
  int filletArcSegments(double radius, double angleRad) => math.max(
    super.filletArcSegments(radius, angleRad),
    (minimum * angleRad.abs() / (2 * math.pi)).ceil(),
  );
}

/// Tessellates a B-Rep into display triangles tagged with [MeshData.faceIds].
class BrepTessellate {
  BrepTessellate._();

  /// [onProgress] counts completed faces, not elapsed processing time.
  static MeshData run(
    BrepSolid solid, {
    MeshTessellation tessellation = MeshTessellation.standard,
    void Function(int completed, int total)? onProgress,
  }) {
    final transformed = BrepXform.unwrapAffine(solid);
    if (transformed != null) {
      final mesh = run(
        transformed.solid,
        tessellation: tessellation,
        onProgress: onProgress,
      );
      final m = transformed.transform;
      final vertices = [for (final p in mesh.vertices) m.point(p)];
      final edgeIds = <String, int>{};
      for (final face in mesh.faces) {
        for (var i = 0; i < face.length; i++) {
          final a = face[i], b = face[(i + 1) % face.length];
          final id =
              mesh.edgeIds?[MeshData.segmentKey(
                mesh.vertices[a],
                mesh.vertices[b],
              )];
          if (id != null)
            edgeIds[MeshData.segmentKey(vertices[a], vertices[b])] = id;
        }
      }
      return MeshData(
        vertices: vertices,
        faces: m.determinant < 0
            ? [for (final f in mesh.faces) f.reversed.toList()]
            : mesh.faces,
        faceIds: mesh.faceIds,
        edgeIds: edgeIds,
      );
    }
    var minimumSegments = tessellation.cylinderSegments;
    for (final face in solid.faces) {
      if (face.surface is! BrepPlane || face.loops.length != 2) continue;
      // A polygonal recess can approach a round rim more closely than its
      // display chord error (e.g. the corners of a cross-head screw socket).
      // Bound each straight side by its endpoints and preserve that clearance.
      final outerUses = face.loops.first.coedges;
      if (outerUses.length == 1 &&
          solid.edges[outerUses.single.edge].curve is BrepCircle &&
          face.loops
              .skip(1)
              .every(
                (l) => l.coedges.every(
                  (co) => solid.edges[co.edge].curve is BrepLine,
                ),
              )) {
        final rim = solid.edges[outerUses.single.edge].curve as BrepCircle;
        if (rim.isClosed) {
          var extent = 0.0;
          for (final l in face.loops.skip(1))
            for (final co in l.coedges) {
              final line = solid.edges[co.edge].curve as BrepLine;
              extent = math.max(
                extent,
                math.max(
                  (line.a - rim.center).length,
                  (line.b - rim.center).length,
                ),
              );
            }
          final clearance = rim.radius - extent;
          if (clearance > 1e-9 && extent > 0) {
            final angle = math.acos(
              (1 - clearance / (4 * rim.radius)).clamp(-1.0, 1.0),
            );
            if (angle > 0)
              minimumSegments = math.max(
                minimumSegments,
                (math.pi / angle).ceil(),
              );
          }
        }
      }
      final circles = <BrepCircle>[];
      var circular = true;
      for (final loop in face.loops) {
        BrepCircle? first;
        var sweep = 0.0;
        for (final co in loop.coedges) {
          final c = solid.edges[co.edge].curve;
          if (c is! BrepCircle) {
            circular = false;
            break;
          }
          first ??= c;
          if ((c.center - first.center).length > 1e-7 ||
              (c.radius - first.radius).abs() > 1e-7 ||
              c.normal.cross(first.normal).length > 1e-7)
            circular = false;
          sweep += (c.t1 - c.t0).abs();
        }
        if (first == null || (sweep - 2 * math.pi).abs() > 1e-6) {
          circular = false;
        } else {
          circles.add(first);
        }
      }
      if (!circular ||
          circles.length != 2 ||
          (circles[0].center - circles[1].center).length > 1e-7)
        continue;
      final outer = math.max(circles[0].radius, circles[1].radius);
      final gap = (circles[0].radius - circles[1].radius).abs();
      if (gap <= 1e-9 || outer <= 0) continue;
      final angle = math.acos((1 - gap / (4 * outer)).clamp(-1.0, 1.0));
      if (angle > 0)
        minimumSegments = math.max(minimumSegments, (math.pi / angle).ceil());
    }
    if (minimumSegments > tessellation.cylinderSegments) {
      tessellation = _BoundaryTessellation(tessellation, minimumSegments);
    }
    final trimmedRuled = solid.edges.any(
      (e) =>
          e.curve is BrepSurfaceIntersectionCurve ||
          e.curve is BrepRuledQuadricCurve ||
          e.curve is BrepPlaneSectionCurve,
    );
    final verts = <Vector3>[];
    final faces = <List<int>>[];
    final faceIds = <int>[];

    void addTri(Vector3 a, Vector3 b, Vector3 c, int faceId) {
      final n = (b - a).cross(c - a);
      if (n.length2 < 1e-28) return;
      final i = verts.length;
      verts.add(a.clone());
      verts.add(b.clone());
      verts.add(c.clone());
      faces.add([i, i + 1, i + 2]);
      faceIds.add(faceId);
    }

    for (var fi = 0; fi < solid.faces.length; fi++) {
      onProgress?.call(fi, solid.faces.length);
      final face = solid.faces[fi];
      final s = face.surface;
      if (IjoLog.verbose)
        IjoLog.line(
          'tessellate face $fi/${solid.faces.length} ${s.runtimeType} triangles=${faces.length}',
        );
      if (s is BrepAffineSurface) {
        final b = BrepBuilder(), vm = <int, int>{}, em = <int, int>{};
        final inverse = BrepAffineTransform.fromFunction(
          s.transform.inversePoint,
        );
        int vertex(int i) => vm.putIfAbsent(
          i,
          () => b.vertex(s.transform.inversePoint(solid.vertices[i])),
        );
        int edge(int i) => em.putIfAbsent(i, () {
          final e = solid.edges[i], c = e.curve;
          return b.edge(
            vertex(e.v0),
            vertex(e.v1),
            c is BrepAffineCurve && c.transform.sameAs(s.transform)
                ? c.base
                : BrepAffineCurve(c, inverse),
          );
        });
        b.face(s.base, [
          for (final l in face.loops)
            BrepLoop([
              for (final co
                  in (s.transform.determinant < 0
                      ? l.coedges.reversed
                      : l.coedges))
                BrepCoedge(
                  edge(co.edge),
                  reversed: co.reversed != (s.transform.determinant < 0),
                ),
            ]),
        ], reversed: face.reversed);
        final mesh = run(b.build(), tessellation: tessellation),
            points = mesh.vertices.map(s.transform.point).toList();
        for (final t in mesh.faces) {
          if (s.transform.determinant < 0) {
            addTri(points[t[0]], points[t[2]], points[t[1]], fi);
          } else {
            addTri(points[t[0]], points[t[1]], points[t[2]], fi);
          }
        }
      } else if (s is BrepSplineSurface || s is BrepOffsetSurface) {
        _parametric(solid, face, fi, tessellation, addTri);
      } else if (s is BrepRuledSurface) {
        // Exact sweep rectangles share the same samples as their cap edges,
        // including affine spline rails and circular rails in mixed profiles.
        if (_tensorRectangle(solid, face, fi, tessellation, addTri)) {
          continue;
        } else if (trimmedRuled ||
            s.lower is BrepSplineCurve ||
            s.upper is BrepSplineCurve ||
            s.lower is BrepAffineCurve ||
            s.upper is BrepAffineCurve ||
            face.loops.any(
              (loop) => loop.coedges.any(
                (co) =>
                    solid.edges[co.edge].curve is BrepSplineCurve ||
                    solid.edges[co.edge].curve is BrepPlaneSectionCurve ||
                    solid.edges[co.edge].curve is BrepRuledQuadricCurve,
              ),
            )) {
          _parametric(solid, face, fi, tessellation, addTri);
        } else {
          _ruled(solid, s, face, fi, tessellation, addTri);
        }
      } else if (s is BrepPlane) {
        _plane(solid, face, fi, tessellation, addTri);
      } else if (s is BrepCylinder) {
        _cylinder(solid, face, fi, tessellation, addTri);
      } else if (s is BrepSphere) {
        _sphere(solid, s, face, fi, tessellation, addTri);
      } else if (s is BrepTorus) {
        _torus(solid, s, face, fi, tessellation, addTri);
      } else if (s is BrepCone) {
        _cone(solid, s, face, fi, tessellation, addTri);
      } else {
        throw BrepUnsupported(
          'Surface ${s.runtimeType} has no display tessellator',
        );
      }
    }

    final mesh = _tagBrepEdges(
      solid,
      _weld(
        MeshData(vertices: verts, faces: faces, faceIds: faceIds),
        eps2:
            solid.edges.any(
              (e) =>
                  e.curve is BrepSurfaceIntersectionCurve ||
                  e.curve is BrepRuledQuadricCurve,
            )
            ? 1e-14
            : 1e-10,
      ),
      tessellation.cylinderSegments,
    );
    onProgress?.call(solid.faces.length, solid.faces.length);
    return mesh;
  }

  /// A ruled strip trimmed at two constant sweep parameters is still ruled.
  /// Zipper its actual boundary samples; independent UV fan refinement creates
  /// unnecessary triangles and can lose the fitted contact boundary.
  static bool _ruledConstantBand(
    BrepSolid solid,
    BrepFace face,
    int fi,
    MeshTessellation tess,
    void Function(Vector3, Vector3, Vector3, int) emit,
  ) {
    final s = face.surface;
    // Helical strips have a curvature-controlled mesher of their own.
    if (s is BrepRuledSurface && (s.lower is BrepHelix || s.upper is BrepHelix))
      return false;
    if (s is! BrepRuledSurface || face.loops.length != 1) return false;
    final sweep = s.dv(0, .5);
    if ([.25, .5, .75, 1.0].any((u) => (s.dv(u, .5) - sweep).length > 1e-7))
      return false;
    if ((s.evaluate(0, .5) - s.evaluate(1, .5)).length < 1e-8) return false;
    final bands = <double, List<(double, Vector3)>>{}, sides = <double>[];
    for (final co in face.loops.single.coedges) {
      final curve = solid.edges[co.edge].curve;
      final probes = [
        for (final t in [0.0, .25, .5, .75, 1.0])
          s.parameters(curve.evaluate(curve.t0 + (curve.t1 - curve.t0) * t)),
      ];
      if (probes.every((p) => (p.$2 - probes.first.$2).abs() < 1e-7)) {
        final v =
            bands.keys
                .where((v) => (v - probes.first.$2).abs() < 1e-7)
                .firstOrNull ??
            probes.first.$2;
        final row = bands[v] ??= [];
        final points = sampleLoop(
          solid,
          BrepLoop([co]),
          tess.cylinderSegments,
          tess: tess,
        );
        points.add(
          solid.vertices[co.reversed
              ? solid.edges[co.edge].v0
              : solid.edges[co.edge].v1],
        );
        for (final p in points) {
          final uv = s.parameters(p);
          if ((s.evaluate(uv.$1, v) - p).length > 1e-5) return false;
          row.add((uv.$1, p));
        }
      } else if (curve is BrepLine &&
          probes.every((p) => (p.$1 - probes.first.$1).abs() < 1e-7)) {
        sides.add(probes.first.$1);
      } else {
        return false;
      }
    }
    if (bands.length != 2 ||
        sides.length != 2 ||
        (sides[0] - sides[1]).abs() < 1e-8)
      return false;
    final levels = bands.keys.toList()..sort();
    List<(double, Vector3)> row(double v) {
      final sorted = bands[v]!..sort((a, b) => a.$1.compareTo(b.$1));
      final out = <(double, Vector3)>[];
      for (final p in sorted) {
        if (out.isEmpty || p.$1 - out.last.$1 > 1e-7) out.add(p);
      }
      return out;
    }

    final a = row(levels[0]), b = row(levels[1]);
    if (a.length < 2 ||
        b.length < 2 ||
        (a.first.$1 - b.first.$1).abs() > 1e-7 ||
        (a.last.$1 - b.last.$1).abs() > 1e-7)
      return false;
    void tri(Vector3 a, Vector3 b, Vector3 c) {
      if (face.reversed)
        emit(a, c, b, fi);
      else
        emit(a, b, c, fi);
    }

    var i = 0, j = 0;
    while (i < a.length - 1 || j < b.length - 1) {
      final ua = i < a.length - 1 ? a[i + 1].$1 : double.infinity,
          ub = j < b.length - 1 ? b[j + 1].$1 : double.infinity;
      if ((ua - ub).abs() < 1e-7) {
        tri(a[i].$2, a[i + 1].$2, b[j].$2);
        tri(a[i + 1].$2, b[j + 1].$2, b[j].$2);
        i++;
        j++;
      } else if (ua < ub) {
        tri(a[i].$2, a[i + 1].$2, b[j].$2);
        i++;
      } else {
        tri(a[i].$2, b[j + 1].$2, b[j].$2);
        j++;
      }
    }
    return true;
  }

  static bool _tensorRectangle(
    BrepSolid solid,
    BrepFace face,
    int fi,
    MeshTessellation tess,
    void Function(Vector3, Vector3, Vector3, int) emit,
  ) {
    if (face.loops.length != 1 || face.loops.single.coedges.length != 4)
      return false;
    final s = face.surface;
    // Helical strips have a curvature-controlled mesher of their own.
    if (s is BrepRuledSurface && (s.lower is BrepHelix || s.upper is BrepHelix))
      return false;
    // Only use the spline grid shortcut for circular-section patches.
    // General freeform interiors retain the adaptive parameter-domain mesher.
    if (s is BrepSplineSurface &&
        (s.degreeV != 2 ||
            s.points.first.length != 3 ||
            s.weights.any(
              (row) =>
                  (row[0] - row[2]).abs() > 1e-12 ||
                  (row[1] / row[0] - math.sqrt(.5)).abs() > 1e-12,
            )))
      return false;
    final u0 = s is BrepSplineSurface ? s.u0 : 0.0,
        u1 = s is BrepSplineSurface ? s.u1 : 1.0;
    final v0 = s is BrepSplineSurface ? s.v0 : 0.0,
        v1 = s is BrepSplineSurface ? s.v1 : 1.0;
    Vector3 side(int index, double t) => switch (index) {
      0 => s.evaluate(u0 + (u1 - u0) * t, v0),
      1 => s.evaluate(u1, v0 + (v1 - v0) * t),
      2 => s.evaluate(u0 + (u1 - u0) * t, v1),
      _ => s.evaluate(u0, v0 + (v1 - v0) * t),
    };
    final counts = <int, int>{};
    for (final co in face.loops.single.coedges) {
      final c = solid.edges[co.edge].curve;
      int? matched;
      for (var i = 0; i < 4; i++) {
        for (final reverse in [false, true]) {
          if ([0.0, .125, .25, .5, .75, .875, 1.0].every(
            (t) =>
                (c.evaluate(c.t0 + (c.t1 - c.t0) * t) -
                        side(i, reverse ? 1 - t : t))
                    .length <
                1e-8,
          ))
            matched = i;
        }
      }
      if (matched == null || counts.containsKey(matched)) return false;
      counts[matched] = sampleLoop(
        solid,
        BrepLoop([co]),
        tess.cylinderSegments,
        tess: tess,
      ).length;
    }
    final transpose = counts[0] != counts[2];
    if (transpose && counts[1] != counts[3]) return false;
    final nu = transpose ? counts[1]! : counts[0]!;
    final startAcross = transpose ? counts[0]! : counts[3]!;
    final endAcross = transpose ? counts[2]! : counts[1]!;
    final nv = math.max(startAcross, endAcross);
    if (s is BrepSplineSurface && (nu < 4 || nv < 4)) {
      // A straight longitudinal blend needs only its endpoint columns. Check
      // the control net, not just boundary samples, before reducing its grid.
      final linearU =
          !transpose &&
          nu == 1 &&
          s.degreeU == 3 &&
          s.points.length == 4 &&
          List.generate(3, (j) => j).every(
            (j) => List.generate(4, (i) => i).every(
              (i) =>
                  (s.points[i][j] -
                              (s.points[0][j] * (1 - i / 3) +
                                  s.points[3][j] * (i / 3)))
                          .length <
                      1e-8 &&
                  (s.weights[i][j] - s.weights[0][j]).abs() < 1e-12,
            ),
          );
      if (!linearU || nv < 4) return false;
    }
    if (nu == 1 &&
        nv == 1 &&
        (s.evaluate((u0 + u1) * .5, (v0 + v1) * .5) -
                    (s.evaluate(u0, v0) + s.evaluate(u1, v1)) * .5)
                .length >
            1e-8)
      return false;
    void triangle(Vector3 a, Vector3 b, Vector3 c) {
      if (face.reversed != transpose) {
        emit(a, c, b, fi);
      } else {
        emit(a, b, c, fi);
      }
    }

    List<Vector3> column(int i) {
      final count = i == 0
          ? startAcross
          : i == nu
          ? endAcross
          : nv;
      return [
        for (var j = 0; j <= count; j++)
          s.evaluate(
            u0 + (u1 - u0) * (transpose ? j / count : i / nu),
            v0 + (v1 - v0) * (transpose ? i / nu : j / count),
          ),
      ];
    }

    // A sharp miter has an elliptic section while a tangent join has a
    // circular one. Zip their different display sample counts without adding
    // unmatched vertices to either shared boundary.
    var left = column(0);
    for (var i = 0; i < nu; i++) {
      final right = column(i + 1);
      var j = 0, k = 0;
      while (j < left.length - 1 || k < right.length - 1) {
        final tj = j < left.length - 1
            ? (j + 1) / (left.length - 1)
            : double.infinity;
        final tk = k < right.length - 1
            ? (k + 1) / (right.length - 1)
            : double.infinity;
        if ((tj - tk).abs() < 1e-12) {
          triangle(left[j], right[k], right[k + 1]);
          triangle(left[j], right[k + 1], left[j + 1]);
          j++;
          k++;
        } else if (tj < tk) {
          triangle(left[j], right[k], left[j + 1]);
          j++;
        } else {
          triangle(left[j], right[k], right[k + 1]);
          k++;
        }
      }
      left = right;
    }
    return true;
  }

  /// Tessellates a bounded parameter domain while preserving shared edge samples.
  static void _parametric(
    BrepSolid solid,
    BrepFace face,
    int fi,
    MeshTessellation tess,
    void Function(Vector3, Vector3, Vector3, int) emit,
  ) {
    final s = face.surface;
    final locator = s is BrepSplineSurface ? BrepSplineSurfaceLocator(s) : null;
    (double, double) locate(Vector3 p) =>
        locator?.parameters(p) ??
        (s.parameters(p) ??
            (throw BrepUnsupported('Cannot resolve surface boundary')));
    if (_tensorRectangle(solid, face, fi, tess, emit)) return;
    if (_ruledConstantBand(solid, face, fi, tess, emit)) return;
    if (face.loops.isEmpty && s is BrepSplineSurface) {
      final n = math.max(24, tess.cylinderSegments);
      for (var i = 0; i < n; i++) {
        for (var j = 0; j < n; j++) {
          final a = s.evaluate(
                s.u0 + (s.u1 - s.u0) * i / n,
                s.v0 + (s.v1 - s.v0) * j / n,
              ),
              b = s.evaluate(
                s.u0 + (s.u1 - s.u0) * (i + 1) / n,
                s.v0 + (s.v1 - s.v0) * j / n,
              ),
              c = s.evaluate(
                s.u0 + (s.u1 - s.u0) * (i + 1) / n,
                s.v0 + (s.v1 - s.v0) * (j + 1) / n,
              ),
              d = s.evaluate(
                s.u0 + (s.u1 - s.u0) * i / n,
                s.v0 + (s.v1 - s.v0) * (j + 1) / n,
              );
          if (face.reversed) {
            emit(a, c, b, fi);
            emit(a, d, c, fi);
          } else {
            emit(a, b, c, fi);
            emit(a, c, d, fi);
          }
        }
      }
      return;
    }
    final chart = s is BrepSplineSurface
        ? s
        : s is BrepOffsetSurface && s.base is BrepSplineSurface
        ? s.base as BrepSplineSurface
        : null;
    // Offset surfaces expose a normalized chart even when their base does not.
    final pu0 = s is BrepSplineSurface ? s.u0 : 0.0;
    final pu1 = s is BrepSplineSurface ? s.u1 : 1.0;
    final pv0 = s is BrepSplineSurface ? s.v0 : 0.0;
    final pv1 = s is BrepSplineSurface ? s.v1 : 1.0;
    final offsetExtrusion =
        s is BrepOffsetSurface &&
        chart != null &&
        chart.degreeV == 1 &&
        chart.points.first.length == 2 &&
        List.generate(chart.points.length, (i) => i).every(
          (i) =>
              (chart.weights[i][0] - chart.weights[i][1]).abs() < 1e-12 &&
              ((chart.points[i][1] - chart.points[i][0]) -
                          (chart.points[0][1] - chart.points[0][0]))
                      .length <
                  1e-8,
        );
    var bandLoops = face.loops;
    if ((s is BrepRuledSurface ||
            s is BrepSplineSurface && s.ruled != null ||
            offsetExtrusion) &&
        face.loops.length == 1 &&
        [0.0, .5, 1.0].every(
          (v) =>
              (s.evaluate(pu0, pv0 + (pv1 - pv0) * v) -
                      s.evaluate(pu1, pv0 + (pv1 - pv0) * v))
                  .length <
              1e-7,
        )) {
      final levels = <double, List<BrepCoedge>>{};
      var valid = true;
      for (final co in face.loops.single.coedges) {
        final c = solid.edges[co.edge].curve;
        final uv = [
          for (final t in [0.0, .5, 1.0])
            locate(c.evaluate(c.t0 + (c.t1 - c.t0) * t)),
        ];
        if ((uv[0].$2 - uv[1].$2).abs() < 1e-7 &&
            (uv[0].$2 - uv[2].$2).abs() < 1e-7) {
          final level =
              levels.keys
                  .where((v) => (v - uv[0].$2).abs() < 1e-7)
                  .firstOrNull ??
              uv[0].$2;
          (levels[level] ??= []).add(co);
        } else if (face.loops.single.coedges
                .where((other) => other.edge == co.edge)
                .length !=
            2) {
          valid = false;
          break;
        }
      }
      if (valid && levels.length == 2)
        bandLoops = [for (final uses in levels.values) BrepLoop(uses)];
    }
    if ((s is BrepSplineSurface || s is BrepRuledSurface || offsetExtrusion) &&
        bandLoops.length == 2 &&
        (s is BrepRuledSurface ||
            s is BrepSplineSurface && s.degreeV == 1 ||
            offsetExtrusion) &&
        [0.0, .25, .5, .75, 1.0].every(
          (t) =>
              (s.evaluate(pu0, pv0 + (pv1 - pv0) * t) -
                      s.evaluate(pu1, pv0 + (pv1 - pv0) * t))
                  .length <
              1e-7,
        )) {
      final rings = <List<(double, Vector3)>>[], levels = <double>[];
      for (final loop in bandLoops) {
        final points = sampleLoop(
              solid,
              loop,
              tess.cylinderSegments,
              tess: tess,
            ),
            ring = <(double, Vector3)>[];
        var v = 0.0;
        for (final p in points) {
          final uv = locate(p);
          var u = uv.$1;
          if ((u - pu1).abs() < 1e-7) u = pu0;
          ring.add((u, p));
          v += uv.$2;
        }
        levels.add(v / points.length);
        ring.sort((a, b) => a.$1.compareTo(b.$1));
        // Close each periodic row with the same position one parameter turn on.
        ring.add((ring.first.$1 + pu1 - pu0, ring.first.$2));
        rings.add(ring);
      }
      if (levels[0] > levels[1]) {
        final r = rings[0];
        rings[0] = rings[1];
        rings[1] = r;
        final v = levels[0];
        levels[0] = levels[1];
        levels[1] = v;
      }
      final knots = {
        for (final ring in rings)
          for (final p in ring) p.$1,
      }.toList()..sort();
      final us = <double>[];
      for (final u in knots) {
        if (us.isEmpty || u - us.last > 1e-7) us.add(u);
      }
      final rows = [
        rings[0],
        for (var row = 1; row < 8; row++)
          [
            for (final u in us)
              (
                u,
                s.evaluate(
                  u > pu1 ? u - (pu1 - pu0) : u,
                  levels[0] + (levels[1] - levels[0]) * row / 8,
                ),
              ),
          ],
        rings[1],
      ];
      void triangle(Vector3 a, Vector3 b, Vector3 c) {
        if (face.reversed)
          emit(a, c, b, fi);
        else
          emit(a, b, c, fi);
      }

      for (var row = 0; row < 8; row++) {
        final a = rows[row], b = rows[row + 1];
        var i = 0, j = 0;
        while (i < a.length - 1 || j < b.length - 1) {
          final ua = i < a.length - 1 ? a[i + 1].$1 : double.infinity,
              ub = j < b.length - 1 ? b[j + 1].$1 : double.infinity;
          if ((ua - ub).abs() < 1e-7) {
            triangle(a[i].$2, a[i + 1].$2, b[j].$2);
            triangle(a[i + 1].$2, b[j + 1].$2, b[j].$2);
            i++;
            j++;
          } else if (ua < ub) {
            triangle(a[i].$2, a[i + 1].$2, b[j].$2);
            i++;
          } else {
            triangle(a[i].$2, b[j + 1].$2, b[j].$2);
            j++;
          }
        }
      }
      return;
    }
    final builder = BrepBuilder(), world = <String, Vector3>{};
    String key(Vector3 p) => '${(p.x * 1e10).round()}:${(p.y * 1e10).round()}';
    final canWrap =
        s is BrepSplineSurface ||
        s is BrepOffsetSurface ||
        s is BrepRuledSurface;
    final periodicU =
        canWrap &&
        [0.0, .25, .5, .75, 1.0].every(
          (t) =>
              (s.evaluate(pu0, pv0 + (pv1 - pv0) * t) -
                      s.evaluate(pu1, pv0 + (pv1 - pv0) * t))
                  .length <
              1e-7,
        );
    final periodicV =
        canWrap &&
        [0.0, .25, .5, .75, 1.0].every(
          (t) =>
              (s.evaluate(pu0 + (pu1 - pu0) * t, pv0) -
                      s.evaluate(pu0 + (pu1 - pu0) * t, pv1))
                  .length <
              1e-7,
        );
    double wrap(double p, double start, double end, bool periodic) =>
        periodic && (p < start || p > end)
        ? start + (p - start) % (end - start)
        : p;
    Vector3 evaluateChart(Vector3 uv) => s.evaluate(
      wrap(uv.x, pu0, pu1, periodicU),
      wrap(uv.y, pv0, pv1, periodicV),
    );
    final loops = <BrepLoop>[];
    for (final loop in face.loops) {
      final points = sampleLoop(solid, loop, tess.cylinderSegments, tess: tess),
          uv = <Vector3>[];
      for (final p in points) {
        final pair = locate(p);
        final q = Vector3(pair.$1, pair.$2, 0);
        for (var axis = 0; axis < 2; axis++) {
          final periodic = axis == 0 ? periodicU : periodicV;
          final period = axis == 0 ? pu1 - pu0 : pv1 - pv0;
          if (!periodic || uv.isEmpty) continue;
          while (q[axis] - uv.last[axis] > period * .5) q[axis] -= period;
          while (q[axis] - uv.last[axis] < -period * .5) q[axis] += period;
        }
        uv.add(q);
      }
      for (var axis = 0; axis < 2; axis++) {
        if (!(axis == 0 ? periodicU : periodicV)) continue;
        final start = axis == 0 ? pu0 : pv0, end = axis == 0 ? pu1 : pv1;
        final mean = uv.fold(0.0, (double sum, p) => sum + p[axis]) / uv.length;
        final shift =
            (((start + end) * .5 - mean) / (end - start)).round() *
            (end - start);
        for (final q in uv) q[axis] += shift;
      }
      for (var i = 0; i < uv.length; i++) world[key(uv[i])] = points[i];
      if (uv.length < 3) continue;
      final ids = uv.map(builder.vertex).toList();
      loops.add(
        BrepLoop([
          for (var i = 0; i < ids.length; i++)
            BrepCoedge(
              builder.edge(
                ids[i],
                ids[(i + 1) % ids.length],
                BrepLine(uv[i], uv[(i + 1) % uv.length]),
              ),
            ),
        ]),
      );
    }
    if (loops.isEmpty)
      throw BrepUnsupported('Surface has no usable parameter boundary');
    final plane = BrepPlane(origin: Vector3.zero(), normal: Vector3(0, 0, 1));
    builder.face(plane, loops);
    final flat = builder.build();
    final helical =
        s is BrepRuledSurface && s.lower is BrepHelix && s.upper is BrepHelix;
    // Axially clipped thread bands are convex in UV. Pair the two boundary
    // chains by angle before refinement: an ear-clipped fan otherwise spans
    // most of a turn, creating thousands of unnecessary interior triangles.
    // Every original boundary sample is retained for watertight shared edges.
    final strip = loops.length == 1
        ? _convexParameterStrip(flat.vertices)
        : null;
    final triangles = strip ?? <List<Vector3>>[];
    if (strip == null) {
      final rings = [
        for (final loop in flat.faces.single.loops)
          [
            for (final co in loop.coedges)
              co.reversed ? flat.edges[co.edge].v1 : flat.edges[co.edge].v0,
          ],
      ];
      var cells = MeshPolygon.triangulateRings(
        flat.vertices,
        rings,
        Vector3(0, 0, 1),
        keepCollinear: true,
      );
      if (cells.isEmpty) {
        cells = MeshPolygon.triangulateRings(
          flat.vertices,
          rings,
          Vector3(0, 0, -1),
          keepCollinear: true,
        );
      }
      if (cells.isEmpty)
        throw BrepUnsupported('Cannot triangulate surface parameter boundary');
      for (final cell in cells) {
        triangles.add([for (final id in cell) flat.vertices[id]]);
      }
    }
    final shortenDiagonals = helical && strip == null;
    if (shortenDiagonals) {
      _shortenParameterDiagonals(triangles);
    }
    final boundary = <String>{};
    String edgeKey(Vector3 a, Vector3 b) {
      final x = key(a), y = key(b);
      return x.compareTo(y) < 0 ? '$x/$y' : '$y/$x';
    }

    for (final loop in loops) {
      for (final co in loop.coedges) {
        final e = flat.edges[co.edge];
        boundary.add(edgeKey(flat.vertices[e.v0], flat.vertices[e.v1]));
      }
    }
    if (!helical) {
      final points = <Vector3>[], ids = <String, int>{};
      int vertex(Vector3 p) => ids.putIfAbsent(key(p), () {
        points.add(p);
        return points.length - 1;
      });
      var indices = [
        for (final t in triangles) [for (final p in t) vertex(p)],
      ]..removeWhere((t) => t[0] == t[1] || t[1] == t[2] || t[2] == t[0]);
      // Parameter welding can collapse tiny boundary slivers to a segment.
      // They carry no area; sending them into edge refinement would request
      // a third vertex from a triangle that only has two distinct vertices.
      final positions = [
        for (final p in points) world[key(p)] ?? evaluateChart(p),
      ];
      final originalCount = points.length;
      (int, int) pair(int a, int b) => a < b ? (a, b) : (b, a);
      final fixed = <(int, int)>{};
      for (final t in indices)
        for (var i = 0; i < 3; i++) {
          final a = t[i], b = t[(i + 1) % 3];
          if (boundary.contains(edgeKey(points[a], points[b])))
            fixed.add(pair(a, b));
        }
      final tolerance =
          tess.surfaceChordToleranceMm ??
          switch (tess.quality) {
            MeshQuality.draft => .02,
            MeshQuality.standard => .005,
            MeshQuality.fine => .001,
          };
      _improveChartTriangles(points, indices);
      for (var pass = 0; pass < tess.maxSurfaceRefinement; pass++) {
        if (IjoLog.verbose)
          IjoLog.line(
            'parametric face $fi refinement $pass: ${indices.length} triangles',
          );
        final splits = <(int, int), int>{};
        for (final t in indices)
          for (var i = 0; i < 3; i++) {
            final a = t[i], b = t[(i + 1) % 3], edge = pair(a, b);
            if (fixed.contains(edge) || splits.containsKey(edge)) continue;
            final uv = (points[a] + points[b]) * .5, p = evaluateChart(uv);
            if (!(a < originalCount && b < originalCount) &&
                (p - (positions[a] + positions[b]) * .5).length <= tolerance)
              continue;
            splits[edge] = points.length;
            points.add(uv);
            positions.add(p);
          }
        if (splits.isEmpty) break;
        final next = <List<int>>[];
        for (final t in indices) {
          var pieces = [t];
          for (var i = 0; i < 3; i++) {
            final a = t[i], b = t[(i + 1) % 3], mid = splits[pair(a, b)];
            if (mid == null) continue;
            final updated = <List<int>>[];
            for (final piece in pieces) {
              if (!piece.contains(a) || !piece.contains(b)) {
                updated.add(piece);
                continue;
              }
              final c = piece.firstWhere((v) => v != a && v != b);
              updated.addAll([
                [a, mid, c],
                [mid, b, c],
              ]);
            }
            pieces = updated;
          }
          next.addAll(pieces);
        }
        indices = next;
        _improveChartTriangles(points, indices);
      }
      for (final t in indices) {
        final a = positions[t[0]], b = positions[t[1]], c = positions[t[2]];
        final uvArea = (points[t[1]] - points[t[0]])
            .cross(points[t[2]] - points[t[0]])
            .z;
        if ((uvArea > 0) != face.reversed)
          emit(a, b, c, fi);
        else
          emit(a, c, b, fi);
      }
      return;
    }
    var current = triangles;
    final chordError = helical
        ? [s.lower as BrepHelix, s.upper as BrepHelix]
                  .map((c) => math.max(c.radiusAt(c.t0), c.radiusAt(c.t1)))
                  .reduce(math.max) *
              (1 - math.cos(math.pi / tess.cylinderSegments))
        : 0.0;
    Vector3 point(Vector3 uv) => world[key(uv)] ?? evaluateChart(uv);
    // A chord joining two original boundary vertices can coincide with an
    // interior chord on the neighbouring face. Give it a surface midpoint so
    // welding cannot merge unrelated diagonals across the cut boundary.
    bool needsSplit(Vector3 a, Vector3 b) =>
        !helical ||
        shortenDiagonals &&
            world.containsKey(key(a)) &&
            world.containsKey(key(b)) ||
        (point((a + b) * .5) - (point(a) + point(b)) * .5).length > chordError;
    // Split all interior edges consistently, leaving shared BREP boundary
    // samples fixed. Adjacent triangles make the same split decision.
    for (var pass = 0; pass < (helical ? 6 : 3); pass++) {
      final next = <List<Vector3>>[];
      for (final t in current) {
        final ring = <Vector3>[];
        for (var i = 0; i < 3; i++) {
          final a = t[i], b = t[(i + 1) % 3];
          ring.add(a);
          if (!boundary.contains(edgeKey(a, b)) && needsSplit(a, b))
            ring.add((a + b) * .5);
        }
        final c = (t[0] + t[1] + t[2]) / 3;
        if (helical &&
            ring.length == 3 &&
            (point(c) - ((point(t[0]) + point(t[1]) + point(t[2])) / 3))
                    .length <=
                chordError) {
          next.add(t);
          continue;
        }
        for (var i = 0; i < ring.length; i++)
          next.add([c, ring[i], ring[(i + 1) % ring.length]]);
      }
      if (next.length == current.length) break;
      current = next;
    }
    for (final t in current) {
      final a = point(t[0]), b = point(t[1]), c = point(t[2]);
      if (face.reversed) {
        emit(a, c, b, fi);
      } else {
        emit(a, b, c, fi);
      }
    }
  }

  /// Replace long angular diagonals without moving any boundary sample.
  /// A valid convex quad flip stays inside the existing trimmed domain, so this
  /// also handles concave boundaries and holes. Curvature along the straight
  /// ruling is zero; favour angular strips over an ear-clipped fan before
  /// adaptive subdivision has to refine that fan many times.
  static void _shortenParameterDiagonals(List<List<Vector3>> triangles) {
    if (triangles.length < 2) return;
    final vertices = <Vector3>[], ids = <String, int>{};
    int vertex(Vector3 p) =>
        ids.putIfAbsent('${(p.x * 1e10).round()}:${(p.y * 1e10).round()}', () {
          vertices.add(p);
          return vertices.length - 1;
        });
    final ts = [
      for (final t in triangles) [for (final p in t) vertex(p)],
    ];
    var minV = double.infinity, maxV = -double.infinity;
    for (final p in vertices) {
      minV = math.min(minV, p.y);
      maxV = math.max(maxV, p.y);
    }
    final vWeight = .001 / math.max(1e-10, maxV - minV);
    double distance(int a, int b) {
      final d = vertices[a] - vertices[b];
      return d.x * d.x + d.y * d.y * vWeight * vWeight;
    }

    double cross(int a, int b, int c) {
      final p = vertices[a], q = vertices[b], r = vertices[c];
      return (q.x - p.x) * (r.y - p.y) - (q.y - p.y) * (r.x - p.x);
    }

    List<int> oriented(int a, int b, int c) =>
        cross(a, b, c) > 0 ? [a, b, c] : [a, c, b];
    for (var pass = 0; pass < math.min(2048, vertices.length * 2); pass++) {
      final uses = <(int, int), List<int>>{};
      for (var i = 0; i < ts.length; i++) {
        final t = ts[i];
        for (var k = 0; k < 3; k++) {
          final a = t[k], b = t[(k + 1) % 3];
          (uses[a < b ? (a, b) : (b, a)] ??= []).add(i);
        }
      }
      final changed = <int>{};
      for (final e in uses.entries) {
        if (e.value.length != 2) continue;
        final i = e.value[0], j = e.value[1];
        if (changed.contains(i) || changed.contains(j)) continue;
        final a = e.key.$1, b = e.key.$2;
        final c = ts[i].firstWhere((v) => v != a && v != b, orElse: () => a);
        final d = ts[j].firstWhere((v) => v != a && v != b, orElse: () => a);
        if (c == d || distance(c, d) >= distance(a, b) * (1 - 1e-12)) continue;
        bool opposite(int p, int q, int r, int t) {
          final x = cross(p, q, r), y = cross(p, q, t);
          final scale =
              (vertices[q] - vertices[p]).length *
              math.max(
                (vertices[r] - vertices[p]).length,
                (vertices[t] - vertices[p]).length,
              );
          return x * y < 0 && math.min(x.abs(), y.abs()) > 1e-10 * scale;
        }

        if (!opposite(a, b, c, d) || !opposite(c, d, a, b)) continue;
        ts[i] = oriented(c, d, a);
        ts[j] = oriented(d, c, b);
        changed.addAll([i, j]);
      }
      if (changed.isEmpty) break;
    }
    for (var i = 0; i < ts.length; i++) {
      triangles[i] = [for (final v in ts[i]) vertices[v]];
    }
  }

  /// Triangulate a convex parameter polygon as a strip along u, keeping the
  /// sampled boundary. Nonconvex domains keep the general polygon tessellator.
  static List<List<Vector3>>? _convexParameterStrip(List<Vector3> polygon) {
    if (polygon.length < 3) return null;
    double cross(Vector3 a, Vector3 b, Vector3 c) =>
        (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
    var sign = 0;
    for (var i = 0; i < polygon.length; i++) {
      final turn = cross(
        polygon[i],
        polygon[(i + 1) % polygon.length],
        polygon[(i + 2) % polygon.length],
      );
      if (turn.abs() <= 1e-12) continue;
      final current = turn > 0 ? 1 : -1;
      if (sign != 0 && current != sign) return null;
      sign = current;
    }
    if (sign == 0) return null;
    int compare(Vector3 a, Vector3 b) {
      if ((a.x - b.x).abs() > 1e-10) return a.x.compareTo(b.x);
      return a.y.compareTo(b.y);
    }

    var first = 0, last = 0;
    for (var i = 1; i < polygon.length; i++) {
      if (compare(polygon[i], polygon[first]) < 0) first = i;
      if (compare(polygon[i], polygon[last]) > 0) last = i;
    }
    List<Vector3> chain(int step) {
      final points = <Vector3>[polygon[first]];
      for (var i = first; i != last;) {
        i = (i + step + polygon.length) % polygon.length;
        points.add(polygon[i]);
      }
      return points;
    }

    final a = chain(1), b = chain(-1);
    for (final side in [a, b]) {
      for (var i = 1; i < side.length; i++) {
        if (side[i].x < side[i - 1].x - 1e-10) return null;
      }
    }
    final triangles = <List<Vector3>>[];
    var i = 0, j = 0;
    while (i + 1 < a.length || j + 1 < b.length) {
      final p = a[i], q = b[j];
      final advanceA =
          j + 1 == b.length || i + 1 < a.length && a[i + 1].x <= b[j + 1].x;
      final r = advanceA ? a[++i] : b[++j];
      final area = cross(p, q, r);
      if (area.abs() > 1e-12) {
        triangles.add(area > 0 ? [p, q, r] : [p, r, q]);
      }
    }
    // A clipped tip can have several collinear boundary samples. Do not turn
    // roundoff-sized UV triangles into real curved flaps in 3D. Split the
    // surviving triangle edge at those samples instead, keeping shared rims.
    for (final p in polygon) {
      for (var ti = 0; ti < triangles.length; ti++) {
        final t = triangles[ti];
        if (t.any((v) => identical(v, p))) continue;
        for (var k = 0; k < 3; k++) {
          final a = t[k], b = t[(k + 1) % 3], c = t[(k + 2) % 3];
          final d = b - a;
          if (d.length2 < 1e-20) continue;
          final fraction = (p - a).dot(d) / d.length2;
          if (fraction <= 1e-10 ||
              fraction >= 1 - 1e-10 ||
              (a + d * fraction - p).length > 1e-10)
            continue;
          triangles[ti] = [a, p, c];
          triangles.add([p, b, c]);
          break;
        }
      }
    }
    // Accept the strip only if every boundary segment is present exactly once.
    final indices = Map<Vector3, int>.identity();
    for (var k = 0; k < polygon.length; k++) {
      indices[polygon[k]] = k;
    }
    final uses = <(int, int), int>{};
    for (final t in triangles) {
      for (var k = 0; k < 3; k++) {
        final a = indices[t[k]]!, b = indices[t[(k + 1) % 3]]!;
        final edge = a < b ? (a, b) : (b, a);
        uses[edge] = (uses[edge] ?? 0) + 1;
      }
    }
    for (var k = 0; k < polygon.length; k++) {
      final j = (k + 1) % polygon.length;
      final edge = k < j ? (k, j) : (j, k);
      if (uses[edge] != 1) return null;
      uses.remove(edge);
    }
    if (uses.values.any((count) => count != 2)) return null;
    return triangles.isEmpty ? null : triangles;
  }

  static void _ruled(
    BrepSolid solid,
    BrepRuledSurface surface,
    BrepFace face,
    int fi,
    MeshTessellation tess,
    void Function(Vector3, Vector3, Vector3, int) emit,
  ) {
    List<double> boundary(BrepCurve curve) {
      final temp = BrepSolid(
        vertices: [curve.evaluate(curve.t0), curve.evaluate(curve.t1)],
        edges: [BrepEdge(v0: 0, v1: 1, curve: curve)],
        faces: [],
      );
      final points = sampleLoop(
        temp,
        const BrepLoop([BrepCoedge(0)]),
        tess.cylinderSegments,
        tess: tess,
      );
      return [for (var i = 0; i <= points.length; i++) i / points.length];
    }

    final lower = boundary(surface.lower), upper = boundary(surface.upper);
    final count = math.max(4, math.max(lower.length, upper.length) - 1);
    // Boundary rows use exactly the same samples as the cap loops. Interior
    // rows may be denser; a zipper joins different counts without T-junctions.
    final middle = [for (var i = 0; i <= count; i++) i / count];
    final steps = surface.lower is BrepHelix && surface.upper is BrepHelix
        ? 1
        : 8;
    final rows = [lower, for (var v = 1; v < steps; v++) middle, upper];
    void triangle(
      double ua,
      double va,
      double ub,
      double vb,
      double uc,
      double vc,
    ) {
      final a = surface.evaluate(ua, va),
          b = surface.evaluate(ub, vb),
          c = surface.evaluate(uc, vc);
      if (face.reversed) {
        emit(a, c, b, fi);
      } else {
        emit(a, b, c, fi);
      }
    }

    for (var k = 0; k < steps; k++) {
      final a = rows[k], b = rows[k + 1], v0 = k / steps, v1 = (k + 1) / steps;
      var i = 0, j = 0;
      while (i < a.length - 1 || j < b.length - 1) {
        final na = i < a.length - 1 ? a[i + 1] : double.infinity;
        final nb = j < b.length - 1 ? b[j + 1] : double.infinity;
        if ((na - nb).abs() < 1e-12) {
          triangle(a[i], v0, na, v0, b[j], v1);
          triangle(na, v0, nb, v1, b[j], v1);
          i++;
          j++;
        } else if (na < nb) {
          triangle(a[i], v0, na, v0, b[j], v1);
          i++;
        } else {
          triangle(a[i], v0, nb, v1, b[j], v1);
          j++;
        }
      }
    }
  }

  /// Map each display segment that lies on a B-Rep curve to that edge id.
  static MeshData _tagBrepEdges(
    BrepSolid solid,
    MeshData mesh,
    int curveSegments,
  ) {
    if (mesh.faces.isEmpty || solid.edges.isEmpty) return mesh;
    final edgeIds = <String, int>{};

    List<int> candidatesForFace(int? fi) {
      if (fi == null || fi < 0 || fi >= solid.faces.length) {
        return [for (var i = 0; i < solid.edges.length; i++) i];
      }
      final ids = <int>{};
      for (final loop in solid.faces[fi].loops) {
        for (final c in loop.coedges) {
          ids.add(c.edge);
        }
      }
      return ids.isEmpty
          ? [for (var i = 0; i < solid.edges.length; i++) i]
          : ids.toList();
    }

    // Only interfaces between BREP faces (and unpaired display edges) can
    // be visible model boundaries. Internal triangles need no curve queries.
    final incidence = <(int, int), List<int>>{};
    for (var t = 0; t < mesh.faces.length; t++) {
      final f = mesh.faces[t], fi = mesh.brepFaceOfTriangle(t) ?? -1;
      for (var k = 0; k < f.length; k++) {
        final a = f[k],
            b = f[(k + 1) % f.length],
            key = a < b ? (a, b) : (b, a);
        (incidence[key] ??= []).add(fi);
      }
    }
    // These IDs label display chords. Reuse the shared sampling used to build
    // those chords instead of repeatedly solving a global nearest point on
    // the native implicit edge for every boundary triangle vertex.
    final displayCurves = <int, List<Vector3>>{};
    bool traced(BrepCurve c) =>
        c is BrepSurfaceIntersectionCurve ||
        c is BrepSplineCurve ||
        c is BrepAffineCurve && traced(c.base);
    double displayDistance(int edge, Vector3 p) {
      final samples = displayCurves.putIfAbsent(edge, () {
        final c = solid.edges[edge].curve,
            n = math.max(
              c is BrepSurfaceIntersectionCurve
                  ? 32
                  : c is BrepSplineCurve
                  ? 16
                  : 8,
              curveSegments,
            );
        return [
          for (var i = 0; i <= n; i++) c.evaluate(c.t0 + (c.t1 - c.t0) * i / n),
        ];
      });
      var best = double.infinity;
      for (var i = 1; i < samples.length; i++) {
        final a = samples[i - 1], d = samples[i] - a;
        final t = d.length2 == 0
            ? 0.0
            : ((p - a).dot(d) / d.length2).clamp(0.0, 1.0);
        best = math.min(best, (p - a - d * t).length2);
      }
      return math.sqrt(best);
    }

    int? matchEdge(Vector3 a, Vector3 b, List<int> candidates) {
      var best = -1;
      var bestD = double.infinity;
      for (final i in candidates) {
        if (i < 0 || i >= solid.edges.length) continue;
        final c = solid.edges[i].curve;
        if (c is BrepPlaneSectionCurve) {
          final tolerance = math.max(1e-3, (a - b).length * .02);
          if (c.plane.height(a).abs() > tolerance ||
              c.plane.height(b).abs() > tolerance)
            continue;
        }
        if (c is BrepIntersectionCurve) {
          final tolerance = math.max(1e-3, (a - b).length * .02);
          if ((c.target.closestPoint(a) - a).length > tolerance ||
              (c.target.closestPoint(b) - b).length > tolerance ||
              (c.carrier.closestPoint(a) - a).length > tolerance ||
              (c.carrier.closestPoint(b) - b).length > tolerance)
            continue;
        }
        final d = traced(c)
            ? displayDistance(i, a) + displayDistance(i, b)
            : (c.closestPoint(a) - a).length + (c.closestPoint(b) - b).length;
        if (d < bestD) {
          bestD = d;
          best = i;
        }
      }
      final chord = (a - b).length;
      if (best >= 0 && bestD < math.max(1e-3, chord * 0.02)) return best;
      return null;
    }

    for (var t = 0; t < mesh.faces.length; t++) {
      final face = mesh.faces[t];
      final fi = mesh.brepFaceOfTriangle(t);
      final candidates = candidatesForFace(fi);
      for (var i = 0; i < face.length; i++) {
        final ia = face[i];
        final ib = face[(i + 1) % face.length];
        if (ia == ib) continue;
        final adjacent = incidence[ia < ib ? (ia, ib) : (ib, ia)]!;
        if (adjacent.length == 2 && adjacent[0] == adjacent[1]) continue;
        final a = mesh.vertices[ia];
        final b = mesh.vertices[ib];
        final key = MeshData.segmentKey(a, b);
        if (edgeIds.containsKey(key)) continue;
        final id = matchEdge(a, b, candidates);
        if (id != null) edgeIds[key] = id;
      }
    }
    if (edgeIds.isEmpty) return mesh;
    return MeshData(
      vertices: mesh.vertices,
      faces: mesh.faces,
      faceIds: mesh.faceIds,
      edgeIds: edgeIds,
    );
  }

  static MeshData _weld(MeshData mesh, {double eps2 = 1e-10}) {
    final verts = <Vector3>[];
    final cellSize = math.sqrt(eps2);
    final cells = <(int, int, int), List<int>>{};
    int idx(Vector3 p) {
      final x = (p.x / cellSize).floor();
      final y = (p.y / cellSize).floor();
      final z = (p.z / cellSize).floor();
      var match = verts.length;
      for (var dx = -1; dx <= 1; dx++) {
        for (var dy = -1; dy <= 1; dy++) {
          for (var dz = -1; dz <= 1; dz++) {
            for (final i in cells[(x + dx, y + dy, z + dz)] ?? const <int>[]) {
              if (i < match && (verts[i] - p).length2 < eps2) match = i;
            }
          }
        }
      }
      // Keep the same first representative as the exhaustive search. Check
      // neighboring cells so points near a cell boundary still weld.
      if (match < verts.length) return match;
      (cells[(x, y, z)] ??= []).add(verts.length);
      verts.add(p.clone());
      return verts.length - 1;
    }

    final faces = <List<int>>[];
    final faceIds = mesh.faceIds == null ? null : <int>[];
    for (var f = 0; f < mesh.faces.length; f++) {
      final mapped = [for (final i in mesh.faces[f]) idx(mesh.vertices[i])];
      if (mapped.length < 3) continue;
      if (mapped[0] == mapped[1] ||
          mapped[1] == mapped[2] ||
          mapped[2] == mapped[0]) {
        continue;
      }
      faces.add(mapped);
      faceIds?.add(mesh.faceIds![f]);
    }
    return MeshData(
      vertices: verts,
      faces: faces,
      faceIds: faceIds,
      edgeIds: mesh.edgeIds,
    );
  }

  static int _splineSamples(
    BrepSplineCurve curve,
    int circleSegs,
    MeshTessellation? tess,
  ) {
    // Polynomial cubic Bezier spans have a control-net bound on their second
    // derivative. Aligned knot breaks ensure no display chord crosses a kink.
    // This avoids sampling every short fitted rim span like an entire circle.
    if (tess != null &&
        curve.degree == 3 &&
        curve.points.length % 3 == 1 &&
        curve.weights.every((w) => (w - 1).abs() < 1e-12) &&
        curve.t0 == curve.knots.first &&
        curve.t1 == curve.knots.last) {
      final breaks = curve.knots.toSet().toList()..sort();
      final spans = (curve.points.length - 1) ~/ 3;
      if (breaks.length == spans + 1 &&
          breaks.asMap().entries.every(
            (e) =>
                curve.knots.where((t) => t == e.value).length ==
                (e.key == 0 || e.key == spans ? 4 : 3),
          )) {
        var secondDerivative = 0.0;
        for (var i = 0; i < spans; i++) {
          final width = (breaks[i + 1] - breaks[i]) / (curve.t1 - curve.t0);
          for (var j = 0; j < 2; j++) {
            final k = i * 3 + j;
            secondDerivative = math.max(
              secondDerivative,
              6 *
                  (curve.points[k + 2] -
                          curve.points[k + 1] * 2 +
                          curve.points[k])
                      .length /
                  (width * width),
            );
          }
        }
        final error = switch (tess.quality) {
          MeshQuality.draft => .02,
          MeshQuality.standard => .005,
          MeshQuality.fine => .001,
        };
        for (
          var n = secondDerivative < 1e-8 ? 1 : 4;
          n <= math.max(16, circleSegs);
          n *= 2
        ) {
          if (secondDerivative / (8 * n * n) <= error &&
              breaks.every((t) {
                final q = (t - curve.t0) / (curve.t1 - curve.t0) * n;
                return (q - q.round()).abs() < 1e-9;
              }))
            return n;
        }
      }
    }
    // A rational quadratic quarter-circle/ellipse is the cross-section of a
    // rolling fillet, even though it is stored as a spline. Match the ordinary
    // circular-arc display density instead of sampling every small arc 32 times.
    if (tess != null &&
        curve.degree == 2 &&
        curve.points.length == 3 &&
        (curve.weights[0] - curve.weights[2]).abs() < 1e-12 &&
        (curve.weights[1] / curve.weights[0] - math.sqrt(.5)).abs() < 1e-12 &&
        curve.t0 == curve.knots.first &&
        curve.t1 == curve.knots.last) {
      final a = curve.points[0],
          b = curve.points[2],
          center = a + b - curve.points[1];
      final x = a - center, y = b - center, r = math.max(x.length, y.length);
      if (r > 1e-9 &&
          math.min(x.length, y.length) > 1e-9 &&
          x.dot(y).abs() < 1e-7 * r * r) {
        return tess.filletArcSegments(r, math.pi / 2);
      }
    }
    return math.max(16, circleSegs);
  }

  static bool _isLinearCurve(BrepCurve curve) =>
      curve is BrepLine ||
      curve is BrepAffineCurve && _isLinearCurve(curve.base);

  static List<Vector3> sampleLoop(
    BrepSolid solid,
    BrepLoop loop,
    int circleSegs, {
    MeshTessellation? tess,
  }) {
    final pts = <Vector3>[];
    for (final c in loop.coedges) {
      final e = solid.edges[c.edge];
      final curve = e.curve;
      final n = switch (curve) {
        BrepHelix() => math.max(
          4,
          ((curve.t1 - curve.t0).abs() * circleSegs / (2 * math.pi)).ceil(),
        ),
        BrepSurfaceIntersectionCurve() => math.max(32, circleSegs),
        BrepRuledQuadricCurve() => math.max(32, circleSegs),
        BrepPlaneSectionCurve() => math.max(16, circleSegs),
        BrepSplineCurve() => _splineSamples(curve, circleSegs, tess),
        BrepAffineCurve() =>
          _isLinearCurve(curve) ? 1 : math.max(8, circleSegs),
        BrepLine() => 1,
        BrepIntersectionCurve() => math.max(
          8,
          ((curve.t1 - curve.t0).abs() * circleSegs / math.pi).ceil(),
        ),
        BrepCircle() =>
          curve.isClosed
              ? (tess?.cylinderSegments ?? circleSegs)
              : (tess != null
                    ? tess.filletArcSegments(
                        curve.radius,
                        (curve.t1 - curve.t0).abs(),
                      )
                    : math.max(
                        4,
                        (circleSegs *
                                (curve.t1 - curve.t0).abs() /
                                (2 * math.pi))
                            .round(),
                      )),
        BrepEllipse() =>
          curve.isClosed
              ? (tess?.cylinderSegments ?? circleSegs)
              : math.max(8, circleSegs ~/ 2),
      };
      final start = solid.vertices[c.reversed ? e.v1 : e.v0];
      if (curve is BrepCircle && curve.isClosed) {
        pts.addAll(_sampleClosedCircle(c, curve, n));
        continue;
      }
      if (curve is BrepEllipse && curve.isClosed) {
        pts.addAll(_sampleClosedEllipse(c, curve, n));
        continue;
      }
      final a0 = curve.evaluate(curve.t0);
      final a1 = curve.evaluate(curve.t1);
      // A closed spline/affine curve has coincident endpoints. Distances
      // cannot identify its traversal direction; use the oriented coedge.
      final forward = e.v0 == e.v1 && (a0 - a1).length < 1e-8
          ? !c.reversed
          : (start - a0).length <= (start - a1).length;
      for (var i = 0; i < n; i++) {
        if (i == 0) {
          pts.add(start.clone());
          continue;
        }
        final frac = i / n;
        final t = forward
            ? curve.t0 + (curve.t1 - curve.t0) * frac
            : curve.t1 - (curve.t1 - curve.t0) * frac;
        pts.add(curve.evaluate(t));
      }
    }
    return pts;
  }

  static List<Vector3> _sampleClosedCircle(
    BrepCoedge c,
    BrepCircle curve,
    int n,
  ) {
    return [
      for (var i = 0; i < n; i++)
        curve.evaluate(
          c.reversed
              ? curve.t1 - (curve.t1 - curve.t0) * (i / n)
              : curve.t0 + (curve.t1 - curve.t0) * (i / n),
        ),
    ];
  }

  static List<Vector3> _sampleClosedEllipse(
    BrepCoedge c,
    BrepEllipse curve,
    int n,
  ) {
    return [
      for (var i = 0; i < n; i++)
        curve.evaluate(
          c.reversed
              ? curve.t1 - (curve.t1 - curve.t0) * (i / n)
              : curve.t0 + (curve.t1 - curve.t0) * (i / n),
        ),
    ];
  }

  static void _plane(
    BrepSolid solid,
    BrepFace face,
    int fi,
    MeshTessellation tess,
    void Function(Vector3, Vector3, Vector3, int) addTri,
  ) {
    if (face.loops.isEmpty) return;
    final segs = tess.cylinderSegments;
    final outer = sampleLoop(solid, face.loops.first, segs, tess: tess);
    if (outer.length < 3) return;
    final plane = face.surface as BrepPlane;
    if (face.loops.length > 1) {
      final holes = [
        for (var i = 1; i < face.loops.length; i++)
          sampleLoop(solid, face.loops[i], segs, tess: tess),
      ]..sort((a, b) => a.length.compareTo(b.length));
      // Use the same complete, boundary-preserving triangulation as planar
      // mesh profiles. The legacy hole bridge/ear path can return partial
      // triangles for several aligned holes, leaving an otherwise closed
      // BREP with gaps in its display mesh.
      final points = <Vector3>[];
      final rings = <List<int>>[];
      for (final loop in [outer, ...holes]) {
        final start = points.length;
        points.addAll(loop);
        rings.add(List.generate(loop.length, (i) => start + i));
      }
      for (final tri in MeshPolygon.triangulateRings(
        points,
        rings,
        face.reversed ? -plane.normal : plane.normal,
        keepCollinear: true,
      )) {
        addTri(points[tri[0]], points[tri[1]], points[tri[2]], fi);
      }
      return;
    }
    // A hole can be encoded in the outer wire using a bridge traversed
    // both ways. Preserve those repeated junctions during ear clipping;
    // deleting them can stall clipping and fan triangles across the hole.
    final edgeUses = <int, int>{};
    for (final c in face.loops.first.coedges) {
      edgeUses[c.edge] = (edgeUses[c.edge] ?? 0) + 1;
    }
    if (edgeUses.values.any((count) => count > 1)) {
      final tris = MeshPolygon.triangulateRings(
        outer,
        [List.generate(outer.length, (i) => i)],
        face.reversed ? -plane.normal : plane.normal,
        keepCollinear: true,
      );
      if (tris.isNotEmpty) {
        for (final t in tris) {
          addTri(outer[t[0]], outer[t[1]], outer[t[2]], fi);
        }
        return;
      }
    }
    // Boolean cuts retain collinear vertices shared with adjacent faces.
    // Remove them from the ear search and restore them on triangle edges;
    // otherwise a thin concave cap can stall and fall back to an invalid fan.
    _earClipPlane(
      outer,
      plane,
      face.reversed,
      fi,
      addTri,
      preserveCollinear: true,
      pruneStraightBoundary: face.loops.single.coedges.every(
        (co) => solid.edges[co.edge].curve is BrepLine,
      ),
    );
  }

  static void _centroidFanPlane(
    List<Vector3> outer,
    BrepPlane plane,
    bool reversed,
    int fi,
    void Function(Vector3, Vector3, Vector3, int) addTri,
  ) {
    final center = outer.reduce((a, b) => a + b) / outer.length.toDouble();
    for (var i = 0; i < outer.length; i++) {
      final a = outer[i], b = outer[(i + 1) % outer.length];
      var n = (a - center).cross(b - center);
      if (reversed) n = -n;
      if (n.dot(plane.normal) < 0) {
        addTri(center, b, a, fi);
      } else {
        addTri(center, a, b, fi);
      }
    }
  }

  static void _earClipPlane(
    List<Vector3> outer,
    BrepPlane plane,
    bool reversed,
    int fi,
    void Function(Vector3, Vector3, Vector3, int) addTri, {
    bool centroidFallback = true,
    bool preserveCollinear = false,
    bool pruneStraightBoundary = false,
  }) {
    final pts = List<Vector3>.from(outer);
    if (pts.length < 3) return;
    final want = reversed ? -plane.normal : plane.normal;
    final tris = <(Vector3, Vector3, Vector3)>[];
    void emit(Vector3 a, Vector3 b, Vector3 c) {
      var n = (b - a).cross(c - a);
      if (n.length2 < 1e-20) return;
      if (n.dot(want) < 0) {
        tris.add((a, c, b));
      } else {
        tris.add((a, b, c));
      }
    }

    if (pts.length == 3) {
      emit(pts[0], pts[1], pts[2]);
      for (final t in tris) {
        addTri(t.$1, t.$2, t.$3, fi);
      }
      return;
    }

    final origin = pts.first;
    var u = plane.xDir.clone();
    var v = plane.yDir.clone();
    if (u.length2 < 1e-16 || v.length2 < 1e-16) {
      u = plane.normal.cross(Vector3(1, 0, 0));
      if (u.length2 < 1e-12) u = plane.normal.cross(Vector3(0, 1, 0));
      u.normalize();
      v = plane.normal.cross(u).normalized();
    }

    var ring = [
      for (final p in pts) Vector2((p - origin).dot(u), (p - origin).dot(v)),
    ];
    var idx = List<int>.generate(pts.length, (i) => i);
    var area = 0.0;
    for (var i = 0; i < ring.length; i++) {
      final a = ring[i];
      final b = ring[(i + 1) % ring.length];
      area += a.x * b.y - b.x * a.y;
    }
    if (area < 0) {
      ring = ring.reversed.toList();
      idx = idx.reversed.toList();
    }
    final origRing = [for (final p in ring) p.clone()];

    bool evenOdd(Vector2 p) {
      var inside = false;
      for (var i = 0, j = origRing.length - 1; i < origRing.length; j = i++) {
        final yi = origRing[i].y;
        final yj = origRing[j].y;
        if ((yi > p.y) == (yj > p.y)) continue;
        final xi = origRing[i].x;
        final xj = origRing[j].x;
        final xHit = (xj - xi) * (p.y - yi) / (yj - yi) + xi;
        if (p.x < xHit) inside = !inside;
      }
      return inside;
    }

    bool convex(int i0, int i1, int i2) {
      final a = ring[i0];
      final b = ring[i1];
      final c = ring[i2];
      return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x) > 1e-12;
    }

    bool pointIn(Vector2 p, Vector2 a, Vector2 b, Vector2 c) {
      final c1 = (b.x - a.x) * (p.y - a.y) - (b.y - a.y) * (p.x - a.x);
      final c2 = (c.x - b.x) * (p.y - b.y) - (c.y - b.y) * (p.x - b.x);
      final c3 = (a.x - c.x) * (p.y - c.y) - (a.y - c.y) * (p.x - c.x);
      final hasNeg = c1 < -1e-12 || c2 < -1e-12 || c3 < -1e-12;
      final hasPos = c1 > 1e-12 || c2 > 1e-12 || c3 > 1e-12;
      return !(hasNeg && hasPos);
    }

    final active = List<int>.generate(ring.length, (i) => i);
    final collinear = <Vector3>[];
    // Remove straight boundary subdivisions before clipping any ears.
    // Otherwise clipping the opposite corners first can leave a long, flat
    // final triangle; restoring its samples then creates overlapping fans.
    // Limit this to linear trims: nearly straight spline samples still carry
    // curvature and must retain their original ear selection order.
    if (pruneStraightBoundary) {
      var changed = true;
      while (changed && active.length > 3) {
        changed = false;
        for (var i = 0; i < active.length; i++) {
          final i0 = active[(i + active.length - 1) % active.length];
          final i1 = active[i];
          final i2 = active[(i + 1) % active.length];
          final a = ring[i0], b = ring[i1], c = ring[i2];
          final d = c - a;
          if (d.length2 < 1e-24) continue;
          final t = (b - a).dot(d) / d.length2;
          if (t >= 0 && t <= 1 && (a + d * t - b).length2 < 1e-18) {
            collinear.add(pts[idx[i1]]);
            active.removeAt(i);
            changed = true;
            break;
          }
        }
      }
    }
    var guard = 0;
    final maxGuard = active.length * active.length * 4;
    while (active.length > 3 && guard < maxGuard) {
      guard++;
      var clipped = false;
      for (var i = 0; i < active.length; i++) {
        final i0 = active[(i + active.length - 1) % active.length];
        final i1 = active[i];
        final i2 = active[(i + 1) % active.length];
        final a = ring[i0];
        final b = ring[i1];
        final c = ring[i2];
        final ac = c - a;
        final along = ac.length2 < 1e-24
            ? double.nan
            : (b - a).dot(ac) / ac.length2;
        if (preserveCollinear &&
            along.isFinite &&
            along >= -1e-12 &&
            along <= 1 + 1e-12 &&
            (a + ac * along - b).length2 < 1e-20) {
          collinear.add(pts[idx[i1]]);
          active.removeAt(i);
          clipped = true;
          break;
        }
        if (!convex(i0, i1, i2)) continue;
        if ((a - b).length2 < 1e-20 ||
            (b - c).length2 < 1e-20 ||
            (c - a).length2 < 1e-20) {
          active.removeAt(i);
          clipped = true;
          break;
        }
        final centroid = Vector2((a.x + b.x + c.x) / 3, (a.y + b.y + c.y) / 3);
        if (!evenOdd(centroid)) continue;
        var inside = false;
        for (final j in active) {
          if (j == i0 || j == i1 || j == i2) continue;
          final p = ring[j];
          if ((p - ring[i0]).length2 < 1e-20 ||
              (p - ring[i1]).length2 < 1e-20 ||
              (p - ring[i2]).length2 < 1e-20) {
            continue;
          }
          if (pointIn(p, ring[i0], ring[i1], ring[i2])) {
            inside = true;
            break;
          }
        }
        if (inside) continue;
        emit(pts[idx[i0]], pts[idx[i1]], pts[idx[i2]]);
        active.removeAt(i);
        clipped = true;
        break;
      }
      if (!clipped) break;
    }
    if (active.length >= 3) {
      emit(pts[idx[active[0]]], pts[idx[active[1]]], pts[idx[active[2]]]);
    }
    if (tris.length < pts.length - collinear.length - 2 && centroidFallback) {
      // Roundoff after a rotated boolean can stall the fixed-tolerance ear
      // search on a concave cap. Retry with scale-aware predicates before
      // considering a fan, which can cross the recess or overlap itself.
      final robust = MeshPolygon.triangulateRings(
        pts,
        [List.generate(pts.length, (i) => i)],
        want,
        keepCollinear: true,
      );
      if (robust.isNotEmpty) {
        for (final t in robust) {
          addTri(pts[t[0]], pts[t[1]], pts[t[2]], fi);
        }
        return;
      }
      _centroidFanPlane(outer, plane, reversed, fi, addTri);
      return;
    }
    // Restore collinear boundary samples after clipping. Removing them only
    // from the ear search avoids deadlocks without introducing T-junctions
    // against the adjoining BREP face.
    var restored = tris;
    for (final point in collinear) {
      final next = <(Vector3, Vector3, Vector3)>[];
      for (final tri in restored) {
        final vertices = [tri.$1, tri.$2, tri.$3];
        var split = false;
        for (var i = 0; i < 3; i++) {
          final a = vertices[i],
              b = vertices[(i + 1) % 3],
              c = vertices[(i + 2) % 3];
          final d = b - a;
          if (d.length2 < 1e-24) continue;
          final t = (point - a).dot(d) / d.length2;
          if (t <= 1e-9 || t >= 1 - 1e-9 || (a + d * t - point).length2 > 1e-18)
            continue;
          next.add((a, point, c));
          next.add((point, b, c));
          split = true;
          break;
        }
        if (!split) next.add(tri);
      }
      restored = next;
    }
    for (final t in restored) {
      addTri(t.$1, t.$2, t.$3, fi);
    }
  }

  static void _planeWithHoles(
    List<Vector3> outer,
    List<List<Vector3>> holes,
    BrepPlane plane,
    bool reversed,
    int fi,
    void Function(Vector3, Vector3, Vector3, int) addTri, {
    bool preferOuterBridges = false,
  }) {
    if (holes.isEmpty) {
      _earClipPlane(
        outer,
        plane,
        reversed,
        fi,
        addTri,
        preserveCollinear: preferOuterBridges,
      );
      return;
    }
    final origin = outer.first;
    var u = plane.xDir.clone();
    var v = plane.yDir.clone();
    if (u.length2 < 1e-16 || v.length2 < 1e-16) {
      u = plane.normal.cross(Vector3(1, 0, 0));
      if (u.length2 < 1e-12) u = plane.normal.cross(Vector3(0, 1, 0));
      u.normalize();
      v = plane.normal.cross(u).normalized();
    }
    Vector2 uv(Vector3 p) => Vector2((p - origin).dot(u), (p - origin).dot(v));

    double signedArea(List<Vector2> ring) {
      var a = 0.0;
      for (var i = 0; i < ring.length; i++) {
        final p = ring[i];
        final q = ring[(i + 1) % ring.length];
        a += p.x * q.y - q.x * p.y;
      }
      return a;
    }

    bool crosses(Vector2 a, Vector2 b, Vector2 c, Vector2 d) {
      final r = b - a;
      final s = d - c;
      final den = r.x * s.y - r.y * s.x;
      if (den.abs() < 1e-14) return false;
      final t = ((c.x - a.x) * s.y - (c.y - a.y) * s.x) / den;
      final uu = ((c.x - a.x) * r.y - (c.y - a.y) * r.x) / den;
      return t > 1e-6 && t < 1 - 1e-6 && uu > 1e-6 && uu < 1 - 1e-6;
    }

    var pts = List<Vector3>.from(outer);
    var ring = [for (final p in pts) uv(p)];
    if (signedArea(ring) < 0) {
      pts = pts.reversed.toList();
      ring = ring.reversed.toList();
    }

    for (final rawHole in holes) {
      if (rawHole.length < 3) continue;
      var holePts = List<Vector3>.from(rawHole);
      var holeUv = [for (final p in holePts) uv(p)];
      if (signedArea(holeUv) > 0) {
        holePts = holePts.reversed.toList();
        holeUv = holeUv.reversed.toList();
      }
      var bestOi = -1;
      var bestHi = -1;
      var bestD = double.infinity;
      bool insideRing(Vector2 p) {
        var inside = false;
        for (var i = 0, j = ring.length - 1; i < ring.length; j = i++) {
          final yi = ring[i].y;
          final yj = ring[j].y;
          if ((yi > p.y) == (yj > p.y)) continue;
          final xi = ring[i].x;
          final xj = ring[j].x;
          final xHit = (xj - xi) * (p.y - yi) / (yj - yi) + xi;
          if (p.x < xHit) inside = !inside;
        }
        return inside;
      }

      bool blocked(int oi, int hi) {
        final a = ring[oi];
        final b = holeUv[hi];
        final mid = Vector2((a.x + b.x) * 0.5, (a.y + b.y) * 0.5);
        if (!insideRing(mid)) return true;
        for (var i = 0; i < ring.length; i++) {
          if (i == oi || (i + 1) % ring.length == oi) continue;
          if (crosses(a, b, ring[i], ring[(i + 1) % ring.length])) {
            return true;
          }
        }
        for (var i = 0; i < holeUv.length; i++) {
          if (i == hi || (i + 1) % holeUv.length == hi) continue;
          if (crosses(a, b, holeUv[i], holeUv[(i + 1) % holeUv.length])) {
            return true;
          }
        }
        return false;
      }

      // Prefer a bridge to the original outer boundary. Chaining through
      // an earlier hole creates coincident bridge junctions when holes align.
      // Retain the general visibility search when no outer bridge exists.
      for (
        var pass = preferOuterBridges ? 0 : 1;
        pass < 2 && bestOi < 0;
        pass++
      ) {
        for (var oi = 0; oi < ring.length; oi++) {
          if (pass == 0 && !outer.contains(pts[oi])) continue;
          for (var hi = 0; hi < holeUv.length; hi++) {
            final d = (ring[oi] - holeUv[hi]).length2;
            if (preferOuterBridges
                ? bestD.isFinite && d >= bestD - 1e-10 * math.max(1.0, bestD)
                : d >= bestD)
              continue;
            if (blocked(oi, hi)) continue;
            bestD = d;
            bestOi = oi;
            bestHi = hi;
          }
        }
      }
      if (bestOi < 0) continue;
      final mergedPts = <Vector3>[
        ...pts.sublist(0, bestOi + 1),
        for (var k = 0; k <= holePts.length; k++)
          holePts[(bestHi + k) % holePts.length],
        pts[bestOi],
        ...pts.sublist(bestOi + 1),
      ];
      pts = mergedPts;
      ring = [for (final p in pts) uv(p)];
    }

    _earClipPlane(
      pts,
      plane,
      reversed,
      fi,
      addTri,
      centroidFallback: false,
      preserveCollinear: preferOuterBridges,
    );
  }

  /// Smallest interval on the circle that covers [raw] (may wrap past 2π).
  static (double, double) _smallestUCover(List<double> raw) {
    const twoPi = math.pi * 2;
    if (raw.isEmpty) return (0.0, twoPi);
    final us = [for (final u in raw) (u % twoPi + twoPi) % twoPi]..sort();
    final pts = <double>[us.first];
    for (var i = 1; i < us.length; i++) {
      if (us[i] - pts.last > 1e-9) pts.add(us[i]);
    }
    if (pts.length == 1) return (pts.first, pts.first);
    var maxGap = (pts.first + twoPi) - pts.last;
    var start = pts.first;
    var end = pts.last;
    for (var i = 0; i < pts.length - 1; i++) {
      final gap = pts[i + 1] - pts[i];
      if (gap > maxGap) {
        maxGap = gap;
        start = pts[i + 1];
        end = pts[i] + twoPi;
      }
    }
    return (start, end);
  }

  /// UV box covering [face] on its cylinder. [u1] may exceed 2π if the
  /// patch wraps the seam. Closed walls use a full 0..2π span.
  static (double, double, double, double) cylinderFaceDomain(
    BrepSolid solid,
    BrepFace face,
  ) {
    return _axialFaceDomain(
      solid,
      face,
      (face.surface as BrepCylinder).origin,
      (face.surface as BrepCylinder).axis,
    );
  }

  /// Angular and axial bounds of a trimmed cone, preserving circular arc spans.
  static (double, double, double, double) coneFaceDomain(
    BrepSolid solid,
    BrepFace face,
  ) => _axialFaceDomain(
    solid,
    face,
    (face.surface as BrepCone).origin,
    (face.surface as BrepCone).axis,
  );

  static (double, double, double, double) _axialFaceDomain(
    BrepSolid solid,
    BrepFace face,
    Vector3 origin,
    Vector3 axis,
  ) {
    final surface = face.surface;
    var z0 = 0.0;
    var z1 = 0.0;
    var got = false;
    final uSamples = <double>[];
    var closedCircle = false;
    final arcDomains = <(double, double)>[];
    void accZ(double z) {
      if (!got) {
        z0 = z1 = z;
        got = true;
      } else {
        z0 = math.min(z0, z);
        z1 = math.max(z1, z);
      }
    }

    for (final loop in face.loops) {
      for (final c in loop.coedges) {
        final e = solid.edges[c.edge];
        for (final p in [solid.vertices[e.v0], solid.vertices[e.v1]]) {
          accZ((p - origin).dot(axis));
          final uv = surface.parameters(p);
          if (uv != null) uSamples.add(uv.$1);
        }
        final curve = e.curve;
        if (curve is BrepLine) {
          for (final p in [curve.a, curve.b]) {
            accZ((p - origin).dot(axis));
          }
        } else if (curve is BrepCircle) {
          if (curve.isClosed) closedCircle = true;
          if (!curve.isClosed && curve.normal.cross(axis).length < 1e-8) {
            final forward = curve.normal.dot(axis) * (curve.t1 - curve.t0) > 0;
            final p = curve.evaluate(forward ? curve.t0 : curve.t1);
            final uv = surface.parameters(p);
            if (uv != null)
              arcDomains.add((uv.$1, (curve.t1 - curve.t0).abs()));
          }
          final ts = curve.isClosed
              ? [
                  curve.t0,
                  curve.t0 + math.pi / 2,
                  curve.t0 + math.pi,
                  curve.t0 + 1.5 * math.pi,
                ]
              : [curve.t0, 0.5 * (curve.t0 + curve.t1), curve.t1];
          for (final t in ts) {
            final p = curve.evaluate(t);
            accZ((p - origin).dot(axis));
            final uv = surface.parameters(p);
            if (uv != null) uSamples.add(uv.$1);
          }
        }
      }
    }
    if (!got || (z1 - z0).abs() < 1e-9) {
      z0 = 0;
      z1 = surface is BrepCone ? surface.height : 1;
    }
    if (closedCircle || uSamples.isEmpty) {
      return (0.0, math.pi * 2, z0, z1);
    }
    // Three samples cannot distinguish a long arc from its complementary
    // sector. Use the exact span for matching circular cap arcs.
    if (arcDomains.isNotEmpty) {
      final first = arcDomains.first;
      bool same((double, double) arc) {
        final delta = arc.$1 - first.$1;
        return math.sin(delta).abs() < 1e-8 &&
            math.cos(delta) > 0 &&
            (arc.$2 - first.$2).abs() < 1e-8;
      }

      if (arcDomains.every(same))
        return (first.$1, first.$1 + first.$2, z0, z1);
    }
    final cover = _smallestUCover(uSamples);
    var u0 = cover.$1;
    var u1 = cover.$2;
    if (u1 <= u0 + 1e-9) u1 += math.pi * 2;
    return (u0, u1, z0, z1);
  }

  static void _cylinder(
    BrepSolid solid,
    BrepFace face,
    int fi,
    MeshTessellation tess,
    void Function(Vector3, Vector3, Vector3, int) addTri,
  ) {
    final cyl = face.surface as BrepCylinder;
    if (face.loops.length == 2 &&
        face.loops.every((l) => l.coedges.length == 1) &&
        _cylinderFromClosedRims(solid, face, cyl, tess, fi, addTri))
      return;
    if (face.loops.any(
      (l) => l.coedges.any(
        (c) =>
            solid.edges[c.edge].curve is BrepSurfaceIntersectionCurve ||
            solid.edges[c.edge].curve is BrepRuledQuadricCurve ||
            solid.edges[c.edge].curve is BrepPlaneSectionCurve ||
            solid.edges[c.edge].curve is BrepSplineCurve ||
            solid.edges[c.edge].curve is BrepAffineCurve ||
            solid.edges[c.edge].curve is BrepIntersectionCurve ||
            solid.edges[c.edge].curve is BrepEllipse,
      ),
    )) {
      _axialTrimmed(solid, face, fi, tess, addTri);
      return;
    }
    var z0 = 0.0;
    var z1 = 0.0;
    var got = false;
    final uSamples = <double>[];
    var closedCircle = false;
    var nLine = 0;
    var nCircle = 0;
    void accZ(double z) {
      if (!got) {
        z0 = z1 = z;
        got = true;
      } else {
        z0 = math.min(z0, z);
        z1 = math.max(z1, z);
      }
    }

    void accU(double u) => uSamples.add(u);

    for (final loop in face.loops) {
      for (final c in loop.coedges) {
        final e = solid.edges[c.edge];
        for (final p in [solid.vertices[e.v0], solid.vertices[e.v1]]) {
          accZ((p - cyl.origin).dot(cyl.axis));
          final uv = cyl.parameters(p);
          if (uv != null) accU(uv.$1);
        }
        final curve = e.curve;
        if (curve is BrepLine) {
          nLine++;
          for (final p in [curve.a, curve.b]) {
            accZ((p - cyl.origin).dot(cyl.axis));
          }
        } else if (curve is BrepCircle) {
          nCircle++;
          if (curve.isClosed) {
            closedCircle = true;
          }
          // Circle t is in the circle's own frame, not cylinder U. Sample
          // 3D points and convert so a 90° fillet is not meshed as 180/270°.
          final ts = curve.isClosed
              ? [
                  curve.t0,
                  curve.t0 + math.pi / 2,
                  curve.t0 + math.pi,
                  curve.t0 + 1.5 * math.pi,
                ]
              : [curve.t0, 0.5 * (curve.t0 + curve.t1), curve.t1];
          for (final t in ts) {
            final p = curve.evaluate(t);
            accZ((p - cyl.origin).dot(cyl.axis));
            final uv = cyl.parameters(p);
            if (uv != null) accU(uv.$1);
          }
        }
      }
    }
    if (!got || (z1 - z0).abs() < 1e-9) {
      z0 = 0;
      z1 = 1;
    }
    final uvRect =
        !closedCircle &&
        face.loops.length == 1 &&
        face.loops.first.coedges.length == 4 &&
        nLine == 2 &&
        nCircle == 2;
    if (uvRect && _cylinderPatchFromArcs(solid, face, cyl, tess, fi, addTri)) {
      return;
    }
    if (closedCircle &&
        _cylinderFromClosedRims(solid, face, cyl, tess, fi, addTri)) {
      return;
    }
    if (!closedCircle && !uvRect && face.loops.isNotEmpty) {
      final outer = sampleLoop(
        solid,
        face.loops.first,
        tess.cylinderSegments,
        tess: tess,
      );
      if (outer.length >= 3) {
        final uvs = axialLoopParameters(cyl, outer);
        final consistentTrims = face.loops.every(
          (loop) => loop.coedges.every((co) {
            final edge = solid.edges[co.edge], curve = edge.curve;
            return (curve.evaluate(curve.t0) - solid.vertices[edge.v0]).length <
                    1e-6 &&
                (curve.evaluate(curve.t1) - solid.vertices[edge.v1]).length <
                    1e-6;
          }),
        );
        if (consistentTrims) {
          // Convex parameter polygons also need curvature-controlled chords:
          // a fan crosses the cavity when a boolean splits circular rims.
          _axialTrimmed(solid, face, fi, tess, addTri);
          return;
        }
        // Legacy blend patches may store endpoint trims that differ from
        // their support curves. Keep their vertex-based boundary here; the
        // general parameter trimmer requires consistent analytic endpoints.
        var cu = 0.0, cv = 0.0;
        for (final uv in uvs) {
          cu += uv.x;
          cv += uv.y;
        }
        cu /= uvs.length;
        cv /= uvs.length;
        final center = cyl.evaluate(cu, cv);
        for (var i = 0; i < outer.length; i++) {
          final a = outer[i], b = outer[(i + 1) % outer.length];
          var n = (a - center).cross(b - center);
          if (face.reversed) n = -n;
          if (n.dot(cyl.normalAt(cu, cv)) < 0) {
            addTri(center, b, a, fi);
          } else {
            addTri(center, a, b, fi);
          }
        }
        return;
      }
    }
    final cover = _smallestUCover(uSamples);
    final uMin = cover.$1;
    final uMax = cover.$2;
    final spanU = (uMax - uMin).abs();
    final fullU =
        closedCircle || uSamples.isEmpty || spanU > math.pi * 2 - 1e-3;
    final n = fullU
        ? tess.cylinderSegments
        : math.max(4, tess.filletArcSegments(cyl.radius, spanU));
    for (var i = 0; i < n; i++) {
      final u0 = fullU ? i / n * math.pi * 2 : uMin + (uMax - uMin) * (i / n);
      final u1 = fullU
          ? (i + 1) / n * math.pi * 2
          : uMin + (uMax - uMin) * ((i + 1) / n);
      final a = cyl.evaluate(u0, z0);
      final b = cyl.evaluate(u1, z0);
      final c = cyl.evaluate(u1, z1);
      final d = cyl.evaluate(u0, z1);
      if (face.reversed) {
        addTri(a, d, c, fi);
        addTri(a, c, b, fi);
      } else {
        addTri(a, b, c, fi);
        addTri(a, c, d, fi);
      }
    }
  }

  static List<Vector3> axialLoopParameters(
    BrepSurface s,
    List<Vector3> points,
  ) {
    if (s is BrepCone) {
      final apex = points.indexWhere(
        (p) => s.radiusAt(s.parameters(p)!.$2).abs() < 1e-8,
      );
      if (apex >= 0) {
        // Start just after the singularity. The regular path establishes the
        // angular branch before its two apex coordinates close the chart.
        points = [...points.skip(apex + 1), ...points.take(apex + 1)];
      }
    }
    final out = <Vector3>[];
    for (var index = 0; index < points.length; index++) {
      final p = points[index];
      final pair = s.parameters(p)!;
      var u = pair.$1;
      if (s is BrepCone &&
          s.radiusAt(pair.$2).abs() < 1e-8 &&
          points.length > 2) {
        // The apex is one 3D vertex but an entire angular edge in a cone's
        // parameter chart. Preserve both incoming and outgoing angles.
        var before = s
            .parameters(points[(index - 1 + points.length) % points.length])!
            .$1;
        var after = s.parameters(points[(index + 1) % points.length])!.$1;
        if (out.isNotEmpty)
          before +=
              ((out.last.x - before) / (2 * math.pi)).round() * 2 * math.pi;
        while (after - before > math.pi + 1e-10) after -= 2 * math.pi;
        while (after - before < -math.pi - 1e-10) after += 2 * math.pi;
        out.add(Vector3(before, pair.$2, 0));
        out.add(Vector3(after, pair.$2, 0));
        continue;
      }
      if ((s is BrepCylinder || s is BrepCone) && out.isNotEmpty) {
        while (u - out.last.x > math.pi) {
          u -= 2 * math.pi;
        }
        while (u - out.last.x < -math.pi) {
          u += 2 * math.pi;
        }
      }
      out.add(Vector3(u, pair.$2, 0));
    }
    return out;
  }

  static void _axialTrimmed(
    BrepSolid solid,
    BrepFace face,
    int fi,
    MeshTessellation tess,
    void Function(Vector3, Vector3, Vector3, int) addTri,
  ) {
    final surface = face.surface;
    final loops = <List<Vector3>>[], sampledLoops = <List<Vector3>>[];
    final boundaryWorld = <String, Vector3>{};
    String worldKey(double u, double v) =>
        '${(math.cos(u) * 1e8).round()}:${(math.sin(u) * 1e8).round()}:${(v * 1e8).round()}';
    Vector3 world(Vector3 uv) =>
        boundaryWorld[worldKey(uv.x, uv.y)] ?? surface.evaluate(uv.x, uv.y);
    double? anchor;
    for (final loop in face.loops) {
      final points = sampleLoop(solid, loop, tess.cylinderSegments, tess: tess);
      for (final p in points) {
        final q = surface.parameters(p)!;
        boundaryWorld[worldKey(q.$1, q.$2)] = p;
      }
      final uv = axialLoopParameters(surface, points);
      if (uv.isEmpty) continue;
      final center = uv.fold(0.0, (double sum, p) => sum + p.x) / uv.length;
      if (anchor == null) {
        anchor = center;
      } else {
        final shift = ((anchor - center) / (2 * math.pi)).round() * 2 * math.pi;
        for (final p in uv) {
          p.x += shift;
        }
      }
      loops.add(uv);
      sampledLoops.add(points);
    }
    if (loops.isEmpty) return;
    // Two closed rims bound a cylindrical band, not a planar outer ring and
    // hole in UV. They may be split across several elliptic/circular edges.
    if (surface is BrepCylinder &&
        loops.length == 2 &&
        loops.every((ring) {
          final span = ring.last.x - ring.first.x;
          if (span.abs() < math.pi) return false;
          for (var i = 1; i < ring.length; i++) {
            if ((ring[i].x - ring[i - 1].x) * span < -1e-8) return false;
          }
          return true;
        })) {
      final rows = [
        for (var k = 0; k < 2; k++)
          [
            for (var i = 0; i < loops[k].length; i++)
              (loops[k][i].x % (2 * math.pi), sampledLoops[k][i]),
          ]..sort((a, b) => a.$1.compareTo(b.$1)),
      ];
      final a = rows[0], b = rows[1];
      void triangle(Vector3 p, Vector3 q, Vector3 r) {
        final uv = surface.parameters((p + q + r) / 3)!;
        final normal =
            surface.normalAt(uv.$1, uv.$2) * (face.reversed ? -1.0 : 1.0);
        if ((q - p).cross(r - p).dot(normal) < 0)
          addTri(p, r, q, fi);
        else
          addTri(p, q, r, fi);
      }

      var i = 0, j = 0;
      while (i < a.length || j < b.length) {
        final ua = i < a.length
            ? a[(i + 1) % a.length].$1 + (i + 1 == a.length ? 2 * math.pi : 0)
            : double.infinity;
        final ub = j < b.length
            ? b[(j + 1) % b.length].$1 + (j + 1 == b.length ? 2 * math.pi : 0)
            : double.infinity;
        if (ua <= ub) {
          triangle(
            a[i % a.length].$2,
            a[(i + 1) % a.length].$2,
            b[j % b.length].$2,
          );
          i++;
        } else {
          triangle(
            a[i % a.length].$2,
            b[(j + 1) % b.length].$2,
            b[j % b.length].$2,
          );
          j++;
        }
      }
      return;
    }

    final plane = BrepPlane(
      origin: Vector3.zero(),
      normal: Vector3(0, 0, 1),
      xDir: Vector3(1, 0, 0),
    );
    final points = <Vector3>[];
    final ids = <String, int>{};
    int vertex(Vector3 p) {
      final key = '${(p.x * 1e9).round()}:${(p.y * 1e9).round()}';
      return ids.putIfAbsent(key, () {
        points.add(p);
        return points.length - 1;
      });
    }

    var triangles = <List<int>>[];
    void capture(Vector3 a, Vector3 b, Vector3 c, int _) {
      triangles.add([vertex(a), vertex(b), vertex(c)]);
    }

    _planeWithHoles(
      loops.first,
      loops.skip(1).toList(),
      plane,
      false,
      fi,
      capture,
      preferOuterBridges: true,
    );
    final initial = [
      for (final t in triangles) [for (final i in t) points[i]],
    ];
    _shortenParameterDiagonals(initial);
    final collapsedVertices = <int>{};
    triangles =
        [
          for (final t in initial) [for (final p in t) vertex(p)],
        ]..removeWhere((t) {
          if (t[0] != t[1] && t[1] != t[2] && t[2] != t[0]) return false;
          collapsedVertices.addAll(t);
          return true;
        });
    // STEP rounding can create slivers whose UV vertices weld together.
    // They have no area after welding and cannot be refined as triangles.
    // Refine interior chords in parameter space without changing the shared
    // boundary samples used by neighbouring faces.
    final maxU = 2 * math.pi / tess.cylinderSegments;
    for (
      var round = 0;
      round < math.max(1, 2 * tess.maxSurfaceRefinement - 4);
      round++
    ) {
      final uses = <String, int>{};
      String key(int a, int b) => a < b ? '$a:$b' : '$b:$a';
      for (final t in triangles) {
        for (var i = 0; i < 3; i++) {
          final k = key(t[i], t[(i + 1) % 3]);
          uses[k] = (uses[k] ?? 0) + 1;
        }
      }
      final split = <String, int>{};
      for (final entry in uses.entries) {
        if (entry.value != 2) continue;
        final ends = entry.key.split(':').map(int.parse).toList();
        final a = points[ends[0]], b = points[ends[1]];
        final joinsBoundary =
            boundaryWorld.containsKey(worldKey(a.x, a.y)) &&
            boundaryWorld.containsKey(worldKey(b.x, b.y));
        if (!joinsBoundary && (a.x - b.x).abs() <= maxU * 1.01) continue;
        if (surface is BrepCone) {
          final error =
              (world((a + b) * .5) - (world(a) + world(b)) * .5).length;
          final tolerance =
              math.max(surface.radius.abs(), surface.topRadius.abs()) *
              (1 - math.cos(maxU / 2));
          if (error <= tolerance) continue;
        }
        final mid = vertex((a + b) * .5);
        // A sub-resolution midpoint can weld to an endpoint. Splitting it
        // would reintroduce a collapsed triangle on the next refinement pass.
        if (mid != ends[0] && mid != ends[1]) split[entry.key] = mid;
      }
      if (split.isEmpty) break;
      final next = <List<int>>[];
      for (final t in triangles) {
        var pieces = [t];
        for (var i = 0; i < 3; i++) {
          final a = t[i], b = t[(i + 1) % 3], mid = split[key(a, b)];
          if (mid == null) continue;
          final updated = <List<int>>[];
          for (final piece in pieces) {
            if (!piece.contains(a) || !piece.contains(b)) {
              updated.add(piece);
              continue;
            }
            final third = piece.firstWhere((v) => v != a && v != b);
            // A collinear sliver can also weld its third vertex to the
            // midpoint; neither of the resulting pieces has any UV area.
            if (third == mid) {
              collapsedVertices.addAll(piece);
              continue;
            }
            updated.add([a, mid, third]);
            updated.add([mid, b, third]);
          }
          pieces = updated;
        }
        next.addAll(pieces);
      }
      triangles = next;
    }
    var worldTriangles = [
      for (final t in triangles) [for (final i in t) world(points[i])],
    ];
    if (surface is BrepCone || collapsedVertices.isNotEmpty) {
      // Cone apex coordinates and welded UV slivers can collapse in 3D.
      // Preserve their edge subdivisions in the surviving neighbours so
      // removing them does not leave T-junctions along the surface.
      final collapsed = [for (final i in collapsedVertices) world(points[i])];
      worldTriangles = worldTriangles.where((t) {
        final ab = t[1] - t[0], ac = t[2] - t[0], bc = t[2] - t[1];
        // On cylinders, compare altitude against the 1e-8 mm point-on-edge
        // tolerance below; absolute area depends on the length of the sliver.
        final threshold = surface is BrepCone
            ? 1e-24
            : 1e-16 * math.max(ab.length2, math.max(ac.length2, bc.length2));
        if (ab.cross(ac).length2 >= threshold) return true;
        collapsed.addAll(t);
        return false;
      }).toList();
      for (final p in collapsed) {
        final next = <List<Vector3>>[];
        for (final t in worldTriangles) {
          var split = false;
          for (var k = 0; k < 3; k++) {
            final a = t[k], b = t[(k + 1) % 3], c = t[(k + 2) % 3], d = b - a;
            if (d.length2 < 1e-24) continue;
            final fraction = (p - a).dot(d) / d.length2;
            if (fraction <= 1e-9 ||
                fraction >= 1 - 1e-9 ||
                (a + d * fraction - p).length > 1e-8)
              continue;
            next.addAll([
              [a, p, c],
              [p, b, c],
            ]);
            split = true;
            break;
          }
          if (!split) next.add(t);
        }
        worldTriangles = next;
      }
    }
    for (final t in worldTriangles) {
      final a = t[0], b = t[1], c = t[2];
      final uv = surface.parameters((a + b + c) / 3)!;
      var normal = surface.normalAt(uv.$1, uv.$2);
      if (face.reversed) normal = -normal;
      if ((b - a).cross(c - a).dot(normal) < 0) {
        addTri(a, c, b, fi);
      } else {
        addTri(a, b, c, fi);
      }
    }
  }

  /// Loft a full cylinder wall from the same closed-circle samples the
  /// planar caps use, so the rims weld even when cylinder U and circle t
  /// frames differ after pack / transform.
  static bool _cylinderFromClosedRims(
    BrepSolid solid,
    BrepFace face,
    BrepCylinder cyl,
    MeshTessellation tess,
    int fi,
    void Function(Vector3, Vector3, Vector3, int) addTri,
  ) {
    final rims = <List<Vector3>>[];
    for (final loop in face.loops) {
      for (final c in loop.coedges) {
        final e = solid.edges[c.edge];
        final curve = e.curve;
        final List<Vector3> pts;
        if (curve is BrepCircle && curve.isClosed) {
          pts = _sampleClosedCircle(c, curve, tess.cylinderSegments);
        } else if (curve is BrepEllipse && curve.isClosed) {
          pts = _sampleClosedEllipse(c, curve, tess.cylinderSegments);
        } else {
          continue;
        }
        if (pts.length >= 3) rims.add(pts);
      }
    }
    if (rims.length < 2) return false;
    final a = rims[0];
    var b = rims[1];
    if (a.length != b.length || a.length < 3) return false;

    var best = 0;
    var bestScore = -1.0;
    for (var i = 0; i < b.length; i++) {
      final d = b[i] - a[0];
      final score = d.dot(cyl.axis).abs() / math.max(d.length, 1e-12);
      if (score > bestScore) {
        bestScore = score;
        best = i;
      }
    }
    if (best != 0) {
      b = [...b.sublist(best), ...b.sublist(0, best)];
    }
    final dFwd = (b[1] - a[1]).length2 + (b[0] - a[0]).length2;
    final dRev = (b[b.length - 1] - a[1]).length2 + (b[0] - a[0]).length2;
    if (dRev < dFwd) {
      b = [b.first, ...b.sublist(1).reversed];
    }

    for (var i = 0; i < a.length; i++) {
      final i1 = (i + 1) % a.length;
      final p00 = a[i];
      final p10 = a[i1];
      final p01 = b[i];
      final p11 = b[i1];
      final mid = (p00 + p10 + p01 + p11) * 0.25;
      final uv = cyl.parameters(mid) ?? (0.0, 0.0);
      var n = (p10 - p00).cross(p01 - p00);
      if (n.length2 < 1e-20) {
        n = (p11 - p00).cross(p01 - p00);
      }
      if (face.reversed) n = -n;
      if (n.dot(cyl.normalAt(uv.$1, uv.$2)) < 0) {
        addTri(p00, p01, p11, fi);
        addTri(p00, p11, p10, fi);
      } else {
        addTri(p00, p10, p11, fi);
        addTri(p00, p11, p01, fi);
      }
    }
    return true;
  }

  /// Mesh a 4-edge cylinder patch from its two circular edges so the rims
  /// share vertices with planar faces that [sampleLoop] the same curves.
  static bool _cylinderPatchFromArcs(
    BrepSolid solid,
    BrepFace face,
    BrepCylinder cyl,
    MeshTessellation tess,
    int fi,
    void Function(Vector3, Vector3, Vector3, int) addTri,
  ) {
    final arcs = <List<Vector3>>[];
    for (final c in face.loops.first.coedges) {
      final e = solid.edges[c.edge];
      final curve = e.curve;
      if (curve is! BrepCircle || curve.isClosed) continue;
      final n = tess.filletArcSegments(
        curve.radius,
        (curve.t1 - curve.t0).abs(),
      );
      final start = solid.vertices[c.reversed ? e.v1 : e.v0];
      final end = solid.vertices[c.reversed ? e.v0 : e.v1];
      final a0 = curve.evaluate(curve.t0);
      final a1 = curve.evaluate(curve.t1);
      final forward = (start - a0).length <= (start - a1).length;
      final pts = <Vector3>[];
      for (var i = 0; i <= n; i++) {
        if (i == 0) {
          pts.add(start.clone());
          continue;
        }
        if (i == n) {
          pts.add(end.clone());
          continue;
        }
        final frac = i / n;
        final t = forward
            ? curve.t0 + (curve.t1 - curve.t0) * frac
            : curve.t1 - (curve.t1 - curve.t0) * frac;
        pts.add(curve.evaluate(t));
      }
      arcs.add(pts);
    }
    if (arcs.length != 2 || arcs[0].length < 2) return false;
    if (arcs[0].length != arcs[1].length) return false;
    final a = arcs[0];
    var b = arcs[1];
    if ((a.first - b.first).length > (a.first - b.last).length) {
      b = b.reversed.toList();
    }
    for (var i = 0; i < a.length - 1; i++) {
      final p00 = a[i];
      final p10 = a[i + 1];
      final p01 = b[i];
      final p11 = b[i + 1];
      final mid = (p00 + p10 + p01 + p11) * 0.25;
      final uv = cyl.parameters(mid) ?? (0.0, 0.0);
      var n = (p10 - p00).cross(p01 - p00);
      if (n.length2 < 1e-20) {
        n = (p11 - p00).cross(p01 - p00);
      }
      if (face.reversed) n = -n;
      if (n.dot(cyl.normalAt(uv.$1, uv.$2)) < 0) {
        addTri(p00, p01, p11, fi);
        addTri(p00, p11, p10, fi);
      } else {
        addTri(p00, p10, p11, fi);
        addTri(p00, p11, p01, fi);
      }
    }
    return true;
  }

  /// Band between two circular rims with extra loops punched as holes, so
  /// the hole vertices match planar faces that sample the same curves.
  static bool _sphereBandWithHoles(
    BrepSolid solid,
    BrepSphere s,
    BrepFace face,
    int fi,
    MeshTessellation tess,
    void Function(Vector3, Vector3, Vector3, int) addTri,
  ) {
    final circleIdx = <int>[];
    for (var i = 0; i < face.loops.length; i++) {
      final loop = face.loops[i];
      if (loop.coedges.length != 1) continue;
      final ei = loop.coedges.first.edge;
      if (ei < 0 || ei >= solid.edges.length) continue;
      final curve = solid.edges[ei].curve;
      if (curve is BrepCircle && curve.isClosed) circleIdx.add(i);
    }
    if (circleIdx.length < 2) return false;
    var i0 = circleIdx[0];
    var i1 = circleIdx[1];
    if (circleIdx.length > 2) {
      var best = -1.0;
      for (var a = 0; a < circleIdx.length; a++) {
        for (var b = a + 1; b < circleIdx.length; b++) {
          final ca =
              solid.edges[face.loops[circleIdx[a]].coedges.first.edge].curve
                  as BrepCircle;
          final cb =
              solid.edges[face.loops[circleIdx[b]].coedges.first.edge].curve
                  as BrepCircle;
          final sep = (ca.center - cb.center).length;
          final aligned = ca.normal.dot(cb.normal).abs();
          final score = sep * (0.5 + 0.5 * aligned);
          if (score > best) {
            best = score;
            i0 = circleIdx[a];
            i1 = circleIdx[b];
          }
        }
      }
    }
    final extra = <int>[
      for (var i = 0; i < face.loops.length; i++)
        if (i != i0 && i != i1) i,
    ];
    if (extra.isEmpty) return false;
    var a = sampleLoop(
      solid,
      face.loops[i0],
      tess.sphereLonSegments,
      tess: tess,
    );
    var b = sampleLoop(
      solid,
      face.loops[i1],
      tess.sphereLonSegments,
      tess: tess,
    );
    if (a.length < 3 || b.length < 3) return false;
    var best = 0;
    var bestD = double.infinity;
    for (var i = 0; i < b.length; i++) {
      final d = (b[i] - a.first).length2;
      if (d < bestD) {
        bestD = d;
        best = i;
      }
    }
    b = [for (var i = 0; i < b.length; i++) b[(best + i) % b.length]];
    if (b.length != a.length) {
      b = [
        for (var i = 0; i < a.length; i++)
          b[((i / a.length) * b.length).floor() % b.length],
      ];
    }
    if (a.length >= 3) {
      final dFwd = (b[1] - a[1]).length2;
      final dRev = (b[1] - a[a.length - 1]).length2;
      if (dRev < dFwd) {
        b = [b.first, ...b.sublist(1).reversed];
      }
    }

    (double, double)? uvOf(Vector3 p, double? prevU) {
      final uv = s.parameters(p);
      if (uv == null) return null;
      var u = uv.$1;
      if (u < 0) u += math.pi * 2;
      if (prevU != null) {
        while (u - prevU > math.pi) {
          u -= math.pi * 2;
        }
        while (prevU - u > math.pi) {
          u += math.pi * 2;
        }
      }
      return (u, uv.$2);
    }

    List<Vector3>? ringUv(List<Vector3> pts) {
      final out = <Vector3>[];
      double? prev;
      for (final p in pts) {
        final uv = uvOf(p, prev);
        if (uv == null) return null;
        prev = uv.$1;
        out.add(Vector3(uv.$1, uv.$2, 0));
      }
      return out.length >= 3 ? out : null;
    }

    final aUv = ringUv(a);
    if (aUv == null) return false;
    final pairedB = <Vector3>[];
    for (var i = 0; i < aUv.length; i++) {
      final bu = uvOf(b[i % b.length], aUv[i].x);
      if (bu == null) return false;
      pairedB.add(Vector3(aUv[i].x, bu.$2, 0));
    }
    final outer = [...aUv, ...pairedB.reversed];
    var minU = outer.first.x, maxU = outer.first.x;
    var minV = outer.first.y, maxV = outer.first.y;
    for (final p in outer) {
      minU = math.min(minU, p.x);
      maxU = math.max(maxU, p.x);
      minV = math.min(minV, p.y);
      maxV = math.max(maxV, p.y);
    }
    final samples = <Vector3>[...a, ...b];
    final holes = <List<Vector3>>[];
    for (final hi in extra) {
      final pts = sampleLoop(
        solid,
        face.loops[hi],
        tess.sphereLonSegments,
        tess: tess,
      );
      samples.addAll(pts);
      final uv = ringUv(pts);
      if (uv == null) continue;
      var meanU = 0.0;
      for (final p in uv) {
        meanU += p.x;
      }
      meanU /= uv.length;
      var shift = 0.0;
      final mid = 0.5 * (minU + maxU);
      while (meanU + shift - mid > math.pi) {
        shift -= math.pi * 2;
      }
      while (mid - (meanU + shift) > math.pi) {
        shift += math.pi * 2;
      }
      if (shift != 0) {
        for (final p in uv) {
          p.x += shift;
        }
      }
      var holeMinV = uv.first.y, holeMaxV = uv.first.y;
      for (final p in uv) {
        holeMinV = math.min(holeMinV, p.y);
        holeMaxV = math.max(holeMaxV, p.y);
      }
      if (holeMaxV < minV - 0.05 || holeMinV > maxV + 0.05) continue;
      holes.add(uv);
    }
    if (holes.isEmpty) return false;
    final dummy = BrepPlane(
      origin: Vector3.zero(),
      normal: Vector3(0, 0, 1),
      xDir: Vector3(1, 0, 0),
    );
    Vector3 snap(Vector3 uv) {
      final p = s.evaluate(uv.x, uv.y);
      var best = p;
      var bestD = 4e-4;
      for (final q in samples) {
        final d = (q - p).length2;
        if (d < bestD) {
          bestD = d;
          best = q;
        }
      }
      return best;
    }

    void mapped(Vector3 A, Vector3 B, Vector3 C, int id) {
      final pa = snap(A);
      final pb = snap(B);
      final pc = snap(C);
      var n = (pb - pa).cross(pc - pa);
      if (face.reversed) n = -n;
      final mid = (pa + pb + pc) / 3;
      if (n.dot(mid - s.center) < 0) {
        addTri(pa, pc, pb, id);
      } else {
        addTri(pa, pb, pc, id);
      }
    }

    var pts = List<Vector3>.from(outer);
    for (final hole in holes) {
      var bestOi = 0;
      var bestHi = 0;
      var bestD = double.infinity;
      for (var oi = 0; oi < pts.length; oi++) {
        for (var hi = 0; hi < hole.length; hi++) {
          final d = (pts[oi] - hole[hi]).length2;
          if (d < bestD) {
            bestD = d;
            bestOi = oi;
            bestHi = hi;
          }
        }
      }
      pts = [
        ...pts.sublist(0, bestOi + 1),
        for (var k = 0; k <= hole.length; k++) hole[(bestHi + k) % hole.length],
        pts[bestOi],
        ...pts.sublist(bestOi + 1),
      ];
    }
    _earClipPlane(pts, dummy, false, fi, mapped, centroidFallback: false);
    return true;
  }

  /// Trimmed fillet sphere: barycentric lattice on a 3-arc spherical
  /// triangle, otherwise concentric rings. A single fan from one interior
  /// point puts a pole on the patch and reads as a pointed corner.
  static bool _spherePatch(
    BrepSolid solid,
    BrepSphere s,
    BrepFace face,
    int fi,
    MeshTessellation tess,
    void Function(Vector3, Vector3, Vector3, int) addTri,
  ) {
    final loop = face.loops.first;
    if (face.loops.length > 2 &&
        _sphereBandWithHoles(solid, s, face, fi, tess, addTri)) {
      return true;
    }
    if (face.loops.length >= 2) {
      final a = sampleLoop(
        solid,
        face.loops[0],
        tess.sphereLonSegments,
        tess: tess,
      );
      var b = sampleLoop(
        solid,
        face.loops[1],
        tess.sphereLonSegments,
        tess: tess,
      );
      if (a.length >= 3 && b.length >= 3) {
        final n = a.length;
        // Rotate [b] so b[0] is nearest a[0].
        var best = 0;
        var bestD = double.infinity;
        for (var i = 0; i < b.length; i++) {
          final d = (b[i] - a.first).length2;
          if (d < bestD) {
            bestD = d;
            best = i;
          }
        }
        b = [for (var i = 0; i < b.length; i++) b[(best + i) % b.length]];
        if (b.length != n) {
          b = [
            for (var i = 0; i < n; i++)
              b[((i / n) * b.length).floor() % b.length],
          ];
        }
        if (n >= 3) {
          final dFwd = (b[1] - a[1]).length2;
          final dRev = (b[1] - a[n - 1]).length2;
          if (dRev < dFwd) {
            b = [b.first, ...b.sublist(1).reversed];
          }
        }
        final rings = math.max(3, tess.sphereLatSegments);
        Vector3 lerp(Vector3 p, Vector3 q, double t) {
          // Interpolate on the sphere itself. Latitude/longitude interpolation
          // folds tilted bands where a rim crosses the parameter-space pole.
          if (t == 0) return p;
          if (t == 1) return q;
          return s.center +
              ((p - s.center) * (1 - t) + (q - s.center) * t).normalized() *
                  s.radius;
        }

        void emit(Vector3 a, Vector3 b, Vector3 d) {
          var nrm = (b - a).cross(d - a);
          if (face.reversed) nrm = -nrm;
          final mid = (a + b + d) / 3;
          if (nrm.dot(mid - s.center) < 0) {
            addTri(a, d, b, fi);
          } else {
            addTri(a, b, d, fi);
          }
        }

        for (var r = 0; r < rings; r++) {
          final t0 = r / rings;
          final t1 = (r + 1) / rings;
          for (var i = 0; i < n; i++) {
            final j = (i + 1) % n;
            final a0 = lerp(a[i], b[i], t0);
            final a1 = lerp(a[j], b[j], t0);
            final b0 = lerp(a[i], b[i], t1);
            final b1 = lerp(a[j], b[j], t1);
            emit(a0, a1, b1);
            emit(a0, b1, b0);
          }
        }
        return true;
      }
    }
    final outer = sampleLoop(solid, loop, tess.sphereLonSegments, tess: tess);
    if (outer.length < 3) return false;
    var c = Vector3.zero();
    for (final p in outer) {
      c += p;
    }
    c = s.closestPoint(c / outer.length.toDouble());
    if (loop.coedges.length == 1) {
      final e = solid.edges[loop.coedges.first.edge];
      final curve = e.curve;
      if (curve is BrepCircle && curve.isClosed) {
        c = s.circularCapPole(
          curve,
          coedgeReversed: loop.coedges.first.reversed,
        );
      } else {
        final small = s.closestPoint(c);
        final large = s.center * 2.0 - small;
        var tang = curve.tangent(curve.t0);
        final p0 = curve.evaluate(curve.t0);
        c = s.loopFillPole(
          p0,
          tang,
          coedgeReversed: loop.coedges.first.reversed,
          faceReversed: face.reversed,
          small: small,
          large: large,
        );
      }
    } else if (loop.coedges.isNotEmpty) {
      final e = solid.edges[loop.coedges.first.edge];
      final curve = e.curve;
      final small = s.closestPoint(c);
      final large = s.center * 2.0 - small;
      c = s.loopFillPole(
        curve.evaluate(curve.t0),
        curve.tangent(curve.t0),
        coedgeReversed: loop.coedges.first.reversed,
        faceReversed: face.reversed,
        small: small,
        large: large,
      );
    }
    final ang = math.acos(
      ((c - s.center).normalized().dot(
        (outer.first - s.center).normalized(),
      )).clamp(-1.0, 1.0),
    );
    final rings = math.max(
      3,
      math.max(
        tess.filletArcSegments(s.radius, math.pi / 2) ~/ 2,
        (tess.sphereLatSegments * ang / math.pi).ceil(),
      ),
    );
    final pts = <List<Vector3>>[
      [c],
      for (var r = 1; r <= rings; r++)
        [for (final p in outer) s.closestPoint(c + (p - c) * (r / rings))],
    ];
    void emit(Vector3 a, Vector3 b, Vector3 d) {
      var n = (b - a).cross(d - a);
      if (face.reversed) n = -n;
      final mid = (a + b + d) / 3;
      if (n.dot(mid - s.center) < 0) {
        addTri(a, d, b, fi);
      } else {
        addTri(a, b, d, fi);
      }
    }

    final inner = pts[1];
    for (var i = 0; i < inner.length; i++) {
      emit(c, inner[i], inner[(i + 1) % inner.length]);
    }
    for (var r = 1; r < rings; r++) {
      final a = pts[r];
      final b = pts[r + 1];
      for (var i = 0; i < a.length; i++) {
        final j = (i + 1) % a.length;
        emit(a[i], b[i], b[j]);
        emit(a[i], b[j], a[j]);
      }
    }
    return true;
  }

  /// A finite, outward-oriented gnomonic chart for a spherical patch contained
  /// in one hemisphere. Circular caps/bands keep their existing analytic path.
  static BrepPlane? sphereFaceFrame(BrepSolid solid, BrepFace face) {
    final s = face.surface;
    if (s is! BrepSphere ||
        face.loops.isEmpty ||
        face.loops.every(
          (l) =>
              l.coedges.length == 1 &&
              solid.edges[l.coedges.single.edge].curve is BrepCircle &&
              (solid.edges[l.coedges.single.edge].curve as BrepCircle).isClosed,
        ))
      return null;
    final points = sampleLoop(solid, face.loops.first, 24);
    if (points.isEmpty) return null;
    var n = Vector3.zero();
    for (final p in points) n += (p - s.center).normalized();
    if (n.length < 1e-8) return null;
    n.normalize();
    if (points.any((p) => (p - s.center).dot(n) < s.radius * 1e-5)) return null;
    final frame = BrepPlane(origin: s.center, normal: n);
    final uv = [for (final p in points) sphereProject(frame, p)];
    var area = 0.0;
    for (var i = 0; i < uv.length; i++) {
      final a = uv[i], b = uv[(i + 1) % uv.length];
      area += a.x * b.y - a.y * b.x;
    }
    // A loop contained in a hemisphere can also bound its large complement.
    // Only use this finite chart when the oriented outer wire fills its inside.
    if ((face.reversed ? -area : area) <= 1e-12) return null;
    return frame;
  }

  static Vector2 sphereProject(BrepPlane frame, Vector3 p) {
    final d = p - frame.origin, den = d.dot(frame.normal);
    return Vector2(d.dot(frame.xDir) / den, d.dot(frame.yDir) / den);
  }

  static void _sphereChartPatch(
    BrepSolid solid,
    BrepFace face,
    BrepPlane frame,
    int fi,
    MeshTessellation tess,
    void Function(Vector3, Vector3, Vector3, int) emit,
  ) {
    final sphere = face.surface as BrepSphere;
    final points = <Vector3>[], world = <Vector3>[], rings = <List<int>>[];
    for (final loop in face.loops) {
      final ring = <int>[];
      for (final p in sampleLoop(
        solid,
        loop,
        tess.cylinderSegments,
        tess: tess,
      )) {
        final uv = sphereProject(frame, p);
        ring.add(points.length);
        points.add(Vector3(uv.x, uv.y, 0));
        world.add(p);
      }
      rings.add(ring);
    }
    var triangles = MeshPolygon.triangulateRings(
      points,
      rings,
      Vector3(0, 0, 1),
      keepCollinear: true,
    );
    if (triangles.isEmpty) {
      // Great-circle arcs become straight lines in this chart. STEP rounding
      // can leave their samples a few trillionths off that line, stranding
      // ear clipping on a nearly flat remainder. Retry with those subdivisions
      // removed, then restore every original boundary sample by edge identity.
      // The world-space samples and exact BREP geometry remain unchanged.
      final simplified = [for (final ring in rings) List<int>.from(ring)];
      final removed = <(int, int, int)>[];
      for (final ring in simplified) {
        var changed = true;
        while (changed && ring.length > 3) {
          changed = false;
          for (var i = 0; i < ring.length; i++) {
            final a = ring[(i + ring.length - 1) % ring.length];
            final b = ring[i], c = ring[(i + 1) % ring.length];
            final d = points[c] - points[a];
            if (d.length2 < 1e-24) continue;
            final t = (points[b] - points[a]).dot(d) / d.length2;
            if (t <= 0 ||
                t >= 1 ||
                (points[a] + d * t - points[b]).length2 > 1e-20) {
              continue;
            }
            removed.add((a, b, c));
            ring.removeAt(i);
            changed = true;
            break;
          }
        }
      }
      triangles = MeshPolygon.triangulateRings(
        points,
        simplified,
        Vector3(0, 0, 1),
        keepCollinear: true,
      );
      for (final (a, b, c) in removed.reversed) {
        final index = triangles.indexWhere(
          (t) => t.contains(a) && t.contains(c),
        );
        if (index < 0) {
          triangles = [];
          break;
        }
        final triangle = triangles.removeAt(index);
        final opposite = triangle.firstWhere((v) => v != a && v != c);
        triangles.addAll([
          [a, b, opposite],
          [b, c, opposite],
        ]);
      }
    }
    final boundary = <(int, int)>{};
    (int, int) key(int a, int b) => a < b ? (a, b) : (b, a);
    for (final ring in rings) {
      for (var i = 0; i < ring.length; i++)
        boundary.add(key(ring[i], ring[(i + 1) % ring.length]));
    }
    // Shared BREP samples remain fixed; split each interior chord on both
    // adjacent triangles, including chords between two boundary vertices.
    final originalCount = points.length;
    final chordError =
        sphere.radius * (1 - math.cos(math.pi / tess.cylinderSegments));
    Vector3 evaluate(Vector3 uv) =>
        sphere.center +
        (frame.normal + frame.xDir * uv.x + frame.yDir * uv.y).normalized() *
            sphere.radius;
    _improveChartTriangles(points, triangles);
    for (var pass = 0; pass < 12; pass++) {
      final splits = <(int, int), int>{};
      for (final t in triangles) {
        for (var k = 0; k < 3; k++) {
          final a = t[k], b = t[(k + 1) % 3], edge = key(a, b);
          if (boundary.contains(edge) || splits.containsKey(edge)) continue;
          final uv = (points[a] + points[b]) * .5, p = evaluate(uv);
          if (!(a < originalCount && b < originalCount) &&
              (p - (world[a] + world[b]) * .5).length <= chordError)
            continue;
          splits[edge] = points.length;
          points.add(uv);
          world.add(p);
        }
      }
      if (splits.isEmpty) break;
      final next = <List<int>>[];
      for (final t in triangles) {
        var pieces = [t];
        for (var k = 0; k < 3; k++) {
          final a = t[k], b = t[(k + 1) % 3], mid = splits[key(a, b)];
          if (mid == null) continue;
          final updated = <List<int>>[];
          for (final piece in pieces) {
            if (!piece.contains(a) || !piece.contains(b)) {
              updated.add(piece);
              continue;
            }
            final c = piece.firstWhere((v) => v != a && v != b);
            updated.addAll([
              [a, mid, c],
              [mid, b, c],
            ]);
          }
          pieces = updated;
        }
        next.addAll(pieces);
      }
      triangles = next;
      _improveChartTriangles(points, triangles);
    }
    for (final t in triangles) {
      final a = world[t[0]], b = world[t[1]], c = world[t[2]];
      final outward =
          (b - a).cross(c - a).dot((a + b + c) / 3 - sphere.center) > 0;
      if (outward != face.reversed)
        emit(a, b, c, fi);
      else
        emit(a, c, b, fi);
    }
  }

  /// Constrained Delaunay flips keep the shared boundary fixed while removing
  /// the long, needle-shaped fans produced by ear clipping and subdivision.
  /// Chord error alone does not control a triangle's normal: a very thin
  /// spherical triangle can be almost perpendicular to the true surface.
  static void _improveChartTriangles(
    List<Vector3> points,
    List<List<int>> triangles,
  ) {
    double cross(int a, int b, int c) {
      final u = points[b] - points[a], v = points[c] - points[a];
      return u.x * v.y - u.y * v.x;
    }

    double cot(int a, int center, int b) {
      final u = points[a] - points[center], v = points[b] - points[center];
      return u.dot(v) / (u.x * v.y - u.y * v.x).abs();
    }

    (int, int) edge(int a, int b) => a < b ? (a, b) : (b, a);
    final uses = <(int, int), Set<int>>{};
    final pending = <(int, int)>[], queued = <(int, int)>{};
    void add(int i, {bool enqueue = false}) {
      final t = triangles[i];
      for (var k = 0; k < 3; k++) {
        final e = edge(t[k], t[(k + 1) % 3]);
        (uses[e] ??= {}).add(i);
        if (enqueue && queued.add(e)) pending.add(e);
      }
    }

    void remove(int i) {
      final t = triangles[i];
      for (var k = 0; k < 3; k++) {
        final e = edge(t[k], t[(k + 1) % 3]);
        uses[e]!.remove(i);
        if (uses[e]!.isEmpty) uses.remove(e);
      }
    }

    for (var i = 0; i < triangles.length; i++) add(i, enqueue: true);
    // Update only the neighbours of a flipped edge. A fixed number of global
    // sweeps left long fan diagonals on slender multi-turn STEP patches;
    // subdivision then multiplied those diagonals into millions of triangles.
    while (pending.isNotEmpty) {
      final e = pending.removeLast();
      queued.remove(e);
      final adjacent = uses[e];
      if (adjacent == null || adjacent.length != 2) continue;
      final i = adjacent.first, j = adjacent.last, a = e.$1, b = e.$2;
      final c = triangles[i].firstWhere((v) => v != a && v != b);
      final d = triangles[j].firstWhere((v) => v != a && v != b);
      if (c == d ||
          cross(a, b, c) * cross(a, b, d) >= 0 ||
          cross(c, d, a) * cross(c, d, b) >= 0)
        continue;
      if (cot(a, c, b) + cot(a, d, b) >= -1e-10) continue;
      remove(i);
      remove(j);
      triangles[i] = [c, d, a];
      triangles[j] = [d, c, b];
      add(i, enqueue: true);
      add(j, enqueue: true);
    }
  }

  static void _sphere(
    BrepSolid solid,
    BrepSphere s,
    BrepFace face,
    int fi,
    MeshTessellation tess,
    void Function(Vector3, Vector3, Vector3, int) addTri,
  ) {
    final frame = sphereFaceFrame(solid, face);
    if (frame != null) {
      _sphereChartPatch(solid, face, frame, fi, tess, addTri);
      return;
    }
    if (face.loops.isNotEmpty && face.loops.first.coedges.isNotEmpty) {
      if (_spherePatch(solid, s, face, fi, tess, addTri)) return;
    }
    final lat = tess.sphereLatSegments;
    final lon = tess.sphereLonSegments;
    final poleN = s.evaluate(0, 0);
    final poleS = s.evaluate(0, math.pi);

    // A sphere's (u, v) is left-handed about its own outward normal: v runs
    // from the +z pole downward, which makes ∂u × ∂v point *into* the sphere.
    // Every other surface in the catalog is right-handed, so this is the one
    // place a triangle wound counter-clockwise in parameter space has to be
    // flipped to end up facing out.
    void emit(Vector3 a, Vector3 b, Vector3 c) {
      if (face.reversed) {
        addTri(a, b, c, fi);
      } else {
        addTri(a, c, b, fi);
      }
    }

    for (var i = 0; i < lat; i++) {
      final v0 = i / lat * math.pi;
      final v1 = (i + 1) / lat * math.pi;
      for (var j = 0; j < lon; j++) {
        final u0 = j / lon * math.pi * 2;
        final u1 = (j + 1) / lon * math.pi * 2;
        if (i == 0) {
          // Same winding as the untrimmed quad's a→c→d triangle
          // (pole, +u, then +v-back), so emit()'s inward-param flip
          // faces the cap outward. The other order left a hole at +Z.
          emit(poleN, s.evaluate(u1, v1), s.evaluate(u0, v1));
          continue;
        }
        if (i == lat - 1) {
          emit(s.evaluate(u0, v0), s.evaluate(u1, v0), poleS);
          continue;
        }
        final a = s.evaluate(u0, v0);
        final b = s.evaluate(u1, v0);
        final c = s.evaluate(u1, v1);
        final d = s.evaluate(u0, v1);
        emit(a, b, c);
        emit(a, c, d);
      }
    }
  }

  static void _torus(
    BrepSolid solid,
    BrepTorus s,
    BrepFace face,
    int fi,
    MeshTessellation tess,
    void Function(Vector3, Vector3, Vector3, int) addTri,
  ) {
    final nu = tess.cylinderSegments;
    var u0 = 0.0, u1 = math.pi * 2;
    var v0 = 0.0, v1 = math.pi * 2;
    final domain = torusFaceDomain(solid, s, face);
    if (domain != null) {
      u0 = domain.$1;
      u1 = domain.$2;
      v0 = domain.$3;
      v1 = domain.$4;
    }
    final uSpan = (u1 - u0).abs();
    final vSpan = (v1 - v0).abs();
    // Partial toroidal rims have different radii at their two ends. Sample
    // each shared circle exactly as its adjacent planar/cylindrical face,
    // then zipper adjacent rows instead of forcing both onto one angular grid.
    if (uSpan < math.pi * 2 - 1e-6 &&
        vSpan < math.pi - 1e-7 &&
        face.loops.length == 1 &&
        face.loops.single.coedges.every((co) {
          final c = solid.edges[co.edge].curve;
          if (c is! BrepCircle) return false;
          final uv = [
            for (final t in [0.0, .25, .5, .75, 1.0])
              s.parameters(c.evaluate(c.t0 + (c.t1 - c.t0) * t))!,
          ];
          bool fixed(bool u) => uv.every(
            (p) =>
                math
                    .sin(
                      ((u ? p.$1 : p.$2) - (u ? uv.first.$1 : uv.first.$2)) *
                          .5,
                    )
                    .abs() <
                1e-7,
          );
          return fixed(true) || fixed(false);
        })) {
      final nv = tess.filletArcSegments(s.minorRadius, vSpan);
      List<Vector3> row(int j) {
        final v = v0 + (v1 - v0) * j / nv;
        final radius = s.majorRadius + s.minorRadius * math.cos(v);
        final count = tess.filletArcSegments(radius.abs(), uSpan);
        return [
          for (var i = 0; i <= count; i++)
            s.evaluate(u0 + (u1 - u0) * i / count, v),
        ];
      }

      void emit(Vector3 a, Vector3 b, Vector3 c) {
        if (face.reversed) {
          addTri(a, c, b, fi);
        } else {
          addTri(a, b, c, fi);
        }
      }

      var a = row(0);
      for (var j = 1; j <= nv; j++) {
        final b = row(j), na = a.length - 1, nb = b.length - 1;
        var ia = 0, ib = 0;
        while (ia < na || ib < nb) {
          if (ia < na && (ib == nb || (ia + 1) / na <= (ib + 1) / nb)) {
            emit(a[ia], a[ia + 1], b[ib]);
            ia++;
          } else {
            emit(a[ia], b[ib + 1], b[ib]);
            ib++;
          }
        }
        a = b;
      }
      return;
    }
    final nU = math.max(4, (nu * uSpan / (math.pi * 2)).ceil());
    final nV = vSpan > math.pi * 1.5
        ? math.max(4, (nu * vSpan / (math.pi * 2)).ceil())
        : math.max(4, tess.filletArcSegments(s.minorRadius, vSpan));
    var phase0 = 0.0, phase1 = 0.0;
    if (vSpan > math.pi * 1.99) {
      for (final loop in face.loops) {
        if (loop.coedges.length != 1) continue;
        final c = solid.edges[loop.coedges.single.edge].curve;
        if (c is! BrepCircle || !c.isClosed) continue;
        final uv = s.parameters(c.evaluate(c.t0))!;
        double distance(double u) => (math.sin((uv.$1 - u) / 2)).abs();
        if (distance(u0) < 1e-6) phase0 = uv.$2;
        if (distance(u1) < 1e-6) phase1 = uv.$2;
      }
      while (phase1 - phase0 > math.pi) {
        phase1 -= 2 * math.pi;
      }
      while (phase1 - phase0 < -math.pi) {
        phase1 += 2 * math.pi;
      }
    }
    for (var i = 0; i < nU; i++) {
      final ua = u0 + (u1 - u0) * (i / nU);
      final ub = u0 + (u1 - u0) * ((i + 1) / nU);
      for (var j = 0; j < nV; j++) {
        final va = v0 + (v1 - v0) * (j / nV);
        final vb = v0 + (v1 - v0) * ((j + 1) / nV);
        // Shared rim circles may use different angular frames at a 3D
        // bend. Interpolate the phase on the same torus, matching both rims.
        final pa = phase0 + (phase1 - phase0) * i / nU;
        final pb = phase0 + (phase1 - phase0) * (i + 1) / nU;
        final a = s.evaluate(ua, va + pa);
        final b = s.evaluate(ub, va + pb);
        final c = s.evaluate(ub, vb + pb);
        final d = s.evaluate(ua, vb + pa);
        if (face.reversed) {
          addTri(a, d, c, fi);
          addTri(a, c, b, fi);
        } else {
          addTri(a, b, c, fi);
          addTri(a, c, d, fi);
        }
      }
    }
  }

  /// UV box covering [face] loops. Empty loops → the whole torus.
  static (double, double, double, double)? torusFaceDomain(
    BrepSolid solid,
    BrepTorus s,
    BrepFace face,
  ) {
    if (face.loops.isEmpty) return null;
    final us = <double>[];
    final vs = <double>[];
    void consider(Vector3 p) {
      final uv = s.parameters(p);
      if (uv == null) return;
      us.add(uv.$1);
      vs.add(uv.$2);
    }

    for (final loop in face.loops) {
      for (final co in loop.coedges) {
        if (co.edge < 0 || co.edge >= solid.edges.length) continue;
        final curve = solid.edges[co.edge].curve;
        final n = curve is BrepCircle ? (curve.isClosed ? 24 : 12) : 4;
        for (var i = 0; i <= n; i++) {
          final t = curve.t0 + (curve.t1 - curve.t0) * (i / n);
          consider(curve.evaluate(t));
        }
      }
    }
    if (us.isEmpty) return (0, math.pi * 2, 0, math.pi / 2);
    final u = _periodicSpan(us);
    final v = _periodicSpan(vs);
    return (u.$1, u.$2, v.$1, v.$2);
  }

  /// Smallest periodic interval covering [raw], or a full turn when samples
  /// wrap all the way around.
  static (double, double) _periodicSpan(
    List<double> raw, {
    double period = math.pi * 2,
  }) {
    if (raw.isEmpty) return (0, period);
    final sorted = <double>[];
    for (final x in raw) {
      var t = x % period;
      if (t < 0) t += period;
      if (sorted.isEmpty || (t - sorted.last).abs() > 1e-6) {
        var inserted = false;
        for (var i = 0; i < sorted.length; i++) {
          if ((t - sorted[i]).abs() <= 1e-6) {
            inserted = true;
            break;
          }
          if (t < sorted[i]) {
            sorted.insert(i, t);
            inserted = true;
            break;
          }
        }
        if (!inserted) sorted.add(t);
      }
    }
    if (sorted.length == 1) {
      return (sorted.first, sorted.first + 1e-6);
    }
    var maxGap = 0.0;
    var maxI = 0;
    for (var i = 0; i < sorted.length; i++) {
      final next = i + 1 < sorted.length
          ? sorted[i + 1]
          : sorted.first + period;
      final gap = next - sorted[i];
      if (gap > maxGap) {
        maxGap = gap;
        maxI = i;
      }
    }
    // ~40° — denser than an octagon chord, looser than a 90° sector gap.
    if (maxGap < 0.7) return (0, period);
    final start = sorted[(maxI + 1) % sorted.length];
    return (start, start + (period - maxGap));
  }

  /// Ruled tessellation from closed rims (ellipse and/or circle) so the
  /// wall shares vertices with planar faces that sample the same curves.
  static bool _coneFromClosedLoops(
    BrepSolid solid,
    BrepCone s,
    BrepFace face,
    int fi,
    MeshTessellation tess,
    void Function(Vector3, Vector3, Vector3, int) addTri,
  ) {
    if (face.loops.isEmpty) return false;
    var hasEllipse = false;
    for (final loop in face.loops) {
      for (final c in loop.coedges) {
        if (c.edge < 0 || c.edge >= solid.edges.length) continue;
        if (solid.edges[c.edge].curve is BrepEllipse) hasEllipse = true;
      }
    }
    if (!hasEllipse) return false;
    final rings = [
      for (final loop in face.loops)
        sampleLoop(solid, loop, tess.cylinderSegments, tess: tess),
    ];
    if (rings.any((r) => r.length < 3)) return false;

    void emitFan(Vector3 tip, List<Vector3> ring) {
      for (var i = 0; i < ring.length; i++) {
        final a = ring[i];
        final b = ring[(i + 1) % ring.length];
        final mid = (a + b) * 0.5;
        final uv = s.parameters(mid) ?? (0.0, 0.0);
        var n = (a - tip).cross(b - tip);
        if (face.reversed) n = -n;
        if (n.dot(s.normalAt(uv.$1, uv.$2)) < 0) {
          addTri(tip, b, a, fi);
        } else {
          addTri(tip, a, b, fi);
        }
      }
    }

    if (rings.length == 1) {
      emitFan(s.origin + s.axis * s.height, rings.first);
      return true;
    }
    if (rings.length != 2) return false;
    var a = rings[0];
    var b = rings[1];
    if (a.length != b.length) return false;

    double uOf(Vector3 p) {
      var u = s.parameters(p)?.$1 ?? 0;
      if (u < 0) u += math.pi * 2;
      return u;
    }

    final u0 = uOf(a[0]);
    var best = 0;
    var bestD = 1e9;
    for (var i = 0; i < b.length; i++) {
      var d = (uOf(b[i]) - u0).abs();
      if (d > math.pi) d = math.pi * 2 - d;
      if (d < bestD) {
        bestD = d;
        best = i;
      }
    }
    if (best != 0) {
      b = [...b.sublist(best), ...b.sublist(0, best)];
    }
    var duA = uOf(a[1]) - uOf(a[0]);
    var duB = uOf(b[1]) - uOf(b[0]);
    while (duA > math.pi) {
      duA -= math.pi * 2;
    }
    while (duA < -math.pi) {
      duA += math.pi * 2;
    }
    while (duB > math.pi) {
      duB -= math.pi * 2;
    }
    while (duB < -math.pi) {
      duB += math.pi * 2;
    }
    if (duA * duB < 0) {
      b = [b.first, ...b.sublist(1).reversed];
    }

    for (var i = 0; i < a.length; i++) {
      final i1 = (i + 1) % a.length;
      final p00 = a[i];
      final p10 = a[i1];
      final p01 = b[i];
      final p11 = b[i1];
      final mid = (p00 + p10 + p01 + p11) * 0.25;
      final uv = s.parameters(mid) ?? (0.0, 0.0);
      var n = (p10 - p00).cross(p01 - p00);
      if (n.length2 < 1e-20) {
        n = (p11 - p00).cross(p01 - p00);
      }
      if (face.reversed) n = -n;
      if (n.dot(s.normalAt(uv.$1, uv.$2)) < 0) {
        addTri(p00, p01, p11, fi);
        addTri(p00, p11, p10, fi);
      } else {
        addTri(p00, p10, p11, fi);
        addTri(p00, p11, p01, fi);
      }
    }
    return true;
  }

  /// Unroll a cone onto a plane before triangulating its trimmed region.
  /// Unlike angle/height coordinates, this chart remains well-conditioned at
  /// the apex and does not turn narrow UV triangles inside out in 3D.
  static bool _coneDevelopedPatch(
    BrepSolid solid,
    BrepCone cone,
    BrepFace face,
    int fi,
    MeshTessellation tess,
    void Function(Vector3, Vector3, Vector3, int) emit,
  ) {
    final slope = (cone.topRadius - cone.radius) / cone.height;
    if (slope.abs() < 1e-8 || face.loops.isEmpty) return false;
    final k = slope.abs() / math.sqrt(1 + slope * slope);
    final sampled = [
      for (final loop in face.loops)
        sampleLoop(solid, loop, tess.cylinderSegments, tess: tess),
    ];
    final charts = [
      for (final points in sampled) axialLoopParameters(cone, points),
    ];
    final angles = charts.first.map((p) => p.x).toList()..sort();
    if (angles.isEmpty || angles.last - angles.first > 2 * math.pi + 1e-7)
      return false;
    final anchor = (angles.first + angles.last) / 2;
    final points = <Vector3>[], world = <Vector3>[], rings = <List<int>>[];
    // Resolve angles in the same branch but retain exact shared edge samples.
    for (var j = 0; j < sampled.length; j++) {
      final ring = <int>[];
      for (final p in sampled[j]) {
        final uv = cone.parameters(p)!;
        var u = uv.$1;
        u += ((anchor - u) / (2 * math.pi)).round() * 2 * math.pi;
        final r = cone.radiusAt(uv.$2);
        if (r < -1e-8) return false;
        final theta = k * (u - anchor), distance = r / k;
        final flat = Vector3(
          distance * math.cos(theta),
          distance * math.sin(theta),
          0,
        );
        if (ring.isNotEmpty && (points[ring.last] - flat).length < 1e-9)
          continue;
        ring.add(points.length);
        points.add(flat);
        world.add(p);
      }
      if (ring.length > 1 &&
          (points[ring.first] - points[ring.last]).length < 1e-9)
        ring.removeLast();
      if (ring.length < 3) return false;
      rings.add(ring);
    }
    var triangles = MeshPolygon.triangulateRings(
      points,
      rings,
      Vector3(0, 0, 1),
      keepCollinear: true,
    );
    if (triangles.isEmpty) return false;
    (int, int) key(int a, int b) => a < b ? (a, b) : (b, a);
    final boundary = {
      for (final ring in rings)
        for (var i = 0; i < ring.length; i++)
          key(ring[i], ring[(i + 1) % ring.length]),
    };
    final originals = points.length;
    Vector3 evaluate(Vector3 p) => cone.evaluate(
      anchor + math.atan2(p.y, p.x) / k,
      (p.length * k - cone.radius) / slope,
    );
    final error =
        math.max(cone.radius.abs(), cone.topRadius.abs()) *
        (1 - math.cos(math.pi / tess.cylinderSegments));
    _improveChartTriangles(points, triangles);
    for (var pass = 0; pass < 10; pass++) {
      final splits = <(int, int), int>{};
      for (final t in triangles) {
        for (var i = 0; i < 3; i++) {
          final a = t[i], b = t[(i + 1) % 3], edge = key(a, b);
          if (boundary.contains(edge) || splits.containsKey(edge)) continue;
          final mid = (points[a] + points[b]) * .5, p = evaluate(mid);
          if (!(a < originals && b < originals) &&
              (p - (world[a] + world[b]) * .5).length <= error)
            continue;
          splits[edge] = points.length;
          points.add(mid);
          world.add(p);
        }
      }
      if (splits.isEmpty) break;
      final next = <List<int>>[];
      for (final t in triangles) {
        var pieces = [t];
        for (var i = 0; i < 3; i++) {
          final a = t[i], b = t[(i + 1) % 3], mid = splits[key(a, b)];
          if (mid == null) continue;
          final updated = <List<int>>[];
          for (final piece in pieces) {
            if (!piece.contains(a) || !piece.contains(b)) {
              updated.add(piece);
              continue;
            }
            final c = piece.firstWhere((v) => v != a && v != b);
            updated.addAll([
              [a, mid, c],
              [mid, b, c],
            ]);
          }
          pieces = updated;
        }
        next.addAll(pieces);
      }
      triangles = next;
      _improveChartTriangles(points, triangles);
    }
    for (final t in triangles) {
      final a = world[t[0]], b = world[t[1]], c = world[t[2]];
      final uv = cone.parameters((a + b + c) / 3)!;
      if ((b - a).cross(c - a).dot(face.outwardNormal(uv.$1, uv.$2)) < 0)
        emit(a, c, b, fi);
      else
        emit(a, b, c, fi);
    }
    return true;
  }

  static bool _coneApexPatch(
    BrepSolid solid,
    BrepCone s,
    BrepFace face,
    int fi,
    MeshTessellation tess,
    void Function(Vector3, Vector3, Vector3, int) emit,
  ) {
    if (face.loops.length != 1) return false;
    final points = sampleLoop(
      solid,
      face.loops.single,
      tess.cylinderSegments,
      tess: tess,
    );
    final tipIndex = points.indexWhere(
      (p) => s.radiusAt(s.parameters(p)!.$2).abs() < 1e-8,
    );
    if (tipIndex < 0) return false;
    final tip = points[tipIndex],
        rim = [...points.skip(tipIndex + 1), ...points.take(tipIndex)];
    if (rim.length < 2) return false;
    var previous = s.parameters(rim.first)!.$1, direction = 0;
    for (final p in rim.skip(1)) {
      var u = s.parameters(p)!.$1;
      while (u - previous > math.pi) u -= 2 * math.pi;
      while (u - previous < -math.pi) u += 2 * math.pi;
      final delta = u - previous;
      // Only a rim that is a single-valued function of angle is visible
      // from the apex. Reentrant trims need the general parameter mesher.
      if (delta.abs() < 1e-10) {
        previous = u;
        continue;
      }
      final sign = delta > 0 ? 1 : -1;
      if (direction != 0 && sign != direction) return false;
      direction = sign;
      previous = u;
    }
    var triangles = <List<Vector3>>[];
    for (var i = 1; i < rim.length; i++) {
      final a = rim[i - 1], b = rim[i];
      if ((a - tip).cross(b - tip).length2 > 1e-24) triangles.add([tip, a, b]);
    }
    // A stepped rim can use the same generator at two heights. Subdivide
    // that generator on both neighboring triangles before welding the tip.
    for (final p in rim) {
      final next = <List<Vector3>>[];
      for (final tri in triangles) {
        var split = false;
        for (var k = 0; k < 3; k++) {
          final a = tri[k],
              b = tri[(k + 1) % 3],
              c = tri[(k + 2) % 3],
              d = b - a;
          if (d.length2 < 1e-24) continue;
          final t = (p - a).dot(d) / d.length2;
          if (t <= 1e-9 || t >= 1 - 1e-9 || (a + d * t - p).length > 1e-7)
            continue;
          next.addAll([
            [a, p, c],
            [p, b, c],
          ]);
          split = true;
          break;
        }
        if (!split) next.add(tri);
      }
      triangles = next;
    }
    for (final tri in triangles) {
      final a = tri[0],
          b = tri[1],
          c = tri[2],
          uv = s.parameters((a + b + c) / 3)!;
      if ((b - a).cross(c - a).dot(face.outwardNormal(uv.$1, uv.$2)) < 0)
        emit(a, c, b, fi);
      else
        emit(a, b, c, fi);
    }
    return true;
  }

  static void _cone(
    BrepSolid solid,
    BrepCone s,
    BrepFace face,
    int fi,
    MeshTessellation tess,
    void Function(Vector3, Vector3, Vector3, int) addTri,
  ) {
    if (face.loops.any(
      (l) => l.coedges.any(
        (c) =>
            solid.edges[c.edge].curve is BrepHelix ||
            solid.edges[c.edge].curve is BrepRuledQuadricCurve ||
            solid.edges[c.edge].curve is BrepPlaneSectionCurve ||
            solid.edges[c.edge].curve is BrepSurfaceIntersectionCurve,
      ),
    )) {
      if (_coneDevelopedPatch(solid, s, face, fi, tess, addTri)) return;
      if (_coneApexPatch(solid, s, face, fi, tess, addTri)) return;
      _axialTrimmed(solid, face, fi, tess, addTri);
      return;
    }
    if (face.loops.length == 1 &&
        face.loops.single.coedges.length == 1 &&
        s.topRadius.abs() < 1e-9) {
      final curve = solid.edges[face.loops.single.coedges.single.edge].curve;
      if (curve is BrepCircle &&
          curve.isClosed &&
          curve.normal.dot(s.axis).abs() > 1 - 1e-9) {
        final rim = sampleLoop(
          solid,
          face.loops.single,
          tess.cylinderSegments,
          tess: tess,
        );
        final tip = s.origin + s.axis * s.height;
        for (var i = 0; i < rim.length; i++) {
          final a = rim[i],
              b = rim[(i + 1) % rim.length],
              uv = s.parameters((a + b + tip) / 3)!;
          if ((b - a).cross(tip - a).dot(face.outwardNormal(uv.$1, uv.$2)) < 0)
            addTri(a, tip, b, fi);
          else
            addTri(a, b, tip, fi);
        }
        return;
      }
    }
    // STEP can place the cone apex at either parameter end (or outside
    // its nominal unit-height representation). Inspect the trim itself.
    if (face.loops.length == 1 &&
        face.loops.single.coedges.length > 1 &&
        _coneApexPatch(solid, s, face, fi, tess, addTri))
      return;
    if (_coneFromClosedLoops(solid, s, face, fi, tess, addTri)) return;
    // Trimmed frusta must share their actual rim samples with adjoining
    // faces. An independent angular grid leaves cracks at those rims.
    final domain = coneFaceDomain(solid, face);
    final slope = (s.topRadius - s.radius) / s.height;
    if (domain.$2 - domain.$1 <= math.pi + 1e-8 &&
        (s.radius + slope * domain.$3).abs() > 1e-8 &&
        (s.radius + slope * domain.$4).abs() > 1e-8) {
      _axialTrimmed(solid, face, fi, tess, addTri);
      return;
    }

    var u0 = 0.0, u1 = math.pi * 2, v0 = 0.0, v1 = s.height;
    if (face.loops.isNotEmpty) {
      final us = <double>[];
      var gotV = false;
      for (final loop in face.loops) {
        for (final co in loop.coedges) {
          if (co.edge < 0 || co.edge >= solid.edges.length) continue;
          final curve = solid.edges[co.edge].curve;
          final n = curve is BrepCircle ? (curve.isClosed ? 24 : 8) : 2;
          for (var i = 0; i <= n; i++) {
            final t = curve.t0 + (curve.t1 - curve.t0) * (i / n);
            final p = curve.evaluate(t);
            final uv = s.parameters(p);
            if (uv == null) continue;
            us.add(uv.$1);
            if (!gotV) {
              v0 = v1 = uv.$2;
              gotV = true;
            } else {
              v0 = math.min(v0, uv.$2);
              v1 = math.max(v1, uv.$2);
            }
          }
        }
      }
      if (us.isNotEmpty) {
        final span = _periodicSpan(us);
        u0 = span.$1;
        u1 = span.$2;
      }
      if (v1 - v0 < 1e-8) {
        v0 = 0;
        v1 = s.height;
      }
    }
    final nU = math.max(4, tess.cylinderSegments);
    final nV = math.max(2, nU ~/ 4);
    final uSpan = (u1 - u0).abs();
    final nu = math.max(4, (nU * uSpan / (math.pi * 2)).ceil());
    for (var i = 0; i < nu; i++) {
      final ua = u0 + (u1 - u0) * (i / nu);
      final ub = u0 + (u1 - u0) * ((i + 1) / nu);
      for (var j = 0; j < nV; j++) {
        final va = v0 + (v1 - v0) * (j / nV);
        final vb = v0 + (v1 - v0) * ((j + 1) / nV);
        final a = s.evaluate(ua, va);
        final b = s.evaluate(ub, va);
        final c = s.evaluate(ub, vb);
        final d = s.evaluate(ua, vb);
        if (face.reversed) {
          addTri(a, d, c, fi);
          addTri(a, c, b, fi);
        } else {
          addTri(a, b, c, fi);
          addTri(a, c, d, fi);
        }
      }
    }
  }
}
