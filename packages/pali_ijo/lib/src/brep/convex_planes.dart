import 'dart:math' as math;
import 'package:vector_math/vector_math_64.dart';
import 'geom.dart';
import 'solid.dart';
import 'heal.dart';
import 'origin.dart';

/// Exact planar clipping for closed convex polyhedra. No tessellated geometry
/// is used: intersections become line edges and each trim retains its plane.
class BrepConvexPlanes {
  BrepConvexPlanes._();

  static bool supports(BrepSolid s) {
    if (s.vertices.length < 4 || !BrepHeal.topologyIsClosed(s)) return false;
    if (s.edges.any((e) => e.curve is! BrepLine)) return false;
    for (final f in s.faces) {
      if (f.surface is! BrepPlane || f.loops.length != 1) return false;
      final p = f.surface as BrepPlane, n = f.outwardNormal(0, 0);
      if (s.vertices.any((v) => (v - p.origin).dot(n) > 1e-7)) return false;
    }
    return true;
  }

  /// Keep the negative side of every cut. [sources] tags new faces with the
  /// original face index for shell openings; -1 denotes a new chamfer face.
  static ({BrepSolid solid, Map<int, int> faceMap}) clip(
    BrepSolid s,
    List<BrepPlane> cuts, {
    List<int>? sources,
  }) {
    const eps = 1e-8;
    var polygons = <_PlanePolygon>[];
    for (var fi = 0; fi < s.faces.length; fi++) {
      final f = s.faces[fi];
      final points = [
        for (final c in f.loops.single.coedges)
          s.vertices[c.reversed ? s.edges[c.edge].v1 : s.edges[c.edge].v0]
              .clone(),
      ];
      final p = f.surface as BrepPlane;
      final plane = BrepPlane(
        origin: p.origin,
        normal: f.outwardNormal(0, 0),
        xDir: p.xDir,
      );
      _orient(points, plane.normal);
      polygons.add(_PlanePolygon(plane, points, fi, f.origin));
    }
    for (var ci = 0; ci < cuts.length; ci++) {
      final plane = cuts[ci];
      final next = <_PlanePolygon>[], rim = <Vector3>[];
      var removed = false;
      void addRim(Vector3 p) {
        if (!rim.any((q) => q.distanceTo(p) < eps)) rim.add(p);
      }

      for (final polygon in polygons) {
        final out = <Vector3>[];
        void add(Vector3 p) {
          if (out.isEmpty || out.last.distanceTo(p) > eps) out.add(p);
        }

        final pts = polygon.points;
        for (var i = 0; i < pts.length; i++) {
          final a = pts[i], b = pts[(i + 1) % pts.length];
          final da = plane.height(a), db = plane.height(b);
          final aIn = da <= eps, bIn = db <= eps;
          if (!aIn) removed = true;
          if (aIn) {
            add(a);
            if (da.abs() <= eps) addRim(a);
          }
          if (aIn != bIn) {
            // Snap endpoints already on the cut rather than making sliver edges.
            final p = da.abs() <= eps
                ? a
                : db.abs() <= eps
                ? b
                : a + (b - a) * (da / (da - db));
            add(p);
            addRim(p);
          }
        }
        if (out.length > 1 && out.first.distanceTo(out.last) <= eps) {
          out.removeLast();
        }
        if (out.length >= 3 && _area(out).length > eps * eps) {
          next.add(
            _PlanePolygon(polygon.plane, out, polygon.source, polygon.origin),
          );
        }
      }
      if (!removed) continue;
      if (next.isEmpty || rim.length < 3) {
        throw BrepUnsupported(
          'Offset or chamfer consumes the remaining interior',
        );
      }
      final center = rim.reduce((a, b) => a + b) / rim.length.toDouble();
      rim.sort((a, b) {
        final av = a - center, bv = b - center;
        return math
            .atan2(av.dot(plane.yDir), av.dot(plane.xDir))
            .compareTo(math.atan2(bv.dot(plane.yDir), bv.dot(plane.xDir)));
      });
      _orient(rim, plane.normal);
      final source = sources?[ci] ?? -1;
      next.add(
        _PlanePolygon(
          plane,
          rim,
          source,
          source < 0 ? null : s.faces[source].origin,
        ),
      );
      polygons = next;
    }
    final b = BrepBuilder(),
        edges = <(int, int), int>{},
        faceMap = <int, int>{};
    int vertex(Vector3 p) {
      for (var i = 0; i < b.vertices.length; i++) {
        if (b.vertices[i].distanceTo(p) < eps) return i;
      }
      return b.vertex(p);
    }

    for (final poly in polygons) {
      final ids = poly.points.map(vertex).toList(), coedges = <BrepCoedge>[];
      for (var i = 0; i < ids.length; i++) {
        final a = ids[i], z = ids[(i + 1) % ids.length];
        if (a == z) continue;
        final key = (math.min(a, z), math.max(a, z));
        final ei = edges.putIfAbsent(
          key,
          () => b.edge(a, z, BrepLine(b.vertices[a], b.vertices[z])),
        );
        coedges.add(BrepCoedge(ei, reversed: b.edges[ei].v0 != a));
      }
      if (coedges.length < 3) continue;
      if (poly.source >= 0) faceMap[poly.source] = b.faces.length;
      b.face(poly.plane, [BrepLoop(coedges)], origin: poly.origin);
    }
    final result = b.build();
    if (!BrepHeal.topologyIsClosed(result)) {
      throw BrepUnsupported('Planar clipping could not form a closed boundary');
    }
    return (solid: result, faceMap: faceMap);
  }

  static Vector3 _area(List<Vector3> p) {
    var area = Vector3.zero();
    for (var i = 1; i + 1 < p.length; i++) {
      area += (p[i] - p[0]).cross(p[i + 1] - p[0]);
    }
    return area;
  }

  static void _orient(List<Vector3> p, Vector3 n) {
    if (_area(p).dot(n) < 0) {
      final reversed = p.reversed.toList();
      p
        ..clear()
        ..addAll(reversed);
    }
  }
}

class _PlanePolygon {
  _PlanePolygon(this.plane, this.points, this.source, this.origin);
  final BrepPlane plane;
  final List<Vector3> points;
  final int source;
  final BrepOrigin? origin;
}
