import 'dart:math' as math;
import 'package:vector_math/vector_math_64.dart';
import 'geom.dart';
import 'solid.dart';
import 'tessellate.dart';

/// A swept outline can visit the same point more than once, or two openings
/// can touch there. Keep separate edges for the occupied material sectors,
/// while sharing their endpoint vertices. No surface or boundary is moved.
BrepSolid pairPlanarLineContacts(BrepSolid source) {
  final positions = <(int, int, int), int>{};
  final vertex = [
    for (var i = 0; i < source.vertices.length; i++)
      positions.putIfAbsent((
        (source.vertices[i].x * 1e9).round(),
        (source.vertices[i].y * 1e9).round(),
        (source.vertices[i].z * 1e9).round(),
      ), () => i),
  ];
  final uses = <int, List<int>>{};
  for (var fi = 0; fi < source.faces.length; fi++) {
    for (final loop in source.faces[fi].loops) {
      for (final co in loop.coedges) (uses[co.edge] ??= []).add(fi);
    }
  }
  final groups = <(int, int), List<int>>{};
  for (final entry in uses.entries) {
    final e = source.edges[entry.key];
    if (e.curve is! BrepLine ||
        entry.value.length != 2 ||
        entry.value.any((fi) => !_lineContactSurface(source.faces[fi].surface)))
      continue;
    final a = vertex[e.v0], b = vertex[e.v1];
    if (a == b) continue;
    (groups[a < b ? (a, b) : (b, a)] ??= []).add(entry.key);
  }
  final replacements = <(int, int), (int, bool)>{};
  final added = <BrepEdge>[];
  for (final group in groups.values.where((g) => g.length > 1)) {
    final faceEdges = <int, int>{};
    for (final eid in group) {
      for (final fi in uses[eid]!) faceEdges[fi] = eid;
    }
    if (faceEdges.length != group.length * 2) continue;
    final pairs = planarContactFacePairs(source, faceEdges);
    if (pairs == null) continue;
    final e = source.edges[group.first];
    for (final pair in pairs) {
      final id = source.edges.length + added.length;
      added.add(
        BrepEdge(
          v0: vertex[e.v0],
          v1: vertex[e.v1],
          curve: e.curve,
          origin: e.origin,
        ),
      );
      for (final fi in pair) {
        final old = faceEdges[fi]!;
        replacements[(fi, old)] = (
          id,
          vertex[source.edges[old].v0] != vertex[e.v0],
        );
      }
    }
  }
  if (added.isEmpty) return source;
  return BrepSolid(
    vertices: source.vertices,
    edges: [
      for (final e in source.edges)
        BrepEdge(
          v0: vertex[e.v0],
          v1: vertex[e.v1],
          curve: e.curve,
          origin: e.origin,
        ),
      ...added,
    ],
    faces: [
      for (var fi = 0; fi < source.faces.length; fi++)
        BrepFace(
          surface: source.faces[fi].surface,
          reversed: source.faces[fi].reversed,
          origin: source.faces[fi].origin,
          loops: [
            for (final loop in source.faces[fi].loops)
              BrepLoop([
                for (final co in loop.coedges)
                  if (replacements[(fi, co.edge)] case final replacement?)
                    BrepCoedge(
                      replacement.$1,
                      reversed: co.reversed != replacement.$2,
                    )
                  else
                    co,
              ]),
          ],
        ),
    ],
  );
}

bool _lineContactSurface(BrepSurface s) =>
    s is BrepPlane || s is BrepRuledSurface;

/// Pair planar sheets and ruled generator walls around a straight contact by
/// their occupied angular sectors. Curved walls use their normal at the rail.
/// Face connectivity alone cannot resolve a contact when the two sheets join
/// elsewhere (for example, a tongue underneath a trailer frame).
List<List<int>>? planarContactFacePairs(BrepSolid s, Map<int, int> faceEdges) {
  if (faceEdges.length < 4) return null;
  final edge = s.edges[faceEdges.values.first];
  if (edge.curve is! BrepLine || faceEdges.length.isOdd) return null;
  final axis = s.vertices[edge.v1] - s.vertices[edge.v0];
  if (axis.length2 < 1e-16) return null;
  axis.normalize();
  var ref = axis.cross(Vector3(0, 0, 1));
  if (ref.length2 < 1e-12) ref = axis.cross(Vector3(0, 1, 0));
  ref.normalize();
  final bin = axis.cross(ref);
  final rays = <(double, int, bool)>[];
  for (final fi in faceEdges.keys) {
    final eid = faceEdges[fi]!;
    final face = s.faces[fi];
    final surface = face.surface;
    if (!_lineContactSurface(surface)) return null;
    final loops = face.loops.where((l) => l.coedges.any((c) => c.edge == eid));
    if (loops.length != 1) return null;
    final loop = loops.single;
    final coedges = loop.coedges.where((c) => c.edge == eid);
    if (coedges.length != 1) return null;
    var u = 0.0, v = 0.0;
    if (surface is BrepRuledSurface) {
      final a = s.vertices[s.edges[eid].v0], b = s.vertices[s.edges[eid].v1];
      final uv = surface.parameters((a + b) * .5);
      u = uv.$1;
      v = uv.$2;
      // Only an actual straight generator on a ruled patch is eligible.
      final pa = surface.parameters(a), pb = surface.parameters(b);
      if ((pa.$1 - pb.$1).abs() > 1e-7 ||
          (surface.evaluate(u, v) - (a + b) * .5).length > 1e-7 ||
          (surface.dv(u, v)).cross(axis).length > 1e-7)
        return null;
    }
    final n = face.outwardNormal(u, v).normalized();
    if (n.dot(axis).abs() > 1e-7) return null;
    final ownEdge = s.edges[eid];
    var tangent = (s.vertices[ownEdge.v1] - s.vertices[ownEdge.v0])
        .normalized();
    if (coedges.single.reversed) tangent = -tangent;
    if (surface is BrepPlane) {
      final points = BrepTessellate.sampleLoop(s, loop, 32);
      var area = Vector3.zero();
      for (var i = 0; i < points.length; i++) {
        area += (points[i] - surface.origin).cross(
          points[(i + 1) % points.length] - surface.origin,
        );
      }
      if (area.dot(n).abs() < 1e-12) return null;
      if ((area.dot(n) < 0) == identical(loop, face.loops.first)) {
        tangent = -tangent;
      }
    } else if (face.reversed) {
      tangent = -tangent;
    }
    final ray = n.cross(tangent).normalized();
    rays.add((
      math.atan2(ray.dot(bin), ray.dot(ref)),
      fi,
      n.dot(axis.cross(ray)) < 0,
    ));
  }
  rays.sort((a, b) => a.$1.compareTo(b.$1));
  final pairs = <List<int>>[];
  for (var i = 0; i < rays.length; i++) {
    final a = rays[i], b = rays[(i + 1) % rays.length];
    final gap = (b.$1 - a.$1 + 2 * math.pi) % (2 * math.pi);
    // Coincident sheets or inconsistent orientations remain an error.
    if (gap < 1e-7 || a.$3 == b.$3) return null;
    if (a.$3) pairs.add([a.$2, b.$2]);
  }
  return pairs.length * 2 == faceEdges.length ? pairs : null;
}

/// Duplicate topological line edges identify sheets touching at a line while
/// keeping separate face pairs. Rebuilding these with coplanar face merging
/// would erase that distinction.
bool hasPairedPlanarContacts(BrepSolid s) {
  final edgeFaces = <int, List<BrepFace>>{};
  for (final f in s.faces) {
    for (final l in f.loops) {
      for (final c in l.coedges) (edgeFaces[c.edge] ??= []).add(f);
    }
  }
  final seen = <(int, int)>{};
  for (final entry in edgeFaces.entries) {
    if (entry.value.length != 2 ||
        entry.value.any((f) => !_lineContactSurface(f.surface)))
      continue;
    final e = s.edges[entry.key];
    if (e.curve is! BrepLine || e.v0 == e.v1) continue;
    final key = e.v0 < e.v1 ? (e.v0, e.v1) : (e.v1, e.v0);
    if (!seen.add(key)) return true;
  }
  return false;
}
