import 'dart:math' as math;
import 'package:vector_math/vector_math_64.dart';

import '../mesh/mesh_data.dart';
import '../ijo_log.dart';
import '../mesh/mesh_edit.dart';
import '../mesh/mesh_seal.dart';
import '../tessellation.dart';
import 'geom.dart';
import 'provenance.dart';
import 'solid.dart';
import 'contact.dart';
import 'tessellate.dart';

/// Repair a solid: merge coplanar splits, drop collinear seams, outward
/// winding, closed tessellation.
class BrepHeal {
  BrepHeal._();

  /// Merges coplanar faces, collapses collinear edges and orients outward.
  ///
  /// Provenance survives. Healing is a repair pass, not a modelling step: it
  /// changes how the same surfaces are divided up, so nothing here should ever
  /// come back marked as newly created. Since every face it emits lies on a
  /// surface that was already present, geometric matching recovers all of it.
  static BrepSolid run(BrepSolid s) {
    if (s.faces.isEmpty) return s;
    var out = mergeCoplanar(sewClosedRims(s));
    out = collapseCollinear(out);
    out = mergeCoplanar(out);
    out = splitTJunctions(out);
    out = splitPlanarIslands(out);
    out = mergeCoplanar(out);
    out = _orientOutward(out);
    return BrepProvenance.inheritFrom(out, s, op: 'heal');
  }

  /// Clipping a concave polygon can leave separate patches joined by an edge
  /// traversed in both directions. Such a zero-width bridge has no surface,
  /// but retaining it makes face adjacency report disconnected bodies as one.
  /// Run after T-junction splitting so overlapping bridge segments share edges.
  /// Also separates islands nested inside holes (such as a solid head closing
  /// a tube). Sampling determines containment only; exact curves are retained.
  static BrepSolid splitPlanarIslands(BrepSolid s) {
    final faces = <BrepFace>[];
    var changed = false;
    for (final face in s.faces) {
      if (face.surface is! BrepPlane) {
        faces.add(face);
        continue;
      }
      var split = false;
      final pending = face.loops.map((l) => l.coedges).toList();
      final rings = <BrepLoop>[];
      while (pending.isNotEmpty) {
        final ring = pending.removeLast();
        final seen = <int, int>{};
        var cut = false;
        for (var j = 0; j < ring.length; j++) {
          final co = ring[j], i = seen[co.edge];
          if (i != null && ring[i].reversed != co.reversed) {
            pending.add(ring.sublist(i + 1, j));
            pending.add([...ring.sublist(j + 1), ...ring.sublist(0, i)]);
            cut = split = true;
            break;
          }
          seen[co.edge] = j;
        }
        if (!cut && ring.isNotEmpty) rings.add(BrepLoop(ring));
      }
      if ((!split && rings.length < 2) || rings.isEmpty) {
        faces.add(face);
        continue;
      }
      final plane = face.surface as BrepPlane;
      final polygons = [
        for (final ring in rings)
          [
            for (final point in BrepTessellate.sampleLoop(s, ring, 128))
              plane.parameters(point)!,
          ],
      ];
      if (polygons.any((p) => p.length < 3)) {
        faces.add(face);
        continue;
      }
      double area(List<(double, double)> polygon) {
        var sum = 0.0;
        for (var i = 0; i < polygon.length; i++) {
          final a = polygon[i], b = polygon[(i + 1) % polygon.length];
          sum += a.$1 * b.$2 - b.$1 * a.$2;
        }
        return sum.abs();
      }

      bool contains(List<(double, double)> polygon, (double, double) p) {
        var inside = false;
        for (var i = 0; i < polygon.length; i++) {
          final a = polygon[i], b = polygon[(i + 1) % polygon.length];
          if ((a.$2 > p.$2) != (b.$2 > p.$2) &&
              p.$1 < a.$1 + (p.$2 - a.$2) * (b.$1 - a.$1) / (b.$2 - a.$2)) {
            inside = !inside;
          }
        }
        return inside;
      }

      // A bridge can also connect an outer wire to a hole. Containment keeps
      // those holes on their enclosing patch rather than filling them in.
      final areas = polygons.map(area).toList();
      final parents = List<int?>.filled(rings.length, null);
      for (var i = 0; i < rings.length; i++) {
        for (var j = 0; j < rings.length; j++) {
          if (areas[j] > areas[i] &&
              (parents[i] == null || areas[j] < areas[parents[i]!]) &&
              contains(polygons[j], polygons[i].first)) {
            parents[i] = j;
          }
        }
      }
      int depth(int i) => parents[i] == null ? 0 : 1 + depth(parents[i]!);
      if (!split &&
          List.generate(rings.length, depth).where((d) => d.isEven).length ==
              1) {
        faces.add(face);
        continue;
      }
      for (var i = 0; i < rings.length; i++) {
        if (depth(i).isOdd) continue;
        faces.add(
          BrepFace(
            surface: face.surface,
            loops: [
              rings[i],
              for (var j = 0; j < rings.length; j++)
                if (parents[j] == i) rings[j],
            ],
            reversed: face.reversed,
            origin: face.origin,
          ),
        );
      }
      changed = true;
    }
    if (!changed) return s;
    // Compact the topology: cancelled bridge edges must not remain as unused
    // edges in an otherwise closed solid. Keep ancestry on every surviving use.
    final b = BrepBuilder(), vertices = <int, int>{}, edges = <int, int>{};
    int vertex(int i) => vertices.putIfAbsent(i, () => b.vertex(s.vertices[i]));
    int edge(int i) => edges.putIfAbsent(i, () {
      final e = s.edges[i];
      return b.edge(vertex(e.v0), vertex(e.v1), e.curve, origin: e.origin);
    });
    for (final face in faces) {
      b.face(
        face.surface,
        [
          for (final loop in face.loops)
            BrepLoop([
              for (final c in loop.coedges)
                BrepCoedge(edge(c.edge), reversed: c.reversed),
            ]),
        ],
        reversed: face.reversed,
        origin: face.origin,
      );
    }
    return b.build();
  }

  /// Recover a straight extrusion whose faceted walls were incorrectly
  /// merged onto a single plane. Both cap wires and their translated edge
  /// pairs must agree exactly; arbitrary nonplanar faces are left untouched.
  static BrepSolid restorePrismaticWalls(BrepSolid s) {
    if (!s.faces.every((f) => f.surface is BrepPlane) ||
        !s.edges.every((e) => e.curve is BrepLine))
      return s;
    bool planar(BrepFace f) => f.loops.every(
      (l) => l.coedges.every((co) {
        final e = s.edges[co.edge], p = f.surface as BrepPlane;
        return p.height(s.vertices[e.v0]).abs() < 1e-7 &&
            p.height(s.vertices[e.v1]).abs() < 1e-7;
      }),
    );
    if (s.faces.every(planar)) return s;
    Set<int> edgeIds(BrepFace f) => {
      for (final l in f.loops)
        for (final c in l.coedges) c.edge,
    };
    final faceEdges = [for (final f in s.faces) edgeIds(f)];
    for (var ci = 0; ci < s.faces.length; ci++) {
      final cap = s.faces[ci], plane = cap.surface as BrepPlane;
      if (!planar(cap)) continue;
      for (var ti = ci + 1; ti < s.faces.length; ti++) {
        final top = s.faces[ti], upper = top.surface as BrepPlane;
        if (!planar(top) ||
            cap.outwardNormal(0, 0).dot(top.outwardNormal(0, 0)) > -1 + 1e-10 ||
            plane.height(upper.origin).abs() < 1e-7)
          continue;
        if (s.vertices.any(
          (v) => plane.height(v).abs() > 1e-7 && upper.height(v).abs() > 1e-7,
        ))
          continue;
        Vector3? direction;
        for (final e in s.edges) {
          final a = s.vertices[e.v0], b = s.vertices[e.v1];
          if (plane.height(a).abs() < 1e-7 && upper.height(b).abs() < 1e-7) {
            direction = b - a;
            break;
          }
          if (plane.height(b).abs() < 1e-7 && upper.height(a).abs() < 1e-7) {
            direction = a - b;
            break;
          }
        }
        if (direction == null) continue;
        final translated = <int, int>{};
        for (final eid in faceEdges[ci]) {
          final e = s.edges[eid];
          for (final vi in [e.v0, e.v1]) {
            final target = s.vertices[vi] + direction;
            final match = s.vertices.indexWhere(
              (v) => (v - target).length < 1e-7,
            );
            if (match >= 0) translated[vi] = match;
          }
        }
        final pairs = <int, (int, int)>{};
        for (final eid in faceEdges[ci]) {
          final e = s.edges[eid], p = translated[e.v0], q = translated[e.v1];
          if (p == null || q == null) break;
          final matches = faceEdges[ti].where((id) {
            final t = s.edges[id];
            return t.v0 == p && t.v1 == q || t.v0 == q && t.v1 == p;
          }).toList();
          if (matches.length != 1) break;
          final walls = [
            for (var fi = 0; fi < s.faces.length; fi++)
              if (fi != ci &&
                  fi != ti &&
                  faceEdges[fi].contains(eid) &&
                  faceEdges[fi].contains(matches.single))
                fi,
          ];
          if (walls.length != 1) break;
          pairs[eid] = (matches.single, walls.single);
        }
        if (pairs.length != faceEdges[ci].length ||
            pairs.values.map((p) => p.$1).toSet().length !=
                faceEdges[ti].length ||
            pairs.values.map((p) => p.$2).toSet().length != s.faces.length - 2)
          continue;
        final builder = BrepBuilder();
        for (final v in s.vertices) builder.vertex(v);
        for (final e in s.edges)
          builder.edge(e.v0, e.v1, e.curve, origin: e.origin);
        for (final f in [cap, top])
          builder.face(
            f.surface,
            f.loops,
            reversed: f.reversed,
            origin: f.origin,
          );
        final generators = <int, int>{};
        int generator(int vi) => generators.putIfAbsent(
          vi,
          () => builder.edge(
            vi,
            translated[vi]!,
            BrepLine(s.vertices[vi], s.vertices[translated[vi]!]),
          ),
        );
        for (final l in cap.loops)
          for (final co in l.coedges) {
            final e = s.edges[co.edge], pair = pairs[co.edge]!;
            final p = co.reversed ? e.v1 : e.v0, q = co.reversed ? e.v0 : e.v1;
            final normal = (s.vertices[p] - s.vertices[q]).cross(direction);
            final te = s.edges[pair.$1];
            builder.face(
              BrepPlane(origin: s.vertices[p], normal: normal),
              [
                BrepLoop([
                  BrepCoedge(co.edge, reversed: !co.reversed),
                  BrepCoedge(generator(p)),
                  BrepCoedge(pair.$1, reversed: te.v0 != translated[p]),
                  BrepCoedge(generator(q), reversed: true),
                ]),
              ],
              reversed: cap.reversed,
              origin: s.faces[pair.$2].origin,
            );
          }
        final repaired = builder.build();
        if (topologyIsClosed(repaired)) return run(repaired);
      }
    }
    return s;
  }

  /// Rebind circular trim parameters to vertices already on the same circle.
  /// Offset junctions can slide along a circle; retaining the original angles
  /// leaves a geometrically open wire even when its edge indices are paired.
  /// This never moves vertices or changes the supporting circle.
  static BrepSolid alignCircularTrims(BrepSolid solid) {
    BrepSolid? out;
    for (var i = 0; i < solid.edges.length; i++) {
      final edge = solid.edges[i], curve = edge.curve;
      if (curve is! BrepCircle || curve.isClosed) continue;
      final a = solid.vertices[edge.v0], b = solid.vertices[edge.v1];
      if ((curve.evaluate(curve.t0) - a).length < 1e-8 &&
          (curve.evaluate(curve.t1) - b).length < 1e-8)
        continue;
      bool onCircle(Vector3 p) {
        final d = p - curve.center;
        return d.dot(curve.normal).abs() < 1e-7 &&
            (d.length - curve.radius).abs() < 1e-7;
      }

      if (!onCircle(a) || !onCircle(b)) continue;
      double angle(Vector3 p, double near) {
        final d = p - curve.center;
        final t = math.atan2(d.dot(curve.yDir), d.dot(curve.xDir));
        return t + ((near - t) / (2 * math.pi)).round() * 2 * math.pi;
      }

      // Prefer parameters on the existing finite arc. Some saved fillet and
      // hollow edges have their vertices in the opposite order to the trim;
      // this also preserves the correct arc across the angle branch cut.
      var t0 = curve.closestParam(a), t1 = curve.closestParam(b);
      if ((curve.evaluate(t0) - a).length > 1e-7 ||
          (curve.evaluate(t1) - b).length > 1e-7) {
        // Offset junctions may extend an arc. Keep the original direction in
        // that case, rather than guessing a reversed complementary arc.
        t0 = angle(a, curve.t0);
        t1 = angle(b, curve.t1);
        if ((t1 - t0) * (curve.t1 - curve.t0) <= 0) continue;
      }
      if ((t1 - t0).abs() < 1e-12 || (t1 - t0).abs() >= 2 * math.pi) continue;
      out ??= solid.clone();
      out.edges[i].curve = BrepCircle(
        center: curve.center,
        normal: curve.normal,
        radius: curve.radius,
        xDir: curve.xDir,
        t0: t0,
        t1: t1,
      );
    }
    return out ?? solid;
  }

  /// Triangle soup → planar n-gons, then [run].
  static BrepSolid fromMesh(MeshData mesh) {
    var m = MeshSeal.unifyCoplanarFaces(mesh);
    m = MeshSeal.sealSolid(m);
    return run(_trianglesToBrep(m));
  }

  /// Remove a planar partition where one face exactly fills another face's
  /// inner trim. Keep edge indices stable so an existing selection survives.
  /// Real holes have no coplanar filling face and remain unchanged.
  static BrepSolid mergeFilledPlanarLoops(BrepSolid s) {
    final faces = s.faces.toList();
    var changed = false;
    var again = true;
    while (again) {
      again = false;
      outer:
      for (var i = 0; i < faces.length; i++) {
        final a = faces[i];
        if (a.surface is! BrepPlane || a.loops.length < 2) continue;
        final plane = a.surface as BrepPlane;
        for (var j = 0; j < faces.length; j++) {
          if (i == j) continue;
          final b = faces[j];
          if (b.surface is! BrepPlane ||
              b.loops.isEmpty ||
              a.outwardNormal(0, 0).dot(b.outwardNormal(0, 0)) < 1 - 1e-10 ||
              plane.height((b.surface as BrepPlane).origin).abs() > 1e-7) {
            continue;
          }
          final fill = b.loops.first.coedges.map((c) => c.edge).toSet();
          if (fill.length != b.loops.first.coedges.length) continue;
          for (var k = 1; k < a.loops.length; k++) {
            final hole = a.loops[k].coedges;
            if (hole.length != fill.length ||
                !hole.every((c) => fill.contains(c.edge))) {
              continue;
            }
            // Do not erase an edge shared with a third face or an overlapping
            // shell. A matching fill must be the only other use of each edge.
            if (fill.any(
              (id) =>
                  faces.fold<int>(
                    0,
                    (count, f) =>
                        count +
                        f.loops
                            .expand((l) => l.coedges)
                            .where((c) => c.edge == id)
                            .length,
                  ) !=
                  2,
            )) {
              continue;
            }
            faces[i] = BrepFace(
              surface: a.surface,
              loops: [
                for (var l = 0; l < a.loops.length; l++)
                  if (l != k) a.loops[l],
                ...b.loops.skip(1),
              ],
              reversed: a.reversed,
              origin: a.origin,
            );
            faces.removeAt(j);
            changed = again = true;
            break outer;
          }
        }
      }
    }
    return changed
        ? BrepSolid(vertices: s.vertices, edges: s.edges, faces: faces)
        : s;
  }

  /// Join adjacent coplanar faces that share an edge into one loop.
  /// Join adjacent planar patches, optionally restricting joins to one source
  /// face. [faceMap] receives every old-to-new face index; vertices and edges
  /// retain their indices so display tags can be updated without retessellation.
  static BrepSolid mergeCoplanar(
    BrepSolid s, {
    Map<int, int>? faceMap,
    bool sameOriginOnly = false,
  }) {
    var groups = [
      for (var i = 0; i < s.faces.length; i++) [i],
    ];
    var faces = [for (final f in s.faces) f];
    var changed = true;
    while (changed) {
      changed = false;
      outer:
      for (var i = 0; i < faces.length; i++) {
        final a = faces[i];
        if (a.surface is! BrepPlane || a.loops.length != 1) continue;
        for (var j = i + 1; j < faces.length; j++) {
          final b = faces[j];
          if (b.surface is! BrepPlane || b.loops.length != 1) continue;
          if (sameOriginOnly && (a.origin == null || a.origin != b.origin)) {
            continue;
          }
          final pa = a.surface as BrepPlane;
          final pb = b.surface as BrepPlane;
          // A reversed +Z plane and a forward -Z plane face the same way.
          // Booleans can retain both encodings on one physical surface.
          if (a.outwardNormal(0, 0).dot(b.outwardNormal(0, 0)) < 1 - 1e-10) {
            continue;
          }
          if (pa.height(pb.origin).abs() > 1e-7) continue;
          if (b.loops.any(
            (loop) => loop.coedges.any((co) {
              final edge = s.edges[co.edge];
              return pa.height(s.vertices[edge.v0]).abs() > 1e-7 ||
                  pa.height(s.vertices[edge.v1]).abs() > 1e-7;
            }),
          ))
            continue;
          final shared = _sharedEdge(a, b);
          if (shared == null) continue;
          // Tangent cylinder generators split a wall into two coplanar
          // halves. Merging them glues the corners back into one solid.
          if (_edgeUsedByNonPlane(faces, shared)) continue;
          final all = _sharedEdges(a, b);
          final merged =
              (all.length > 1 &&
                      !all.any(
                        (e) => e != shared && _edgeUsedByNonPlane(faces, e),
                      )
                  ? _mergeSharedRun(a, b, all)
                  : null) ??
              _mergeLoops(a, b, shared);
          if (merged == null) continue;
          faces = [
            for (var k = 0; k < faces.length; k++)
              if (k != i && k != j) faces[k],
            BrepFace(
              surface: pa,
              loops: [merged],
              reversed: a.reversed,
              origin: a.origin,
            ),
          ];
          groups = [
            for (var k = 0; k < groups.length; k++)
              if (k != i && k != j) groups[k],
            [...groups[i], ...groups[j]],
          ];
          changed = true;
          break outer;
        }
      }
    }
    if (faceMap != null) {
      for (var i = 0; i < groups.length; i++) {
        for (final old in groups[i]) {
          faceMap[old] = i;
        }
      }
    }
    if (faces.length == s.faces.length) return s;
    return BrepSolid(vertices: s.vertices, edges: s.edges, faces: faces);
  }

  /// Drop vertices that sit on a straight run of line edges (extrude waist).
  static BrepSolid collapseCollinear(BrepSolid s) {
    final b = BrepBuilder();
    int weld(Vector3 p) {
      for (var i = 0; i < b.vertices.length; i++) {
        if ((b.vertices[i] - p).length2 < 1e-16) return i;
      }
      return b.vertex(p);
    }

    final edgeOf = <String, int>{};
    (int, bool) sew(int i, int j, BrepCurve curve) {
      final lo = i <= j ? i : j;
      final hi = i <= j ? j : i;
      final String k;
      if (curve is BrepLine) {
        k = '$lo-$hi-line';
      } else {
        final mid = curve.evaluate(0.5 * (curve.t0 + curve.t1));
        final q =
            '${mid.x.toStringAsFixed(5)}_${mid.y.toStringAsFixed(5)}_${mid.z.toStringAsFixed(5)}';
        k = '$lo-$hi-$q';
      }
      final existing = edgeOf[k];
      if (existing != null) {
        final e = b.edges[existing];
        return (existing, e.v0 != i);
      }
      final id = b.edge(i, j, curve);
      edgeOf[k] = id;
      return (id, false);
    }

    var changed = false;
    for (final face in s.faces) {
      if (face.loops.isEmpty) {
        b.face(face.surface, const [], reversed: face.reversed);
        continue;
      }
      final loops = <BrepLoop>[];
      for (final loop in face.loops) {
        final segs = <_HealSeg>[];
        for (final c in loop.coedges) {
          final e = s.edges[c.edge];
          final a = s.vertices[c.reversed ? e.v1 : e.v0];
          final end = s.vertices[c.reversed ? e.v0 : e.v1];
          segs.add(_HealSeg(a, end, e.curve, reversed: c.reversed));
        }
        final collapsed = _collapseSegs(segs);
        if (collapsed.length != segs.length) changed = true;
        // A circle may have two arc edges, or a D profile one arc and one line.
        if (collapsed.length < 3 &&
            !collapsed.any((seg) => seg.curve is! BrepLine)) {
          continue;
        }
        final coedges = <BrepCoedge>[];
        for (final seg in collapsed) {
          final i0 = weld(seg.a);
          final i1 = weld(seg.b);
          final closedCurve =
              i0 == i1 &&
              seg.curve is BrepCircle &&
              (seg.curve as BrepCircle).isClosed;
          if (i0 == i1 && !closedCurve) continue;
          final reversed = seg.curve is! BrepLine && seg.reversed;
          final sewn = reversed
              ? sew(i1, i0, seg.curve)
              : sew(i0, i1, seg.curve);
          coedges.add(BrepCoedge(sewn.$1, reversed: sewn.$2 != reversed));
        }
        if (coedges.length >= 3 ||
            coedges.isNotEmpty &&
                collapsed.any((seg) => seg.curve is! BrepLine)) {
          loops.add(BrepLoop(coedges));
        }
      }
      if (loops.isEmpty && face.loops.isNotEmpty) continue;
      b.face(face.surface, loops, reversed: face.reversed);
    }
    if (!changed || b.faces.isEmpty) return s;
    return b.build();
  }

  /// Split line edges at vertices that sit in their interior. Overlapping-box
  /// unions leave a long wall edge beside a split neighbour (a T-junction);
  /// sewing those pieces shut is what makes the stepped solid manifold.
  static BrepSolid splitTJunctions(BrepSolid s) {
    if (s.faces.isEmpty || s.edges.isEmpty) return s;
    var needed = false;
    for (final e in s.edges) {
      if (e.curve is! BrepLine) continue;
      final a = s.vertices[e.v0];
      final b = s.vertices[e.v1];
      for (final v in s.vertices) {
        if (_pointOnSegmentInterior(v, a, b)) {
          needed = true;
          break;
        }
      }
      if (needed) break;
    }
    if (!needed) return s;

    final b = BrepBuilder();
    int weld(Vector3 p) {
      for (var i = 0; i < b.vertices.length; i++) {
        if ((b.vertices[i] - p).length < 1e-6) return i;
      }
      return b.vertex(p);
    }

    final edgeOf = <String, int>{};
    (int, bool) sew(int i, int j, BrepCurve curve) {
      final lo = i <= j ? i : j;
      final hi = i <= j ? j : i;
      final String k;
      if (curve is BrepLine) {
        k = '$lo-$hi-line';
      } else {
        final mid = curve.evaluate(0.5 * (curve.t0 + curve.t1));
        final q =
            '${mid.x.toStringAsFixed(5)}_${mid.y.toStringAsFixed(5)}_${mid.z.toStringAsFixed(5)}';
        k = '$lo-$hi-$q';
      }
      final existing = edgeOf[k];
      if (existing != null) {
        final e = b.edges[existing];
        return (existing, e.v0 != i);
      }
      final id = b.edge(i, j, curve);
      edgeOf[k] = id;
      return (id, false);
    }

    for (final face in s.faces) {
      final loops = <BrepLoop>[];
      for (final loop in face.loops) {
        final coedges = <BrepCoedge>[];
        for (final c in loop.coedges) {
          final e = s.edges[c.edge];
          final start = s.vertices[c.reversed ? e.v1 : e.v0];
          final end = s.vertices[c.reversed ? e.v0 : e.v1];
          if (e.curve is! BrepLine) {
            final i0 = weld(s.vertices[e.v0]);
            final i1 = weld(s.vertices[e.v1]);
            final closedCurve =
                i0 == i1 &&
                e.curve is BrepCircle &&
                (e.curve as BrepCircle).isClosed;
            if (i0 == i1 && !closedCurve) continue;
            final sewn = sew(i0, i1, e.curve);
            coedges.add(BrepCoedge(sewn.$1, reversed: sewn.$2 != c.reversed));
            continue;
          }
          final mids = <(double, Vector3)>[];
          for (final v in s.vertices) {
            final t = _segmentInteriorParam(v, start, end);
            if (t == null) continue;
            mids.add((t, v));
          }
          mids.sort((a, b) => a.$1.compareTo(b.$1));
          var prev = start;
          for (final mid in mids) {
            if ((mid.$2 - prev).length < 1e-6) continue;
            if ((mid.$2 - end).length < 1e-6) continue;
            final sewn = sew(weld(prev), weld(mid.$2), BrepLine(prev, mid.$2));
            coedges.add(BrepCoedge(sewn.$1, reversed: sewn.$2));
            prev = mid.$2;
          }
          if ((prev - end).length < 1e-6) continue;
          final sewn = sew(weld(prev), weld(end), BrepLine(prev, end));
          coedges.add(BrepCoedge(sewn.$1, reversed: sewn.$2));
        }
        if (coedges.length >= 3 ||
            coedges.isNotEmpty &&
                loop.coedges.any((c) => s.edges[c.edge].curve is! BrepLine)) {
          loops.add(BrepLoop(coedges));
        }
      }
      if (loops.isEmpty && face.loops.isNotEmpty) continue;
      b.face(face.surface, loops, reversed: face.reversed);
    }
    if (b.faces.isEmpty) return s;
    return b.build();
  }

  static bool _pointOnSegmentInterior(Vector3 p, Vector3 a, Vector3 b) =>
      _segmentInteriorParam(p, a, b) != null;

  /// Parameter in (0, 1) when [p] lies on the open segment [a]–[b].
  static double? _segmentInteriorParam(Vector3 p, Vector3 a, Vector3 b) {
    final ab = b - a;
    final len2 = ab.length2;
    if (len2 < 1e-18) return null;
    final t = (p - a).dot(ab) / len2;
    if (t <= 1e-8 || t >= 1 - 1e-8) return null;
    final dist2 = (p - (a + ab * t)).length2;
    if (dist2 > 1e-12) return null;
    return t;
  }

  static List<_HealSeg> _collapseSegs(List<_HealSeg> segs) {
    if (segs.length < 4) return segs;
    var cur = List<_HealSeg>.from(segs);
    for (var pass = 0; pass < 8; pass++) {
      if (cur.length < 4) break;
      final out = <_HealSeg>[];
      var i = 0;
      var mergedAny = false;
      while (i < cur.length) {
        final a = cur[i];
        final b = cur[(i + 1) % cur.length];
        if (i + 1 < cur.length && _collinearLines(a, b)) {
          out.add(_HealSeg(a.a, b.b, BrepLine(a.a, b.b)));
          i += 2;
          mergedAny = true;
        } else {
          out.add(a);
          i++;
        }
      }
      if (!mergedAny) {
        final first = out.first;
        final last = out.last;
        if (out.length >= 4 && _collinearLines(last, first)) {
          out
            ..removeLast()
            ..removeAt(0);
          out.insert(0, _HealSeg(last.a, first.b, BrepLine(last.a, first.b)));
          mergedAny = true;
        }
      }
      if (!mergedAny || out.length < 3) break;
      cur = out;
    }
    return cur.length >= 3 ? cur : segs;
  }

  static bool _collinearLines(_HealSeg a, _HealSeg b) {
    if (a.curve is! BrepLine || b.curve is! BrepLine) return false;
    if ((a.b - b.a).length2 > 1e-16) return false;
    final e1 = a.b - a.a;
    final e2 = b.b - b.a;
    final l1 = e1.length;
    final l2 = e2.length;
    if (l1 < 1e-9 || l2 < 1e-9) return true;
    final cr = e1.cross(e2).length / (l1 * l2);
    final dot = e1.dot(e2) / (l1 * l2);
    return cr <= 1e-6 && dot >= 0.999999;
  }

  static bool _edgeUsedByNonPlane(List<BrepFace> faces, int eid) {
    for (final f in faces) {
      if (f.surface is BrepPlane) continue;
      for (final loop in f.loops) {
        if (loop.coedges.any((c) => c.edge == eid)) return true;
      }
    }
    return false;
  }

  static int? _sharedEdge(BrepFace a, BrepFace b) {
    final ea = {for (final c in a.loops.first.coedges) c.edge};
    for (final c in b.loops.first.coedges) {
      if (ea.contains(c.edge)) return c.edge;
    }
    return null;
  }

  static Set<int> _sharedEdges(BrepFace a, BrepFace b) {
    final ea = {for (final c in a.loops.first.coedges) c.edge};
    return {
      for (final c in b.loops.first.coedges)
        if (ea.contains(c.edge)) c.edge,
    };
  }

  /// Start index and length of the one contiguous cyclic run of [shared]
  /// coedges in [loop], or null when they do not form a single run.
  static (int, int)? _sharedRun(List<BrepCoedge> loop, Set<int> shared) {
    final n = loop.length;
    if (shared.length >= n) return null;
    var start = -1;
    for (var i = 0; i < n; i++) {
      final prev = loop[(i - 1 + n) % n];
      if (shared.contains(loop[i].edge) && !shared.contains(prev.edge)) {
        if (start >= 0) return null;
        start = i;
      }
    }
    if (start < 0) return null;
    var len = 0;
    while (len < n && shared.contains(loop[(start + len) % n].edge)) {
      len++;
    }
    return len == shared.length ? (start, len) : null;
  }

  /// Join two coplanar faces along a contiguous run of shared edges.
  ///
  /// Splicing out only one of several shared edges would leave the rest in
  /// the merged loop twice, which is a pinched, non-manifold face.
  static BrepLoop? _mergeSharedRun(BrepFace a, BrepFace b, Set<int> shared) {
    final la = a.loops.first.coedges;
    var lb = b.loops.first.coedges;
    // Face normals and trim winding are independent. A subtract keeps the
    // cutter's wire but reverses its face normal, so equal-facing patches can
    // traverse their common edge in the same direction. Reverse the whole
    // second wire before splicing, otherwise the result jumps between vertices.
    final first = shared.first;
    if (la.firstWhere((c) => c.edge == first).reversed ==
        lb.firstWhere((c) => c.edge == first).reversed) {
      lb = [
        for (final c in lb.reversed) BrepCoedge(c.edge, reversed: !c.reversed),
      ];
    }
    if (shared.any(
      (e) =>
          la.firstWhere((c) => c.edge == e).reversed ==
          lb.firstWhere((c) => c.edge == e).reversed,
    )) {
      return null;
    }
    final ra = _sharedRun(la, shared);
    final rb = _sharedRun(lb, shared);
    if (ra == null || rb == null) return null;

    List<BrepCoedge> rest(List<BrepCoedge> loop, (int, int) run) {
      final n = loop.length;
      return [
        for (var k = 0; k < n - run.$2; k++) loop[(run.$1 + run.$2 + k) % n],
      ];
    }

    final out = [...rest(la, ra), ...rest(lb, rb)];
    if (out.length < 3) return null;
    return BrepLoop(out);
  }

  static BrepLoop? _mergeLoops(BrepFace a, BrepFace b, int shared) {
    final la = a.loops.first.coedges;
    var lb = b.loops.first.coedges;
    final ia = la.indexWhere((c) => c.edge == shared);
    var ib = lb.indexWhere((c) => c.edge == shared);
    if (ia < 0 || ib < 0) return null;
    if (la[ia].reversed == lb[ib].reversed) {
      lb = [
        for (final c in lb.reversed) BrepCoedge(c.edge, reversed: !c.reversed),
      ];
      ib = lb.indexWhere((c) => c.edge == shared);
    }
    final out = <BrepCoedge>[
      ...la.sublist(ia + 1),
      ...la.sublist(0, ia),
      ...lb.sublist(ib + 1),
      ...lb.sublist(0, ib),
    ];
    if (out.length < 3) return null;
    return BrepLoop(out);
  }

  static BrepSolid _trianglesToBrep(MeshData mesh) {
    final b = BrepBuilder();
    int weld(Vector3 p) {
      for (var i = 0; i < b.vertices.length; i++) {
        if ((b.vertices[i] - p).length2 < 1e-16) return i;
      }
      return b.vertex(p);
    }

    final edgeOf = <String, int>{};
    int line(int i, int j) {
      final lo = i < j ? i : j;
      final hi = i < j ? j : i;
      final k = '$lo-$hi';
      final existing = edgeOf[k];
      if (existing != null) return existing;
      final id = b.edge(i, j, BrepLine(b.vertices[i], b.vertices[j]));
      edgeOf[k] = id;
      return id;
    }

    for (final face in mesh.faces) {
      if (face.length < 3) continue;
      final pts = [for (final i in face) mesh.vertices[i]];
      var n = Vector3.zero();
      for (var i = 0; i < pts.length; i++) {
        n += pts[i].cross(pts[(i + 1) % pts.length]);
      }
      if (n.length2 < 1e-20) continue;
      n.normalize();
      final ids = [for (final p in pts) weld(p)];
      final unique = <int>[];
      for (final id in ids) {
        if (unique.isEmpty || unique.last != id) unique.add(id);
      }
      if (unique.length >= 2 && unique.first == unique.last)
        unique.removeLast();
      if (unique.length < 3) continue;
      final eids = <int>[];
      for (var i = 0; i < unique.length; i++) {
        final a = unique[i];
        final c = unique[(i + 1) % unique.length];
        if (a == c) continue;
        eids.add(line(a, c));
      }
      if (eids.length < 3) continue;
      b.face(BrepPlane(origin: pts.first, normal: n), [
        BrepLoop([
          for (var i = 0; i < eids.length; i++)
            BrepCoedge(eids[i], reversed: b.edges[eids[i]].v0 != unique[i]),
        ]),
      ]);
    }
    if (b.faces.isEmpty) {
      throw BrepUnsupported('Mesh has no triangles to convert to B-Rep');
    }
    return b.build();
  }

  /// Every used edge has exactly two coedges. Unused edges are ignored.
  ///
  /// A complete circular boundary has no distinguished starting point. Sew
  /// coincident full rims by their geometry, including opposite circle frames.
  /// This also removes redundant chart seams from complete revolutions.
  static BrepSolid sewClosedRims(BrepSolid s, {bool removeSeams = false}) {
    final rims = <int>[];
    final replacements = <int, (int, bool)>{};
    for (var i = 0; i < s.edges.length; i++) {
      final c = s.edges[i].curve;
      if (c is! BrepCircle || !c.isClosed || c.radius < 1e-7) continue;
      for (final id in rims) {
        final other = s.edges[id].curve as BrepCircle;
        if ((c.center - other.center).length > 1e-7 ||
            (c.radius - other.radius).abs() > 1e-7 ||
            c.normal.cross(other.normal).length > 1e-9)
          continue;
        replacements[i] = (id, c.normal.dot(other.normal) < 0);
        break;
      }
      if (!replacements.containsKey(i)) rims.add(i);
    }
    final faces = <BrepFace>[];
    var changed = replacements.isNotEmpty;
    for (final face in s.faces) {
      var loops = <BrepLoop>[];
      for (final sourceLoop in face.loops) {
        var coedges = sourceLoop.coedges.where((co) {
          final curve = s.edges[co.edge].curve;
          return curve is! BrepLine || (curve.b - curve.a).length > 1e-9;
        }).toList();
        if (coedges.length != sourceLoop.coedges.length) changed = true;
        var cancelled = true;
        while (cancelled && coedges.length >= 2) {
          cancelled = false;
          for (var i = 0; i < coedges.length; i++) {
            final j = (i + 1) % coedges.length, a = coedges[i], b = coedges[j];
            if (a.edge == b.edge && a.reversed != b.reversed) {
              coedges = [
                for (var k = 0; k < coedges.length; k++)
                  if (k != i && k != j) coedges[k],
              ];
              changed = true;
              cancelled = true;
              break;
            }
          }
        }
        if (coedges.isEmpty && sourceLoop.coedges.isNotEmpty) continue;
        final loop = BrepLoop(coedges);
        final full = loop.coedges.where((co) {
          final c = s.edges[co.edge].curve;
          return c is BrepCircle && c.isClosed;
        }).toList();
        final seams = loop.coedges.where((co) => !full.contains(co)).toList();
        if ((face.surface is BrepCylinder || face.surface is BrepCone) &&
            full.isNotEmpty &&
            seams.length == 2 &&
            seams[0].edge == seams[1].edge &&
            seams[0].reversed != seams[1].reversed &&
            (removeSeams ||
                full.any((c) => replacements.containsKey(c.edge)))) {
          loops.addAll([
            for (final c in full) BrepLoop([c]),
          ]);
          changed = true;
        } else {
          loops.add(loop);
        }
      }
      if (loops.isEmpty && face.loops.isNotEmpty && face.surface is BrepPlane)
        continue;
      if (face.surface is BrepPlane) {
        final circles = [
          for (final loop in loops)
            if (loop.coedges.length == 1 &&
                s.edges[loop.coedges.single.edge].curve is BrepCircle)
              s.edges[loop.coedges.single.edge].curve as BrepCircle,
        ];
        loops = loops.where((loop) {
          if (loop.coedges.length < 8 ||
              !loop.coedges.every((co) => s.edges[co.edge].curve is BrepLine))
            return true;
          final points = [
            for (final co in loop.coedges) s.vertices[s.edges[co.edge].v0],
          ];
          for (final circle in circles) {
            if (circle.isClosed &&
                points.every(
                  (p) => (p - circle.closestPoint(p)).length < 1e-7,
                )) {
              changed = true;
              return false;
            }
          }
          return true;
        }).toList();
      }
      faces.add(
        BrepFace(
          surface: face.surface,
          reversed: face.reversed,
          origin: face.origin,
          loops: [
            for (final loop in loops)
              BrepLoop([
                for (final co in loop.coedges)
                  if (replacements[co.edge] case final replacement?)
                    BrepCoedge(
                      replacement.$1,
                      reversed: co.reversed != replacement.$2,
                    )
                  else
                    co,
              ]),
          ],
        ),
      );
    }
    if (!changed) return s;
    return BrepSolid(vertices: s.vertices, edges: s.edges, faces: faces);
  }

  static bool topologyIsClosed(BrepSolid s) {
    if (s.faces.isEmpty) return false;
    if (s.vertices.any((p) => !p.x.isFinite || !p.y.isFinite || !p.z.isFinite))
      return false;
    bool closedSurface(BrepSurface surface) =>
        surface is BrepSplineSurface && surface.isClosed ||
        surface is BrepSphere ||
        surface is BrepTorus ||
        (surface is BrepAffineSurface && closedSurface(surface.base));
    if (s.edges.isEmpty)
      return s.faces.every((f) => f.loops.isEmpty && closedSurface(f.surface));
    final uses = List<int>.filled(s.edges.length, 0);
    for (final f in s.faces) {
      if (f.loops.isEmpty && !closedSurface(f.surface)) return false;
      for (final loop in f.loops) {
        if (loop.coedges.isEmpty) return false;
        for (final c in loop.coedges) {
          if (c.edge < 0 || c.edge >= uses.length) return false;
          final e = s.edges[c.edge];
          if (e.v0 < 0 ||
              e.v1 < 0 ||
              e.v0 >= s.vertices.length ||
              e.v1 >= s.vertices.length)
            return false;
          uses[c.edge]++;
        }
      }
    }
    for (var i = 0; i < uses.length; i++) {
      final u = uses[i], curve = s.edges[i].curve;
      // A zero-radius rim is a collapsed parameter boundary at a cone apex.
      if (curve is BrepCircle && curve.isClosed && curve.radius <= 1e-7)
        continue;
      if (u == 1 || u > 2) return false;
    }
    return true;
  }

  static String brokenSolidMessage(String op) =>
      '$op produced a broken solid and was not applied.';

  /// Tessellate a closed solid, including explicitly paired line contacts.
  /// Throws instead of returning a hollow / open display mesh.
  static MeshData tessellateClosed(
    BrepSolid s, {
    MeshTessellation tessellation = MeshTessellation.standard,
    String op = 'This change',
  }) {
    if (!topologyIsClosed(s)) throw BrepUnsupported(brokenSolidMessage(op));
    // A closed edge graph is not sufficient when a face has been moved off
    // its trim. The planar mesher uses the shared edge positions, which could
    // otherwise conceal a missing or misplaced cap behind a closed display.
    for (final face in s.faces) {
      final plane = face.surface;
      if (plane is! BrepPlane) continue;
      for (final loop in face.loops) {
        for (final co in loop.coedges) {
          final curve = s.edges[co.edge].curve;
          for (final fraction in [0.0, .25, .5, .75, 1.0]) {
            final p = curve.evaluate(
              curve.t0 + (curve.t1 - curve.t0) * fraction,
            );
            if ((p - plane.origin).dot(plane.normal).abs() > 1e-5) {
              throw BrepUnsupported(brokenSolidMessage(op));
            }
          }
        }
      }
    }
    final mesh = BrepTessellate.run(s, tessellation: tessellation);
    if (s.faces.isEmpty || mesh.faces.isEmpty) {
      throw BrepUnsupported(brokenSolidMessage(op));
    }
    if (MeshSeal.isClosedSolid(mesh) || _closedPlanarContacts(s, mesh))
      return mesh;
    final topo = topologyIsClosed(s);
    IjoLog.line(
      'tessellateClosed open=${MeshSeal.openBoundaryEdgeCount(mesh)} '
      'topo=$topo faces=${s.faces.length}',
    );
    if (topo) {
      // Repair sampling seams only. Capping a missing BREP surface invents
      // geometry and can conceal an invalid offset or inverted wall.
      final sealed = MeshSeal.preserveClosed(mesh);
      if (MeshSeal.isClosedSolid(sealed) || _closedPlanarContacts(s, sealed)) {
        return sealed;
      }
    }
    throw BrepUnsupported(brokenSolidMessage(op));
  }

  /// A position-welded display can have four triangles on a contact line,
  /// even though the BREP has two distinct, closed edges there. Accept only
  /// exact straight contacts between planar or ruled generator walls whose
  /// oriented triangle pairs match those BREP edges. Gaps and missing faces
  /// remain invalid.
  static bool _closedPlanarContacts(BrepSolid s, MeshData mesh) {
    if (!hasPairedPlanarContacts(s)) return false;
    final faceIds = mesh.faceIds;
    if (faceIds == null || faceIds.length != mesh.faces.length) return false;
    final positions = <(int, int, int), int>{};
    final points = <Vector3>[];
    final ids = <int>[];
    for (final p in mesh.vertices) {
      if (!p.x.isFinite || !p.y.isFinite || !p.z.isFinite) return false;
      ids.add(
        positions.putIfAbsent(
          ((p.x * 1e6).round(), (p.y * 1e6).round(), (p.z * 1e6).round()),
          () {
            points.add(p);
            return points.length - 1;
          },
        ),
      );
    }
    final uses = <(int, int), List<(int, bool)>>{};
    for (var i = 0; i < mesh.faces.length; i++) {
      final f = mesh.faces[i], fi = faceIds[i];
      if (fi < 0 ||
          fi >= s.faces.length ||
          f.length != 3 ||
          f.any((v) => v < 0 || v >= ids.length))
        return false;
      final ring = f.map((v) => ids[v]).toList();
      if (ring.toSet().length < 3) continue;
      for (var j = 0; j < 3; j++) {
        final a = ring[j], b = ring[(j + 1) % 3];
        (uses[a < b ? (a, b) : (b, a)] ??= []).add((fi, a < b));
      }
    }
    final crowded = <(int, int), List<(int, bool)>>{};
    for (final entry in uses.entries) {
      final fs = entry.value;
      if (fs.length == 2 && fs[0].$2 != fs[1].$2) continue;
      if (fs.length < 4 ||
          fs.length.isOdd ||
          fs.map((f) => f.$1).toSet().length != fs.length)
        return false;
      crowded[entry.key] = fs;
    }
    if (crowded.isEmpty) return false;
    final edgeFaces = <int, Set<int>>{};
    for (var fi = 0; fi < s.faces.length; fi++) {
      for (final l in s.faces[fi].loops) {
        for (final c in l.coedges) (edgeFaces[c.edge] ??= {}).add(fi);
      }
    }
    for (final entry in crowded.entries) {
      final a = points[entry.key.$1], b = points[entry.key.$2];
      final fs = {for (final use in entry.value) use.$1: use.$2};
      final faceEdges = <int, int>{};
      for (final e in edgeFaces.entries) {
        if (e.value.length != 2 || !e.value.every(fs.containsKey)) continue;
        final curve = s.edges[e.key].curve;
        if (curve is! BrepLine) continue;
        final d = curve.b - curve.a;
        if (d.length2 < 1e-16) continue;
        bool onSegment(Vector3 p) {
          final t = (p - curve.a).dot(d) / d.length2;
          return t >= -1e-8 &&
              t <= 1 + 1e-8 &&
              (p - (curve.a + d * t)).length < 1e-6;
        }

        if (!onSegment(a) || !onSegment(b)) continue;
        for (final fi in e.value) {
          if (faceEdges.containsKey(fi)) return false;
          faceEdges[fi] = e.key;
        }
      }
      if (faceEdges.length != fs.length) return false;
      final pairs = planarContactFacePairs(s, faceEdges);
      if (pairs == null) return false;
      for (final pair in pairs) {
        if (faceEdges[pair[0]] != faceEdges[pair[1]] ||
            fs[pair[0]] == fs[pair[1]])
          return false;
      }
    }
    return true;
  }

  static BrepSolid _orientOutward(BrepSolid s) {
    final mesh = BrepTessellate.run(s);
    if (mesh.faces.isEmpty) return s;
    if (MeshEdit.signedVolume(mesh) >= 0) return s;
    return BrepSolid(
      vertices: s.vertices,
      edges: s.edges,
      faces: [
        for (final f in s.faces)
          BrepFace(
            surface: f.surface,
            loops: f.loops,
            reversed: !f.reversed,
            origin: f.origin,
          ),
      ],
    );
  }
}

class _HealSeg {
  _HealSeg(this.a, this.b, this.curve, {this.reversed = false});
  final bool reversed;
  final Vector3 a;
  final Vector3 b;
  final BrepCurve curve;
}
