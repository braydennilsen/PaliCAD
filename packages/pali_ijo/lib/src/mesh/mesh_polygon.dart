import 'dart:math' as math;

import 'package:vector_math/vector_math_64.dart';

/// Planar polygon triangulation, including holes.
///
/// Triangulate the 2-D profile first, then
/// sweep. Interior diagonals are shared by two triangles; only the outer and
/// hole boundaries are 1-use edges (those become the walls).
class MeshPolygon {
  MeshPolygon._();

  /// [rings].first is the outer loop; the rest are holes. Indices into [verts].
  ///
  /// When [keepCollinear] is true, vertices that lie on a straight edge stay
  /// in the triangulation (plane-cut T-junctions). The default drops them so
  /// a sampled line in a sketch profile does not become a sliver fan.
  static List<List<int>> triangulateRings(
    List<Vector3> verts,
    List<List<int>> rings,
    Vector3 normal, {
    bool keepCollinear = false,
  }) {
    if (rings.isEmpty) return const [];
    final nHat = normal.length2 < 1e-20
        ? Vector3(0, 0, 1)
        : normal.normalized();
    var u = nHat.cross(Vector3(1, 0, 0));
    if (u.length2 < 1e-12) u = nHat.cross(Vector3(0, 1, 0));
    u.normalize();
    final v = nHat.cross(u).normalized();
    final origin = verts[rings.first.first];

    final projected = [
      for (final p in verts) Vector2((p - origin).dot(u), (p - origin).dot(v)),
    ];
    Vector2 uv(int i) => projected[i];

    List<int> clean(List<int> ring) {
      final out = <int>[];
      for (final i in ring) {
        if (i < 0 || i >= verts.length) continue;
        if (out.isEmpty || out.last != i) out.add(i);
      }
      if (out.length > 1 && out.first == out.last) out.removeLast();
      return out;
    }

    double areaOf(List<int> ring) {
      var a = 0.0;
      for (var i = 0; i < ring.length; i++) {
        final p = uv(ring[i]);
        final q = uv(ring[(i + 1) % ring.length]);
        a += p.x * q.y - q.x * p.y;
      }
      return a;
    }

    final cleaned = [for (final r in rings) clean(r)];
    if (cleaned.first.length < 3) return const [];

    var outer = List<int>.from(cleaned.first);
    if (areaOf(outer) < 0) outer = outer.reversed.toList();

    final holes = <List<int>>[
      for (var i = 1; i < cleaned.length; i++)
        if (cleaned[i].length >= 3) List<int>.from(cleaned[i]),
    ];
    for (var i = 0; i < holes.length; i++) {
      if (areaOf(holes[i]) > 0) holes[i] = holes[i].reversed.toList();
    }
    holes.sort((a, b) {
      var ax = -1e300;
      var bx = -1e300;
      for (final i in a) {
        ax = math.max(ax, uv(i).x);
      }
      for (final i in b) {
        bx = math.max(bx, uv(i).x);
      }
      return bx.compareTo(ax);
    });

    for (var i = 0; i < holes.length; i++) {
      final merged = _bridge(outer, holes[i], uv, holes.skip(i + 1));
      if (merged == null) return const [];
      outer = merged;
    }
    return _earClip(verts, outer, nHat, uv, keepCollinear: keepCollinear);
  }

  /// Convenience: 3-D rings on one plane → triangles as point triples.
  static List<List<Vector3>> triangulate(
    List<Vector3> outer,
    List<List<Vector3>> holes,
    Vector3 normal,
  ) {
    final verts = <Vector3>[...outer];
    final rings = <List<int>>[
      [for (var i = 0; i < outer.length; i++) i],
    ];
    for (final hole in holes) {
      if (hole.length < 3) continue;
      final start = verts.length;
      verts.addAll(hole);
      rings.add([for (var i = 0; i < hole.length; i++) start + i]);
    }
    return [
      for (final t in triangulateRings(verts, rings, normal))
        [verts[t[0]], verts[t[1]], verts[t[2]]],
    ];
  }

  // Rotation/projection of metre-scale profiles leaves roundoff in
  // collinear turns. Use the same scale-aware test for hole bridges and
  // ears, so bridges do not run along an existing hole's boundary.
  static double _turnTolerance(Vector2 a, Vector2 b, Vector2 c) =>
      1e-12 + 4e-15 * (b - a).length * (c - a).length;

  static List<int>? _bridge(
    List<int> outer,
    List<int> hole,
    Vector2 Function(int) uv,
    Iterable<List<int>> pendingHoles,
  ) {
    // Openings can meet at a vertex. Join their boundary walks there; adding
    // a second bridge leaves a zero-area cycle that cannot be ear-clipped.
    for (var oi = 0; oi < outer.length; oi++) {
      for (var hi = 0; hi < hole.length; hi++) {
        if ((uv(outer[oi]) - uv(hole[hi])).length2 > 1e-20) continue;
        return [
          ...outer.take(oi + 1),
          for (var k = 1; k <= hole.length; k++) hole[(hi + k) % hole.length],
          ...outer.skip(oi + 1),
        ];
      }
    }
    var hi = 0;
    var hx = uv(hole[0]).x;
    for (var i = 1; i < hole.length; i++) {
      final x = uv(hole[i]).x;
      if (x > hx) {
        hx = x;
        hi = i;
      }
    }
    final h = hole[hi];
    final hp = uv(h);

    // A bridge must enter the material sector at both ends. Merely testing
    // proper segment crossings misses a bridge through a reflex vertex, or
    // one attached to the wrong occurrence of an earlier bridge junction.
    // Both the CCW outer wire and CW hole wire have material on their left.
    bool entersMaterial(List<int> ring, int index, Vector2 point) {
      final p = uv(ring[(index + ring.length - 1) % ring.length]);
      final q = uv(ring[index]);
      final r = uv(ring[(index + 1) % ring.length]);
      double cross(Vector2 a, Vector2 b, Vector2 c) =>
          (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
      final after = cross(q, r, point) > _turnTolerance(q, r, point);
      final before = cross(q, point, p) > _turnTolerance(q, point, p);
      return cross(p, q, r) >= -_turnTolerance(p, q, r)
          ? after && before
          : after || before;
    }

    bool crosses(Vector2 a, Vector2 b, Vector2 c, Vector2 d) {
      final r = b - a;
      final s = d - c;
      final den = r.x * s.y - r.y * s.x;
      if (den.abs() < 1e-14) {
        // A bridge along an earlier hole edge creates an overlapping wire,
        // even though there is no proper crossing. Aligned lettering and
        // repeated holes commonly have exactly these collinear candidates.
        if ((r.x * (c.y - a.y) - r.y * (c.x - a.x)).abs() >
            _turnTolerance(a, b, c))
          return false;
        final length2 = r.length2;
        if (length2 < 1e-24) return true;
        final t0 = (c - a).dot(r) / length2;
        final t1 = (d - a).dot(r) / length2;
        return math.max(0, math.min(t0, t1)) <
            math.min(1, math.max(t0, t1)) - 1e-10;
      }
      final t = ((c.x - a.x) * s.y - (c.y - a.y) * s.x) / den;
      final u = ((c.x - a.x) * r.y - (c.y - a.y) * r.x) / den;
      // Touching a third vertex in the interior of the bridge also splits
      // the polygon unless explicitly represented as a junction.
      return t > 1e-8 && t < 1 - 1e-8 && u >= -1e-10 && u <= 1 + 1e-10;
    }

    bool blocked(int oi) {
      final a = hp;
      final b = uv(outer[oi]);
      bool hitRing(List<int> ring) {
        for (var i = 0; i < ring.length; i++) {
          final c = uv(ring[i]);
          final d = uv(ring[(i + 1) % ring.length]);
          if (crosses(a, b, c, d)) return true;
        }
        return false;
      }

      return hitRing(outer) || hitRing(hole) || pendingHoles.any(hitRing);
    }

    var best = -1;
    var bestScore = double.infinity;
    for (var i = 0; i < outer.length; i++) {
      if (!entersMaterial(outer, i, hp) ||
          !entersMaterial(hole, hi, uv(outer[i]))) {
        continue;
      }
      final d = (uv(outer[i]) - hp).length2;
      final preferRight = uv(outer[i]).x + 1e-9 < hx ? 1e6 : 0.0;
      final score = d + preferRight;
      if (score >= bestScore || blocked(i)) continue;
      if (score < bestScore) {
        bestScore = score;
        best = i;
      }
    }
    if (best < 0) return null;

    final merged = <int>[
      ...outer.sublist(0, best + 1),
      for (var k = 0; k <= hole.length; k++) hole[(hi + k) % hole.length],
      outer[best],
      ...outer.sublist(best + 1),
    ];
    return merged;
  }

  static List<List<int>> _earClip(
    List<Vector3> verts,
    List<int> loop,
    Vector3 normal,
    Vector2 Function(int) uv, {
    bool keepCollinear = false,
  }) {
    if (loop.length < 3) return const [];
    if (loop.length == 3) {
      final n = (verts[loop[1]] - verts[loop[0]]).cross(
        verts[loop[2]] - verts[loop[0]],
      );
      return n.dot(normal) < 0
          ? [
              [loop[0], loop[2], loop[1]],
            ]
          : [List<int>.from(loop)];
    }

    double cross2(Vector2 a, Vector2 b, Vector2 c) =>
        (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);

    bool inside(Vector2 p, Vector2 a, Vector2 b, Vector2 c) {
      final c1 = cross2(a, b, p);
      final c2 = cross2(b, c, p);
      final c3 = cross2(c, a, p);
      final e1 = math.max(1e-12, _turnTolerance(a, b, p));
      final e2 = math.max(1e-12, _turnTolerance(b, c, p));
      final e3 = math.max(1e-12, _turnTolerance(c, a, p));
      final hasNeg = c1 < -e1 || c2 < -e2 || c3 < -e3;
      final hasPos = c1 > e1 || c2 > e2 || c3 > e3;
      return !(hasNeg && hasPos);
    }

    final active = [for (var i = 0; i < loop.length; i++) i];
    final tris = <List<int>>[];
    var guard = 0;
    final maxGuard = loop.length * loop.length + 8;
    while (active.length > 3 && guard < maxGuard) {
      guard++;
      var clipped = false;
      for (var i = 0; i < active.length; i++) {
        final i0 = active[(i + active.length - 1) % active.length];
        final i1 = active[i];
        final i2 = active[(i + 1) % active.length];
        final a = uv(loop[i0]);
        final b = uv(loop[i1]);
        final c = uv(loop[i2]);
        final cr = cross2(a, b, c);
        final epsilon = _turnTolerance(a, b, c);
        if (cr.abs() <= epsilon) {
          if (keepCollinear) continue;
          active.removeAt(i);
          clipped = true;
          break;
        }
        if (cr <= epsilon) continue;
        var blocked = false;
        for (final j in active) {
          if (j == i0 || j == i1 || j == i2) continue;
          final p = uv(loop[j]);
          // A hole bridge reuses the same vertex twice in the loop. That
          // second copy sits on the ear and is not an interior blocker.
          if ((p - a).length2 < 1e-20 ||
              (p - b).length2 < 1e-20 ||
              (p - c).length2 < 1e-20) {
            continue;
          }
          if (inside(p, a, b, c)) {
            blocked = true;
            break;
          }
        }
        if (blocked) continue;
        if (active.length == 4) {
          final remaining = active.where((j) => j != i1).toList();
          final x = uv(loop[remaining[0]]),
              y = uv(loop[remaining[1]]),
              z = uv(loop[remaining[2]]);
          // Do not leave a collinear final triangle: that loses the middle
          // boundary vertex when zero-area triangles are dropped downstream.
          if (cross2(x, y, z) <= _turnTolerance(x, y, z)) continue;
        }
        tris.add([loop[i0], loop[i1], loop[i2]]);
        active.removeAt(i);
        clipped = true;
        break;
      }
      if (!clipped) break;
    }
    if (active.length != 3) return const [];
    tris.add([loop[active[0]], loop[active[1]], loop[active[2]]]);
    if (tris.isEmpty) return const [];
    final t0 = tris.first;
    final tn = (verts[t0[1]] - verts[t0[0]]).cross(verts[t0[2]] - verts[t0[0]]);
    if (tn.dot(normal) < 0) {
      return [for (final t in tris) t.reversed.toList()];
    }
    return tris;
  }
}
