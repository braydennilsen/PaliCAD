import 'dart:math' as math;
import 'package:vector_math/vector_math_64.dart';
import 'geom.dart';
import 'solid.dart';
import 'chamfer_options.dart';
import 'classify.dart';

/// Blends a closed extrusion rim. Smooth curved joins and sharp line-line
/// miters share contact curves. Contacts follow the source surfaces;
/// cubic spans approximate the longitudinal offset to [tolerance], while each
/// fillet cross-section is a rational circular arc. No display mesh is used.
class BrepProfileBlend {
  static const tolerance = 1e-5;

  /// Null means the selection is not a supported extrusion rim. Invalid offsets
  /// on a recognized rim throw before any source topology is changed.
  static BrepSolid? tryRun(
    BrepSolid s,
    List<int> selected,
    double distance, {
    required bool round,
    ChamferOptions options = const ChamferOptions(),
  }) {
    if (selected.isEmpty) return null;
    // Older edge blends sometimes stored a reused line backwards under its
    // topological endpoints. Reorient only exact endpoint swaps, on a copy;
    // this repairs parameter direction without changing the geometric edge.
    final oriented = <int, BrepEdge>{};
    for (var i = 0; i < s.edges.length; i++) {
      final e = s.edges[i], c = e.curve;
      if (c is BrepLine &&
          e.v0 != e.v1 &&
          c.a.distanceTo(s.vertices[e.v1]) < 1e-7 &&
          c.b.distanceTo(s.vertices[e.v0]) < 1e-7 &&
          c.a.distanceTo(s.vertices[e.v0]) > 1e-7) {
        oriented[i] = BrepEdge(
          v0: e.v0,
          v1: e.v1,
          curve: BrepLine(c.b, c.a),
          origin: e.origin,
        );
      }
    }
    if (oriented.isNotEmpty) {
      s = BrepSolid(
        vertices: s.vertices,
        faces: s.faces,
        edges: [
          for (var i = 0; i < s.edges.length; i++) oriented[i] ?? s.edges[i],
        ],
      );
    }
    for (var capId = 0; capId < s.faces.length; capId++) {
      final cap = s.faces[capId];
      if (cap.surface is! BrepPlane) continue;
      for (final loop in cap.loops) {
        final ids = loop.coedges.map((c) => c.edge).toSet();
        if (!selected.every(ids.contains)) continue;
        if (ids.length == 1 &&
            s.edges[ids.first].curve is BrepCircle &&
            s.faces.any(
              (f) =>
                  f.surface is BrepCylinder &&
                  f.loops.any((l) => l.coedges.any((c) => c.edge == ids.first)),
            ))
          continue;
        // Tangent propagation is permitted along a smooth closed rim. Sharp
        // line-line joins require the whole rim to be explicitly selected.
        if (!ids.any((id) => s.edges[id].curve is! BrepLine)) {
          // Convex caps retain analytic corner spheres. Concave/annular rims
          // need a coordinated wire offset rather than independent corners.
          final pts = [
            for (final co in loop.coedges)
              s.vertices[co.reversed
                  ? s.edges[co.edge].v1
                  : s.edges[co.edge].v0],
          ];
          final n = cap.outwardNormal(0, 0);
          var positive = false, negative = false;
          for (var i = 0; i < pts.length; i++) {
            final turn = (pts[(i + 1) % pts.length] - pts[i])
                .cross(pts[(i + 2) % pts.length] - pts[(i + 1) % pts.length])
                .dot(n);
            positive |= turn > 1e-8;
            negative |= turn < -1e-8;
          }
          if (!(positive && negative) && cap.loops.length == 1) continue;
          if (!ids.every(selected.contains)) continue;
        }
        final walls = <int, int>{};
        var eligible = true;
        for (final id in ids) {
          final adj = <int>[
            for (var fi = 0; fi < s.faces.length; fi++)
              if (fi != capId &&
                  s.faces[fi].loops.any(
                    (l) => l.coedges.any((c) => c.edge == id),
                  ))
                fi,
          ];
          if (adj.length != 1) {
            eligible = false;
            break;
          }
          walls[id] = adj.single;
        }
        if (!eligible) continue;
        final n = cap.outwardNormal(0, 0);
        final depth = round
            ? distance
            : options.secondDistance(distance, math.pi / 2);
        final stations = <int, _RimSection Function(double)>{};
        for (final id in ids) {
          final curve = s.edges[id].curve;
          final wall = s.faces[walls[id]!];
          final use = loop.coedges.firstWhere((c) => c.edge == id);
          final cache = <double, _RimSection>{};
          _RimSection station(double u) => cache.putIfAbsent(u, () {
            final t = curve.t0 + (curve.t1 - curve.t0) * u;
            final p = curve.evaluate(t);
            final tangent =
                curve.tangent(t) *
                (curve.t1 < curve.t0 ? -1.0 : 1.0) *
                (use.reversed ? -1.0 : 1.0);
            var normal = tangent.cross(n).normalized();
            // Saved/affine caps can wind either way. The adjacent wall, not
            // the cap's wire convention, determines the material side.
            final uv = wall.surface.parameters(p);
            if (uv != null &&
                normal.dot(wall.outwardNormal(uv.$1, uv.$2)) < 0) {
              normal = -normal;
            }
            final a = p - normal * distance, b = p - n * depth;
            final center = a - n * depth;
            return _RimSection(
              a,
              round ? a + b - center : (a + b) * .5,
              b,
              normal,
            );
          });
          for (var i = 0; i <= 16; i++) {
            final u = i / 16,
                p = curve.evaluate(curve.t0 + (curve.t1 - curve.t0) * u);
            final uv = wall.surface.parameters(p);
            if (uv == null ||
                (wall.outwardNormal(uv.$1, uv.$2) - station(u).normal).length >
                    1e-5 ||
                (wall.surface.closestPoint(station(u).b) - station(u).b)
                        .length >
                    tolerance) {
              eligible = false;
              break;
            }
          }
          if (!eligible) break;
          stations[id] = station;
        }
        if (!eligible) continue;
        try {
          final shifts = <int, Map<double, Vector3>>{};
          for (var i = 0; i < loop.coedges.length; i++) {
            final a = loop.coedges[i],
                b = loop.coedges[(i + 1) % loop.coedges.length];
            final ua = a.reversed ? 0.0 : 1.0, ub = b.reversed ? 1.0 : 0.0;
            final sa = stations[a.edge]!(ua), sb = stations[b.edge]!(ub);
            if ((sa.normal - sb.normal).length <= 1e-5) continue;
            if (!ids.every(selected.contains) ||
                s.edges[a.edge].curve is! BrepLine ||
                s.edges[b.edge].curve is! BrepLine) {
              throw const _NotExtrusionRim();
            }
            // Intersect the two offset lines. The correction is tangential
            // to each source line, so its blend still lies on the same
            // cylinder/plane. At the corner both patches share one rational
            // conic section, rather than leaving a gap between quarter arcs.
            final denominator = 1 + sa.normal.dot(sb.normal);
            if (denominator < 1e-8) {
              throw BrepUnsupported('Blend meets a folded-back profile corner');
            }
            final vertex = sa.b + n * depth;
            final contact =
                vertex - (sa.normal + sb.normal) * (distance / denominator);
            (shifts[a.edge] ??= {})[ua] = contact - sa.a;
            (shifts[b.edge] ??= {})[ub] = contact - sb.a;
          }
          for (final entry in shifts.entries) {
            final base = stations[entry.key]!;
            final start = entry.value[0] ?? Vector3.zero(),
                end = entry.value[1] ?? Vector3.zero();
            stations[entry.key] = (u) {
              final section = base(u), shift = start * (1 - u) + end * u;
              return _RimSection(
                section.a + shift,
                round ? section.middle : section.middle + shift * .5,
                section.b,
                section.normal,
              );
            };
          }
          for (final id in ids) {
            for (var i = 0; i <= 16; i++) {
              stations[id]!(i / 16);
            }
          }
          // Smooth joins and sharp miters must share all three control points.
          for (var i = 0; i < loop.coedges.length; i++) {
            final a = loop.coedges[i],
                b = loop.coedges[(i + 1) % loop.coedges.length];
            final sa = stations[a.edge]!(a.reversed ? 0 : 1);
            final sb = stations[b.edge]!(b.reversed ? 1 : 0);
            if ((sa.a - sb.a).length > tolerance ||
                (sa.b - sb.b).length > tolerance ||
                (sa.middle - sb.middle).length > tolerance)
              throw const _NotExtrusionRim();
          }
          final limiting = round
              ? _limitingCircularCorners(s, loop, stations, distance)
              : <int>{};
          _checkClearance(s, capId, loop, stations, distance, depth, limiting);
          return _rebuild(s, capId, loop, stations, round, limiting);
        } on _NotExtrusionRim {
          continue;
        }
      }
    }
    return null;
  }

  /// At exactly the convex corner radius, the inset circular arc is a point
  /// and the blend is spherical. This is a topology transition, not a reversed
  /// offset. Only analytic circular arcs with a proven constant inset qualify.
  static Set<int> _limitingCircularCorners(
    BrepSolid s,
    BrepLoop rim,
    Map<int, _RimSection Function(double)> stations,
    double radius,
  ) => {
    for (final co in rim.coedges)
      if (() {
        final c = s.edges[co.edge].curve;
        if (c is! BrepCircle ||
            c.isClosed ||
            (c.radius - radius).abs() > 1e-9 * math.max(1.0, radius)) {
          return false;
        }
        return [
          0.0,
          .25,
          .5,
          .75,
          1.0,
        ].every((u) => stations[co.edge]!(u).a.distanceTo(c.center) < 1e-8);
      }())
        co.edge,
  };

  static void _checkClearance(
    BrepSolid s,
    int capId,
    BrepLoop rim,
    Map<int, _RimSection Function(double)> stations,
    double d,
    double depth,
    Set<int> limiting,
  ) {
    final cap = s.faces[capId], plane = cap.surface as BrepPlane;
    List<Vector2> flatten(BrepCurve c) {
      final points = <Vector2>[];
      Vector2 project(Vector3 p) {
        final uv = plane.parameters(p)!;
        return Vector2(uv.$1, uv.$2);
      }

      void span(double a, double b, int depth) {
        final x = c.evaluate(a), y = c.evaluate(b), delta = y - x;
        var error = 0.0;
        for (final f in [.25, .5, .75]) {
          final p = c.evaluate(a + (b - a) * f),
              t = delta.length2 < 1e-24
                  ? 0.0
                  : (p - x).dot(delta) / delta.length2;
          error = math.max(error, (p - (x + delta * t.clamp(0, 1))).length);
        }
        if (error > math.min(tolerance * .1, d * .1) && depth < 16) {
          span(a, (a + b) * .5, depth + 1);
          span((a + b) * .5, b, depth + 1);
        } else {
          points.add(project(x));
        }
      }

      final breaks = {
        for (var i = 0; i <= 8; i++) c.t0 + (c.t1 - c.t0) * i / 8,
        if (c is BrepSplineCurve) ...c.knots.where((t) => t > c.t0 && t < c.t1),
      }.toList()..sort();
      for (var i = 1; i < breaks.length; i++) {
        span(breaks[i - 1], breaks[i], 0);
      }
      return points;
    }

    final loops = [
      for (final loop in cap.loops)
        <Vector2>[
          for (final co in loop.coedges)
            ...(() {
              final c = s.edges[co.edge].curve, points = flatten(c);
              if (!co.reversed) return points;
              final end = plane.parameters(c.evaluate(c.t1))!;
              return [
                Vector2(end.$1, end.$2),
                ...points.skip(1).toList().reversed,
              ];
            })(),
        ],
    ];
    bool inside(Vector3 p) {
      final exact = BrepClassify.containsPlaneFace(s, cap, p);
      if (exact != null) return exact;
      final uv = plane.parameters(p)!;
      final x = uv.$1, y = uv.$2;
      var inside = false;
      for (final ring in loops) {
        for (var i = 0, j = ring.length - 1; i < ring.length; j = i++) {
          final a = ring[i], b = ring[j];
          if ((a.y > y) != (b.y > y) &&
              x < (b.x - a.x) * (y - a.y) / (b.y - a.y) + a.x)
            inside = !inside;
        }
      }
      return inside;
    }

    for (final co in rim.coedges) {
      final edge = s.edges[co.edge], c = edge.curve;
      for (var i = 0; i <= 32; i++) {
        final section = stations[co.edge]!(i / 32);
        if (!inside(section.a))
          throw BrepUnsupported('Blend setback reaches another cap boundary');
        // A shrinking offset must retain the source traversal direction.
        final lo = math.max(0.0, i / 32 - 1e-4),
            hi = math.min(1.0, i / 32 + 1e-4);
        final tangent =
            c.tangent(c.t0 + (c.t1 - c.t0) * i / 32) * (c.t1 - c.t0);
        if (!limiting.contains(co.edge) &&
            (stations[co.edge]!(hi).a - stations[co.edge]!(lo).a).dot(
                  tangent,
                ) <=
                1e-12) {
          throw BrepUnsupported(
            'Blend would collapse or reverse a cap boundary; use a smaller radius',
          );
        }
      }
    }
    if (limiting.isNotEmpty) {
      // A sphere corner may consume an arc, but must not erase the whole cap
      // (for example a stadium inset collapsing to a line).
      final points = [
        for (final co in rim.coedges)
          for (var i = 0; i < 16; i++)
            stations[co.edge]!(co.reversed ? 1 - i / 16 : i / 16).a,
      ];
      final origin = points.first;
      var area = Vector3.zero();
      for (var i = 0; i < points.length; i++) {
        area += (points[i] - origin).cross(
          points[(i + 1) % points.length] - origin,
        );
      }
      if (area.length < 1e-10) {
        throw BrepUnsupported('Blend would collapse the entire cap boundary');
      }
    }
    // Moving the upper endpoints of a wall must not pass its opposite edge.
    final n = cap.outwardNormal(0, 0);
    for (final co in rim.coedges) {
      final e = s.edges[co.edge];
      for (final v in {e.v0, e.v1}) {
        for (final other in s.edges) {
          if (other.v0 != v && other.v1 != v ||
              rim.coedges.any((c) => identical(s.edges[c.edge], other)))
            continue;
          final far = s.vertices[other.v0 == v ? other.v1 : other.v0];
          final available = (s.vertices[v] - far).dot(n);
          if (available > 1e-7 && depth >= available - 1e-7)
            throw BrepUnsupported(
              'Blend setback reaches the opposite end of the wall',
            );
        }
      }
    }
  }

  static List<List<Vector3>> _controls(
    _RimSection Function(double) f,
    double a,
    double b,
  ) {
    final fa = f(a), fb = f(b), h = math.min(1e-5, (b - a) * .001);
    List<Vector3> derivative(double t) {
      final lo = math.max(0.0, t - h), hi = math.min(1.0, t + h);
      final x = f(lo).points, y = f(hi).points;
      return [for (var i = 0; i < 3; i++) (y[i] - x[i]) / (hi - lo)];
    }

    final da = derivative(a), db = derivative(b), x = fa.points, y = fb.points;
    return [
      x,
      [for (var i = 0; i < 3; i++) x[i] + da[i] * ((b - a) / 3)],
      [for (var i = 0; i < 3; i++) y[i] - db[i] * ((b - a) / 3)],
      y,
    ];
  }

  static BrepSolid _rebuild(
    BrepSolid s,
    int capId,
    BrepLoop rim,
    Map<int, _RimSection Function(double)> stations,
    bool round,
    Set<int> limiting,
  ) {
    final b = BrepBuilder();
    for (final v in s.vertices) {
      b.vertex(v);
    }
    final capChains = <int, List<int>>{}, wallChains = <int, List<int>>{};
    final capVertices = <int, int>{}, crossAtVertex = <int, int>{};
    final moved = <int, Vector3>{};
    final weight = round ? math.sqrt(.5) : 1.0;
    final patches = <BrepFace>[];
    final capParents = <int, int>{};
    int capRoot(int v) =>
        capParents.containsKey(v) ? capRoot(capParents[v]!) : v;
    final limitingVertices = <int>{};
    for (final id in limiting) {
      final e = s.edges[id];
      final a = capRoot(e.v0), z = capRoot(e.v1);
      if (a != z) capParents[z] = a;
      limitingVertices.addAll([e.v0, e.v1]);
    }
    int capVertex(int original, _RimSection section) {
      moved[original] = section.b;
      return capVertices.putIfAbsent(
        capRoot(original),
        () => b.vertex(section.a),
      );
    }

    int cross(int av, int bv, _RimSection section, {bool analytic = false}) =>
        b.edge(
          av,
          bv,
          round && analytic
              ? BrepCircle(
                  center: section.a + section.b - section.middle,
                  normal: (section.middle - section.b).cross(
                    section.middle - section.a,
                  ),
                  xDir: section.middle - section.b,
                  radius: (section.middle - section.b).length,
                  t0: 0,
                  t1: math.pi / 2,
                )
              : round
              ? BrepSplineCurve(
                  degree: 2,
                  points: section.points,
                  knots: [0, 0, 0, 1, 1, 1],
                  weights: [1, weight, 1],
                )
              : BrepLine(section.a, section.b),
        );
    for (final co in rim.coedges) {
      final e = s.edges[co.edge], f = stations[co.edge]!;
      if (limiting.contains(co.edge)) {
        final circle = e.curve as BrepCircle;
        final first = f(0), last = f(1);
        final apex = capVertex(e.v0, first);
        capVertex(e.v1, last);
        final center = first.a + first.b - first.middle;
        final sphere = BrepSphere(center: center, radius: circle.radius);
        final count = ((circle.t1 - circle.t0).abs() / (math.pi / 2)).ceil();
        var from = e.v0;
        var cross0 = crossAtVertex.putIfAbsent(
          e.v0,
          () => cross(apex, from, first, analytic: true),
        );
        final wall = <int>[];
        for (var i = 0; i < count; i++) {
          final u = (i + 1) / count, end = f(u);
          final to = i == count - 1 ? e.v1 : b.vertex(end.b);
          final cross1 = i == count - 1
              ? crossAtVertex.putIfAbsent(
                  e.v1,
                  () => cross(apex, to, end, analytic: true),
                )
              : cross(apex, to, end, analytic: true);
          final arc = BrepCircle(
            center: center,
            normal: circle.normal,
            xDir: circle.xDir,
            radius: circle.radius,
            t0: circle.t0 + (circle.t1 - circle.t0) * i / count,
            t1: circle.t0 + (circle.t1 - circle.t0) * u,
          );
          final edge = b.edge(from, to, arc);
          wall.add(edge);
          var loop = [
            BrepCoedge(cross0),
            BrepCoedge(edge),
            BrepCoedge(cross1, reversed: true),
          ];
          // Sphere trim winding selects the finite patch, independently of
          // the old cap's wire convention (including mirrored instances).
          final mid = arc.evaluate((arc.t0 + arc.t1) / 2);
          final sign = (first.a - center)
              .cross(mid - center)
              .dot(arc.tangent((arc.t0 + arc.t1) / 2) * (arc.t1 - arc.t0));
          if (sign < 0) {
            loop = [
              for (final c in loop.reversed)
                BrepCoedge(c.edge, reversed: !c.reversed),
            ];
          }
          patches.add(BrepFace(surface: sphere, loops: [BrepLoop(loop)]));
          from = to;
          cross0 = cross1;
        }
        capChains[co.edge] = [];
        wallChains[co.edge] = wall;
        continue;
      }
      final spans = <(double, double, List<List<Vector3>>)>[];
      void fit(double a, double z, int depth) {
        final net = _controls(f, a, z);
        var error = 0.0;
        for (final t in [.125, .25, .375, .5, .625, .75, .875]) {
          final expected = f(a + (z - a) * t).points;
          for (var j = 0; j < 3; j++) {
            final actual =
                net[0][j] * math.pow(1 - t, 3).toDouble() +
                net[1][j] * (3 * t * (1 - t) * (1 - t)) +
                net[2][j] * (3 * t * t * (1 - t)) +
                net[3][j] * (t * t * t);
            error = math.max(error, (actual - expected[j]).length);
          }
        }
        if (error > tolerance * .25) {
          if (depth >= 14 || spans.length > 4096)
            throw BrepUnsupported(
              'Blend offset could not meet the curve tolerance',
            );
          fit(a, (a + z) / 2, depth + 1);
          fit((a + z) / 2, z, depth + 1);
        } else {
          spans.add((a, z, net));
        }
      }

      // Keep adaptive spans in a shared spline patch instead of multiplying
      // topological faces (and their display sampling cost) at every fit split.
      final groups = <(double, double, List<List<Vector3>>, List<double>)>[];
      final seeds = e.curve is BrepLine ? 1 : 8;
      for (var i = 0; i < seeds; i++) {
        spans.clear();
        final curve = e.curve;
        final breaks = {
          i / seeds,
          (i + 1) / seeds,
          if (curve is BrepSplineCurve)
            ...curve.knots
                .map((t) => (t - curve.t0) / (curve.t1 - curve.t0))
                .where((u) => u > i / seeds && u < (i + 1) / seeds),
        }.toList()..sort();
        for (var j = 1; j < breaks.length; j++)
          fit(breaks[j - 1], breaks[j], 0);
        final net = <List<Vector3>>[], knots = <double>[0, 0, 0, 0];
        for (var j = 0; j < spans.length; j++) {
          net.addAll(j == 0 ? spans[j].$3 : spans[j].$3.skip(1));
          if (j < spans.length - 1)
            knots.addAll(List.filled(3, (spans[j].$2 - i / seeds) * seeds));
        }
        knots.addAll([1, 1, 1, 1]);
        groups.add((i / seeds, (i + 1) / seeds, net, knots));
      }
      final ca = <int>[], wa = <int>[];
      var av = capVertex(e.v0, f(0)), bv = e.v0;
      var startCross = crossAtVertex.putIfAbsent(
        e.v0,
        () => cross(av, bv, f(0), analytic: limitingVertices.contains(e.v0)),
      );
      for (var si = 0; si < groups.length; si++) {
        final span = groups[si],
            last = si == groups.length - 1,
            end = f(span.$2);
        final nextA = last ? capVertex(e.v1, end) : b.vertex(end.a);
        final nextB = last ? e.v1 : b.vertex(end.b);
        final endCross = last
            ? crossAtVertex.putIfAbsent(
                e.v1,
                () => cross(
                  nextA,
                  nextB,
                  end,
                  analytic: limitingVertices.contains(e.v1),
                ),
              )
            : cross(nextA, nextB, end);
        final net = span.$3;
        final linear = e.curve is BrepLine;
        final BrepCurve ac = linear
            ? BrepLine(net.first[0], net.last[0])
            : BrepSplineCurve(
                degree: 3,
                points: [for (final row in net) row[0]],
                knots: span.$4,
              );
        final BrepCurve bc = linear
            ? BrepLine(net.first[2], net.last[2])
            : BrepSplineCurve(
                degree: 3,
                points: [for (final row in net) row[2]],
                knots: span.$4,
              );
        final ae = b.edge(av, nextA, ac), be = b.edge(bv, nextB, bc);
        ca.add(ae);
        wa.add(be);
        final section = f((span.$1 + span.$2) * .5);
        final axis = e.curve
            .tangent((e.curve.t0 + e.curve.t1) * .5)
            .normalized();
        final radius = (section.middle - section.a).dot(section.normal).abs();
        final BrepSurface surface = linear
            ? (round
                  ? BrepCylinder(
                      origin:
                          section.a -
                          s.faces[capId].outwardNormal(0, 0) * radius,
                      axis: axis,
                      radius: radius,
                    )
                  : BrepPlane(
                      origin: section.a,
                      normal: axis.cross(section.b - section.a),
                    ))
            : round
            ? BrepSplineSurface(
                degreeU: 3,
                degreeV: 2,
                points: net,
                knotsU: span.$4,
                knotsV: [0, 0, 0, 1, 1, 1],
                weights: [
                  for (var i = 0; i < net.length; i++) [1, weight, 1],
                ],
              )
            : BrepRuledSurface(ac, bc);
        final loop = [
          BrepCoedge(ae),
          BrepCoedge(endCross),
          BrepCoedge(be, reversed: true),
          BrepCoedge(startCross, reversed: true),
        ];
        final expected =
            s.faces[capId].outwardNormal(0, 0) +
            f((span.$1 + span.$2) / 2).normal;
        final uv = linear
            ? surface.parameters((section.a + section.b) * .5)!
            : (.5, .5);
        patches.add(
          BrepFace(
            surface: surface,
            loops: [
              BrepLoop(
                co.reversed
                    ? loop
                    : [
                        for (final c in loop.reversed)
                          BrepCoedge(c.edge, reversed: !c.reversed),
                      ],
              ),
            ],
            reversed: surface.normalAt(uv.$1, uv.$2).dot(expected) < 0,
          ),
        );
        av = nextA;
        bv = nextB;
        startCross = endCross;
      }
      capChains[co.edge] = ca;
      wallChains[co.edge] = wa;
    }
    for (final entry in moved.entries) {
      b.vertices[entry.key] = entry.value;
    }
    final copied = <int, int>{};
    for (var i = 0; i < s.edges.length; i++) {
      if (capChains.containsKey(i)) continue;
      final e = s.edges[i];
      var curve = e.curve;
      if (moved.containsKey(e.v0) || moved.containsKey(e.v1)) {
        if (curve is! BrepLine)
          throw BrepUnsupported(
            'Rim endpoint meets a curved boundary that needs a vertex blend',
          );
        curve = BrepLine(b.vertices[e.v0], b.vertices[e.v1]);
      }
      copied[i] = b.edge(e.v0, e.v1, curve, origin: e.origin);
    }
    for (var fi = 0; fi < s.faces.length; fi++) {
      final face = s.faces[fi];
      b.face(
        face.surface,
        [
          for (final loop in face.loops)
            BrepLoop([
              for (final co in loop.coedges)
                ...(() {
                  final chain = fi == capId
                      ? capChains[co.edge]
                      : wallChains[co.edge];
                  if (chain == null)
                    return [
                      BrepCoedge(copied[co.edge]!, reversed: co.reversed),
                    ];
                  return [
                    for (final id in co.reversed ? chain.reversed : chain)
                      BrepCoedge(id, reversed: co.reversed),
                  ];
                })(),
            ]),
        ],
        reversed: face.reversed,
        origin: face.origin,
      );
    }
    b.faces.addAll(patches);
    return b.build();
  }
}

class _RimSection {
  _RimSection(this.a, this.middle, this.b, this.normal);
  final Vector3 a, middle, b, normal;
  List<Vector3> get points => [a, middle, b];
}

class _NotExtrusionRim implements Exception {
  const _NotExtrusionRim();
}
