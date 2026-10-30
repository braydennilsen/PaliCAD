import 'dart:math' as math;
import 'package:vector_math/vector_math_64.dart';
import 'geom.dart';
import 'solid.dart';
import 'rational_ruled.dart';
import 'spline_convert.dart';

/// Bounded rational surface intersections. Positive homogeneous weights give
/// conservative convex-hull boxes for subdivision; continuation corrects onto
/// both original surfaces. Singular or unresolved cells fail explicitly.
class BrepSplineIntersect {
  static List<BrepSurfaceIntersectionCurve> intersect(
    BrepSurface a,
    BrepSurface b,
  ) {
    BrepSplineSurface? rational(BrepSurface s) {
      if (s is BrepSplineSurface) return s;
      if (s is BrepTorus) return BrepSplineConvert.surface(s);
      if (s is BrepRuledSurface) {
        final exact = BrepSplineConvert.exactRuled(s);
        if (exact != null) return exact;
        if (s.lower is BrepHelix && s.upper is BrepHelix) {
          // Helices are not rational curves. Conversion has an analytic
          // Hermite error bound; keep it below the intersection tolerance.
          return BrepSplineConvert.surface(s, tolerance: 1e-10);
        }
      }
      if (s is BrepAffineSurface) {
        final base = rational(s.base);
        return base?.transformed(s.transform);
      }
      return null;
    }

    var x = rational(a), y = rational(b);
    BrepSplineSurface plane(BrepPlane p, BrepSplineSurface other) {
      final uv = [
        for (final row in other.points)
          for (final q in row) p.parameters(q)!,
      ];
      final u0 = uv.map((q) => q.$1).reduce(math.min),
          u1 = uv.map((q) => q.$1).reduce(math.max),
          v0 = uv.map((q) => q.$2).reduce(math.min),
          v1 = uv.map((q) => q.$2).reduce(math.max);
      final pad = math.max(1e-4, math.max(u1 - u0, v1 - v0) * .01);
      return BrepSplineSurface(
        degreeU: 1,
        degreeV: 1,
        points: [
          for (final u in [u0 - pad, u1 + pad])
            [
              for (final v in [v0 - pad, v1 + pad]) p.evaluate(u, v),
            ],
        ],
        knotsU: [0, 0, 1, 1],
        knotsV: [0, 0, 1, 1],
      );
    }

    if (x == null && a is BrepPlane && y != null) x = plane(a, y);
    if (y == null && b is BrepPlane && x != null) y = plane(b, x);
    if (x == null &&
        (a is BrepCylinder ||
            a is BrepCone ||
            a is BrepSphere ||
            a is BrepTorus) &&
        y != null)
      x = BrepSplineConvert.surface(
        a,
        boundary: y.points.expand((r) => r).toList(),
      );
    if (y == null &&
        (b is BrepCylinder ||
            b is BrepCone ||
            b is BrepSphere ||
            b is BrepTorus) &&
        x != null)
      y = BrepSplineConvert.surface(
        b,
        boundary: x.points.expand((r) => r).toList(),
      );
    if (x == null || y == null)
      throw BrepUnsupported(
        'Surface pair has no bounded intersection representation',
      );
    return surfaces(x, y);
  }

  /// Certifies that a positive-weight spline patch stays on one side of a
  /// quadric, allowing boundary-only contact within the topology tolerance.
  /// Bernstein coefficients bound the whole patch, not just sampled points.
  static bool boundaryContactOnly(
    BrepSurface a,
    BrepSurface b, {
    double tolerance = 1e-6,
  }) {
    final patch = a is BrepSplineSurface
        ? a
        : b is BrepSplineSurface
        ? b
        : null;
    final quadric = identical(patch, a) ? b : a;
    if (patch == null || (quadric is! BrepCylinder && quadric is! BrepSphere)) {
      return false;
    }
    final radius = quadric is BrepCylinder
        ? quadric.radius
        : (quadric as BrepSphere).radius;
    // For spheres/cylinders, |F| <= 2*r*eps - eps² is strictly inside
    // the geometric distance tolerance on either side of the surface.
    if (radius <= tolerance) return false;
    final implicitTolerance = 2 * radius * tolerance - tolerance * tolerance;
    double choose(int n, int k) {
      var out = 1.0;
      for (var i = 1; i <= k; i++) {
        out *= (n - i + 1) / i;
      }
      return out;
    }

    (double, double, double) range(_SplineCell cell) {
      final u = cell.control.length - 1, v = cell.control.first.length - 1;
      final coefficients = List.generate(
        2 * u + 1,
        (_) => List.filled(2 * v + 1, 0.0),
      );
      final points = [
        for (final row in cell.control)
          [for (final h in row) Vector3(h.x, h.y, h.z) / h.w],
      ];
      var minWeight = double.infinity;
      for (var i = 0; i <= u; i++) {
        for (var j = 0; j <= v; j++) {
          final p = points[i][j], w = cell.control[i][j].w;
          minWeight = math.min(minWeight, w);
          for (var k = 0; k <= u; k++) {
            for (var l = 0; l <= v; l++) {
              final q = points[k][l];
              final bilinear =
                  (BrepQuadric.value(quadric, p) +
                      BrepQuadric.value(quadric, q) -
                      BrepQuadric.line(quadric, p, q - p).$1) *
                  .5;
              coefficients[i + k][j + l] +=
                  bilinear *
                  w *
                  cell.control[k][l].w *
                  choose(u, i) *
                  choose(u, k) /
                  choose(2 * u, i + k) *
                  choose(v, j) *
                  choose(v, l) /
                  choose(2 * v, j + l);
            }
          }
        }
      }
      final values = coefficients.expand((r) => r);
      return (
        values.reduce(math.min),
        values.reduce(math.max),
        implicitTolerance * minWeight * minWeight,
      );
    }

    final cells = _SplineCell.spans(patch);
    // A coincident patch is not a boundary-only contact certificate.
    final ranges = cells.map(range).toList();
    final high = ranges.map((r) => r.$2).reduce(math.max),
        low = ranges.map((r) => r.$1).reduce(math.min);
    final sign = high.abs() >= low.abs() ? 1 : -1;
    if (ranges.every((r) => math.max(r.$1.abs(), r.$2.abs()) <= r.$3)) {
      return false;
    }
    var visited = 0;
    bool prove(_SplineCell cell, int depth) {
      if (++visited > 4096) return false;
      final r = range(cell);
      if (sign > 0 ? r.$1 >= -r.$3 : r.$2 <= r.$3) return true;
      if (depth >= 16 || (sign > 0 ? r.$2 < -r.$3 : r.$1 > r.$3)) return false;
      return cell.split().every((c) => prove(c, depth + 1));
    }

    return cells.every((c) => prove(c, 0));
  }

  /// Verified ray/surface roots in the bounded spline domain. Ambiguous
  /// tangencies return null so material classification can choose another ray.
  static List<double>? ray(
    BrepSplineSurface s,
    Vector3 origin,
    Vector3 direction, {
    double tolerance = 1e-8,
  }) {
    final n = brepFrameU(direction), m = direction.cross(n), out = <double>[];
    var visited = 0, ambiguous = false;
    void search(_SplineCell cell) {
      if (ambiguous || !cell.rayBox(origin, direction, tolerance)) return;
      if (++visited > 20000) {
        ambiguous = true;
        return;
      }
      if (cell.width <= 1 / 16) {
        var u = (cell.u0 + cell.u1) / 2, v = (cell.v0 + cell.v1) / 2;
        var converged = false;
        for (var i = 0; i < 24; i++) {
          final j = BrepSurfaceIntersectionCurve.normalizedJet(s, u, v),
              p = j.$1 - origin;
          final x = p.dot(n), y = p.dot(m);
          if (x * x + y * y < tolerance * tolerance) {
            converged = true;
            break;
          }
          final a = j.$2.dot(n),
              b = j.$3.dot(n),
              c = j.$2.dot(m),
              d = j.$3.dot(m),
              det = a * d - b * c;
          if (det.abs() < 1e-14 * j.$2.length * j.$3.length) break;
          final du = (d * x - b * y) / det, dv = (a * y - c * x) / det;
          var step = 1.0, accepted = false;
          while (step >= 1 / 1024) {
            final nu = u - du * step, nv = v - dv * step;
            if (nu >= -1e-8 &&
                nu <= 1 + 1e-8 &&
                nv >= -1e-8 &&
                nv <= 1 + 1e-8) {
              final q = BrepSurfaceIntersectionCurve.point(s, nu, nv) - origin;
              if (q.dot(n) * q.dot(n) + q.dot(m) * q.dot(m) < x * x + y * y) {
                u = nu;
                v = nv;
                accepted = true;
                break;
              }
            }
            step *= .5;
          }
          if (!accepted) break;
        }
        if (converged && cell.contains(u, v, padding: cell.width)) {
          final j = BrepSurfaceIntersectionCurve.normalizedJet(s, u, v);
          if (j.$2.cross(j.$3).normalized().dot(direction).abs() < 1e-8) {
            ambiguous = true;
            return;
          }
          final t = (j.$1 - origin).dot(direction);
          if (!out.any((x) => (x - t).abs() < tolerance * 4)) out.add(t);
          return;
        }
        if (cell.width < 1 / 1024) {
          ambiguous = true;
          return;
        }
      }
      for (final child in cell.split()) search(child);
    }

    for (final cell in _SplineCell.spans(s)) search(cell);
    return ambiguous ? null : out;
  }

  static List<BrepSurfaceIntersectionCurve> surfaces(
    BrepSplineSurface a,
    BrepSplineSurface b, {
    double tolerance = 1e-9,
  }) {
    if (!tolerance.isFinite || tolerance <= 0)
      throw ArgumentError.value(tolerance, 'tolerance');
    for (final surface in [a, b]) {
      for (final domain in [
        (surface.knotsU, surface.degreeU, surface.u0, surface.u1),
        (surface.knotsV, surface.degreeV, surface.v0, surface.v1),
      ]) {
        final counts = <double, int>{};
        for (final knot in domain.$1)
          if (knot > domain.$3 && knot < domain.$4)
            counts[knot] = (counts[knot] ?? 0) + 1;
        if (counts.values.any((n) => n > domain.$2))
          throw BrepUnsupported(
            'Discontinuous spline spans must be separate faces before intersection',
          );
      }
    }

    final out = <BrepSurfaceIntersectionCurve>[];
    var visited = 0;
    bool covered(Vector4 x) {
      for (final c in out)
        for (var i = 1; i < c.anchors.length; i++) {
          final p = c.anchors[i - 1], d = c.anchors[i] - p;
          final t = ((x - p).dot(d) / d.length2).clamp(0.0, 1.0);
          if ((p + d * t - x).length < d.length * .5 + 1e-7) {
            final y = BrepSurfaceIntersectionCurve.correct(
              a,
              b,
              p + d * t,
              x,
              d.normalized(),
              tolerance,
            );
            if (y != null && (y - x).length < 1e-7) return true;
          }
        }
      return false;
    }

    void search(_SplineCell p, _SplineCell q, int depth) {
      if (++visited > 200000)
        throw BrepUnsupported(
          'Surface intersection subdivision budget exhausted',
        );
      if (!p.overlaps(q, tolerance)) return;
      final width = math.max(p.width, q.width);
      if (width <= 1 / 16) {
        final seed = Vector4(
          (p.u0 + p.u1) / 2,
          (p.v0 + p.v1) / 2,
          (q.u0 + q.u1) / 2,
          (q.v0 + q.v1) / 2,
        );
        Vector4? hit;
        try {
          final dir = BrepSurfaceIntersectionCurve.direction(a, b, seed);
          hit = BrepSurfaceIntersectionCurve.correct(
            a,
            b,
            seed,
            seed,
            dir,
            tolerance,
          );
        } on BrepUnsupported {
          /* Subdivide before diagnosing a singular cell. */
        }
        if (hit != null &&
            p.contains(hit.x, hit.y, padding: width) &&
            q.contains(hit.z, hit.w, padding: width)) {
          if (!covered(hit)) {
            final points = _trace(a, b, hit, tolerance);
            out.add(
              BrepSurfaceIntersectionCurve(
                first: a,
                second: b,
                anchors: points,
                tolerance: tolerance,
              ),
            );
          }
          return;
        }
        // A trace passing through this box also establishes a local branch.
        if (out.any(
          (c) => c.anchors.any(
            (x) => p.contains(x.x, x.y) && q.contains(x.z, x.w),
          ),
        ))
          return;
        if (width < 1 / 1024 || depth > 48)
          throw BrepUnsupported(
            'Surface intersection contains an unresolved contact',
          );
      }
      if (p.width >= q.width) {
        for (final child in p.split()) search(child, q, depth + 1);
      } else {
        for (final child in q.split()) search(p, child, depth + 1);
      }
    }

    for (final p in _SplineCell.spans(a))
      for (final q in _SplineCell.spans(b)) search(p, q, 0);
    return out;
  }

  static List<Vector4> _trace(
    BrepSplineSurface a,
    BrepSplineSurface b,
    Vector4 start,
    double tolerance,
  ) {
    List<Vector4> walk(int sign) {
      final points = <Vector4>[start];
      var x = start,
          dir =
              BrepSurfaceIntersectionCurve.direction(a, b, start) *
              sign.toDouble();
      for (var k = 0; k < 8000; k++) {
        var step = .02;
        int? boundary;
        for (var i = 0; i < 4; i++) {
          if (dir[i].abs() < 1e-12) continue;
          final h = ((dir[i] > 0 ? 1.0 : 0.0) - x[i]) / dir[i];
          if (h < step) {
            step = math.max(0.0, h);
            boundary = i;
          }
        }
        if (step < 1e-9) return points;
        Vector4? y;
        for (var attempt = 0; attempt < 12; attempt++) {
          final prediction = x + dir * step;
          final normal = boundary == null ? dir : Vector4.zero();
          if (boundary != null) normal[boundary] = 1;
          y = BrepSurfaceIntersectionCurve.correct(
            a,
            b,
            prediction,
            prediction,
            normal,
            tolerance,
          );
          if (y != null && (y - prediction).length < step * .3 + 1e-9) break;
          y = null;
          step *= .5;
          boundary = null;
        }
        if (y == null)
          throw BrepUnsupported('Surface intersection continuation failed');
        if (points.length > 8 && (y - start).length < step * 1.1) {
          points.add(start);
          return points;
        }
        if ((y - x).length < 1e-10)
          throw BrepUnsupported('Surface intersection continuation stalled');
        points.add(y);
        x = y;
        if (boundary != null) return points;
        var next = BrepSurfaceIntersectionCurve.direction(a, b, x);
        if (next.dot(dir) < 0) next = -next;
        dir = next;
      }
      throw BrepUnsupported(
        'Surface intersection continuation budget exhausted',
      );
    }

    final forward = walk(1);
    if (forward.length > 2 && (forward.last - start).length < 1e-9)
      return forward;
    final backward = walk(-1);
    return [...backward.reversed.skipLastOne(), ...forward];
  }
}

extension<T> on Iterable<T> {
  Iterable<T> skipLastOne() sync* {
    final a = toList();
    for (var i = 0; i < a.length - 1; i++) yield a[i];
  }
}

class _SplineCell {
  _SplineCell(this.control, this.u0, this.u1, this.v0, this.v1) {
    lo = Vector3.all(double.infinity);
    hi = Vector3.all(-double.infinity);
    for (final row in control)
      for (final h in row)
        for (var i = 0; i < 3; i++) {
          final p = h[i] / h.w;
          lo[i] = math.min(lo[i], p);
          hi[i] = math.max(hi[i], p);
        }
  }
  final List<List<Vector4>> control;
  final double u0, u1, v0, v1;
  late final Vector3 lo, hi;
  double get width => math.max(u1 - u0, v1 - v0);
  bool contains(double u, double v, {double padding = 1e-7}) =>
      u >= u0 - padding &&
      u <= u1 + padding &&
      v >= v0 - padding &&
      v <= v1 + padding;
  bool overlaps(_SplineCell q, double tolerance) => List.generate(
    3,
    (i) => i,
  ).every((i) => hi[i] >= q.lo[i] - tolerance && q.hi[i] >= lo[i] - tolerance);
  bool rayBox(Vector3 p, Vector3 d, double tolerance) {
    var near = double.negativeInfinity, far = double.infinity;
    for (var i = 0; i < 3; i++) {
      if (d[i].abs() < 1e-15) {
        if (p[i] < lo[i] - tolerance || p[i] > hi[i] + tolerance) return false;
        continue;
      }
      final a = (lo[i] - tolerance - p[i]) / d[i],
          b = (hi[i] + tolerance - p[i]) / d[i];
      near = math.max(near, math.min(a, b));
      far = math.min(far, math.max(a, b));
      if (near > far) return false;
    }
    return true;
  }

  List<_SplineCell> split() {
    final u = u1 - u0 >= v1 - v0;
    final grid = u
        ? control
        : [
            for (var j = 0; j < control.first.length; j++)
              [for (final row in control) row[j]],
          ];
    final left = List.generate(grid.length, (_) => <Vector4>[]),
        right = List.generate(grid.length, (_) => <Vector4>[]);
    for (var j = 0; j < grid.first.length; j++) {
      var row = [for (final r in grid) r[j]];
      final l = <Vector4>[row.first], r = <Vector4>[row.last];
      while (row.length > 1) {
        row = [for (var i = 1; i < row.length; i++) (row[i - 1] + row[i]) * .5];
        l.add(row.first);
        r.add(row.last);
      }
      for (var i = 0; i < grid.length; i++) {
        left[i].add(l[i]);
        right[i].add(r[r.length - 1 - i]);
      }
    }
    List<List<Vector4>> transpose(List<List<Vector4>> x) => [
      for (var j = 0; j < x.first.length; j++) [for (final row in x) row[j]],
    ];
    if (u)
      return [
        _SplineCell(left, u0, (u0 + u1) / 2, v0, v1),
        _SplineCell(right, (u0 + u1) / 2, u1, v0, v1),
      ];
    return [
      _SplineCell(transpose(left), u0, u1, v0, (v0 + v1) / 2),
      _SplineCell(transpose(right), u0, u1, (v0 + v1) / 2, v1),
    ];
  }

  static List<_SplineCell> spans(BrepSplineSurface s) {
    final us = s.knotsU.where((k) => k >= s.u0 && k <= s.u1).toSet().toList()
      ..sort();
    final vs = s.knotsV.where((k) => k >= s.v0 && k <= s.v1).toSet().toList()
      ..sort();
    final out = <_SplineCell>[];
    for (var ui = 1; ui < us.length; ui++) {
      final columns = [
        for (var j = 0; j < s.points.first.length; j++)
          BrepRationalRuled.homogeneousSpan(
            BrepSplineCurve(
              degree: s.degreeU,
              points: [for (final r in s.points) r[j]],
              knots: s.knotsU,
              weights: [for (final r in s.weights) r[j]],
              start: us[ui - 1],
              end: us[ui],
            ),
            0,
            1,
          ),
      ];
      for (var vi = 1; vi < vs.length; vi++) {
        final rows = [
          for (var i = 0; i <= s.degreeU; i++)
            BrepRationalRuled.homogeneousSpan(
              BrepSplineCurve(
                degree: s.degreeV,
                points: [
                  for (final col in columns)
                    Vector3(col[i].x, col[i].y, col[i].z) / col[i].w,
                ],
                weights: [for (final col in columns) col[i].w],
                knots: s.knotsV,
                start: vs[vi - 1],
                end: vs[vi],
              ),
              0,
              1,
            ),
        ];
        out.add(
          _SplineCell(
            rows,
            (us[ui - 1] - s.u0) / (s.u1 - s.u0),
            (us[ui] - s.u0) / (s.u1 - s.u0),
            (vs[vi - 1] - s.v0) / (s.v1 - s.v0),
            (vs[vi] - s.v0) / (s.v1 - s.v0),
          ),
        );
      }
    }
    return out;
  }
}
