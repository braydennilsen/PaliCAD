part of 'geom.dart';

/// Rational B-spline geometry, independent of display triangles.
void _splineCheck(int degree, List<double> knots, int count) {
  if (degree < 1 ||
      count <= degree ||
      knots.length != count + degree + 1 ||
      knots.any((v) => !v.isFinite))
    throw BrepUnsupported('Invalid spline degree or knots');
  for (var i = 1; i < knots.length; i++) {
    if (knots[i] < knots[i - 1])
      throw BrepUnsupported('Spline knots must be ordered');
  }
  if (knots[degree] >= knots[count])
    throw BrepUnsupported('Collapsed spline domain');
}

List<double> _splineBasis(
  int degree,
  List<double> knots,
  int count,
  double value,
) {
  final t = value.clamp(knots[degree], knots[count]);
  var span = degree;
  if (t >= knots[count]) {
    span = count - 1;
  } else {
    while (span < count - 1 && t >= knots[span + 1]) {
      span++;
    }
  }
  var basis = List<double>.filled(knots.length - 1, 0);
  basis[span] = 1;
  for (var d = 1; d <= degree; d++) {
    final next = List<double>.filled(knots.length - d - 1, 0);
    for (var i = math.max(0, span - d); i <= span && i < next.length; i++) {
      final a = knots[i + d] - knots[i], b = knots[i + d + 1] - knots[i + 1];
      next[i] =
          (a == 0 ? 0 : (t - knots[i]) * basis[i] / a) +
          (b == 0 ? 0 : (knots[i + d + 1] - t) * basis[i + 1] / b);
    }
    basis = next;
  }
  return basis.take(count).toList();
}

List<double> _splineBasisDerivative(
  int degree,
  List<double> knots,
  int count,
  double t,
) {
  final lower = _splineBasis(
    degree - 1,
    knots.sublist(1, knots.length - 1),
    count - 1,
    t,
  );
  final out = List<double>.filled(count, 0);
  for (var i = 0; i < count - 1; i++) {
    final span = knots[i + degree + 1] - knots[i + 1];
    final value = span == 0 ? 0.0 : degree * lower[i] / span;
    out[i] -= value;
    out[i + 1] += value;
  }
  return out;
}

void _splinePoints(Iterable<Vector3> points, Iterable<double> weights) {
  if (points.any((p) => ![p.x, p.y, p.z].every((v) => v.isFinite)) ||
      weights.any((w) => !w.isFinite || w <= 0))
    throw BrepUnsupported('Spline points and positive weights must be finite');
}

double _curveNearest(
  BrepCurve c,
  Vector3 p,
  int count, {
  List<Vector3>? samples,
}) {
  final lo = math.min(c.t0, c.t1), hi = math.max(c.t0, c.t1);
  final distances = [
    for (var i = 0; i <= count; i++)
      ((samples?[i] ?? c.evaluate(lo + (hi - lo) * i / count)) - p).length2,
  ];
  var best = distances.first < distances.last ? lo : hi,
      bestD = math.min(distances.first, distances.last);
  for (var i = 0; i <= count; i++) {
    if (i > 0 && distances[i] > distances[i - 1] ||
        i < count && distances[i] > distances[i + 1])
      continue;
    var a = lo + (hi - lo) * math.max(0, i - 1) / count,
        b = lo + (hi - lo) * math.min(count, i + 1) / count;
    for (var k = 0; k < 48; k++) {
      final x = a + (b - a) / 3, y = b - (b - a) / 3;
      if ((c.evaluate(x) - p).length2 < (c.evaluate(y) - p).length2) {
        b = y;
      } else {
        a = x;
      }
    }
    final t = (a + b) / 2, d = (c.evaluate(t) - p).length2;
    if (d < bestD) {
      best = t;
      bestD = d;
    }
  }
  return best;
}

/// Reuses only the seed samples for nearest-curve searches. Refinement still
/// evaluates the original curve. Keep a locator within one operation so edits
/// to caller-owned geometry cannot leave a stale cross-operation cache.
class BrepCurveLocator {
  BrepCurveLocator(this.curve)
    : count = curve is BrepSplineCurve
          ? math.max(24, curve.points.length * 3)
          : curve is BrepSurfaceIntersectionCurve
          ? math.max(32, curve.anchors.length)
          : 32;
  final BrepCurve curve;
  final int count;
  late final List<Vector3> samples = [
    for (var i = 0; i <= count; i++)
      curve.evaluate(
        math.min(curve.t0, curve.t1) + (curve.t1 - curve.t0).abs() * i / count,
      ),
  ];
  double closestParam(Vector3 p) =>
      _curveNearest(curve, p, count, samples: samples);
}

double _splineEndpoint(int degree, List<double> knots, int count, bool end) {
  _splineCheck(degree, knots, count);
  return knots[end ? count : degree];
}

class BrepSplineCurve implements BrepCurve {
  BrepSplineCurve({
    required this.degree,
    required List<Vector3> points,
    required List<double> knots,
    List<double>? weights,
    double? start,
    double? end,
  }) : points = List.unmodifiable(points.map((p) => p.clone())),
       knots = List.unmodifiable(knots),
       weights = List.unmodifiable(weights ?? List.filled(points.length, 1.0)),
       t0 = start ?? _splineEndpoint(degree, knots, points.length, false),
       t1 = end ?? _splineEndpoint(degree, knots, points.length, true) {
    _splineCheck(degree, knots, points.length);
    if (this.weights.length != points.length)
      throw BrepUnsupported('Spline weights do not match points');
    _splinePoints(points, this.weights);
    final lo = knots[degree], hi = knots[points.length], period = hi - lo;
    if (!t0.isFinite ||
        !t1.isFinite ||
        (isClosed
            ? math.min(t0, t1) < lo - period - 1e-10 ||
                  math.max(t0, t1) > hi + period + 1e-10 ||
                  (t1 - t0).abs() > period + 1e-10
            : math.min(t0, t1) < lo - 1e-10 || math.max(t0, t1) > hi + 1e-10)) {
      throw BrepUnsupported('Spline trim outside knot domain');
    }
  }
  final int degree;
  final List<Vector3> points;
  final List<double> knots, weights;
  @override
  final double t0, t1;
  bool get isClosed => (points.first - points.last).length < 1e-9;
  double _parameter(double t) {
    final lo = knots[degree], hi = knots[points.length];
    return isClosed && (t < lo || t > hi) ? lo + (t - lo) % (hi - lo) : t;
  }

  @override
  Vector3 evaluate(double t) {
    final b = _splineBasis(degree, knots, points.length, _parameter(t));
    var p = Vector3.zero(), w = 0.0;
    for (var i = 0; i < points.length; i++) {
      if (b[i] == 0) continue;
      final a = b[i] * weights[i];
      p += points[i] * a;
      w += a;
    }
    if (w <= 0) throw BrepUnsupported('Undefined rational curve');
    return p / w;
  }

  Vector3 derivative(double t) {
    t = _parameter(t);
    final b = _splineBasis(degree, knots, points.length, t),
        d = _splineBasisDerivative(degree, knots, points.length, t);
    var h = Vector3.zero(), dh = Vector3.zero(), w = 0.0, dw = 0.0;
    for (var i = 0; i < points.length; i++) {
      h += points[i] * (b[i] * weights[i]);
      dh += points[i] * (d[i] * weights[i]);
      w += b[i] * weights[i];
      dw += d[i] * weights[i];
    }
    return (dh - h * (dw / w)) / w;
  }

  @override
  Vector3 tangent(double t) => derivative(t).normalized();
  @override
  double closestParam(Vector3 p) =>
      _curveNearest(this, p, math.max(24, points.length * 3));
  @override
  Vector3 closestPoint(Vector3 p) => evaluate(closestParam(p));
  BrepSplineCurve trimmed(double a, double b) => BrepSplineCurve(
    degree: degree,
    points: points,
    knots: knots,
    weights: weights,
    start: a,
    end: b,
  );
  BrepSplineCurve transformed(BrepAffineTransform m) => BrepSplineCurve(
    degree: degree,
    points: points.map(m.point).toList(),
    knots: knots,
    weights: weights,
    start: t0,
    end: t1,
  );
  @override
  Map<String, dynamic> toJson() => {
    't': 'spline',
    'degree': degree,
    'points': points.map(_v).toList(),
    'knots': knots,
    'weights': weights,
    'a': t0,
    'b': t1,
  };
  factory BrepSplineCurve.fromJson(Map<String, dynamic> j) => BrepSplineCurve(
    degree: j['degree'] as int,
    points: (j['points'] as List).map(_p).toList(),
    knots: _splineDoubles(j['knots']),
    weights: _splineDoubles(j['weights']),
    start: (j['a'] as num).toDouble(),
    end: (j['b'] as num).toDouble(),
  );
}

List<double> _splineDoubles(dynamic xs) =>
    (xs as List).map((v) => (v as num).toDouble()).toList();

class BrepSplineSurface implements BrepSurface {
  BrepSplineSurface({
    required this.degreeU,
    required this.degreeV,
    required List<List<Vector3>> points,
    required List<double> knotsU,
    required List<double> knotsV,
    List<List<double>>? weights,
  }) : points = List.unmodifiable(
         points.map((r) => List<Vector3>.unmodifiable(r.map((p) => p.clone()))),
       ),
       knotsU = List.unmodifiable(knotsU),
       knotsV = List.unmodifiable(knotsV),
       weights = List.unmodifiable(
         (weights ??
                 [for (final row in points) List<double>.filled(row.length, 1)])
             .map((r) => List<double>.unmodifiable(r)),
       ) {
    if (points.isEmpty ||
        points.first.isEmpty ||
        points.any((r) => r.length != points.first.length))
      throw BrepUnsupported('Spline grid must be rectangular');
    _splineCheck(degreeU, knotsU, points.length);
    _splineCheck(degreeV, knotsV, points.first.length);
    if (this.weights.length != points.length ||
        this.weights.any((r) => r.length != points.first.length))
      throw BrepUnsupported('Spline weights do not match grid');
    _splinePoints(points.expand((r) => r), this.weights.expand((r) => r));
  }
  final int degreeU, degreeV;
  final List<List<Vector3>> points;
  final List<List<double>> weights;
  final List<double> knotsU, knotsV;
  double get u0 => knotsU[degreeU];
  double get u1 => knotsU[points.length];
  double get v0 => knotsV[degreeV];
  double get v1 => knotsV[points.first.length];
  BrepRuledSurface? get ruled {
    if (degreeV != 1 ||
        points.first.length != 2 ||
        List.generate(points.length, (i) => i).any(
          (i) => (weights[i][0] - weights[i][1]).abs() > 1e-12 * weights[i][0],
        ))
      return null;
    return BrepRuledSurface(
      BrepSplineCurve(
        degree: degreeU,
        points: [for (final r in points) r.first],
        knots: knotsU,
        weights: [for (final r in weights) r.first],
      ),
      BrepSplineCurve(
        degree: degreeU,
        points: [for (final r in points) r.last],
        knots: knotsU,
        weights: [for (final r in weights) r.last],
      ),
    );
  }

  bool get isClosed {
    bool same(Vector3 a, Vector3 b) => (a - b).length < 1e-7;
    final uClosed = List.generate(9, (i) => i / 8).every(
      (t) => same(
        evaluate(u0, v0 + (v1 - v0) * t),
        evaluate(u1, v0 + (v1 - v0) * t),
      ),
    );
    final vClosed = List.generate(9, (i) => i / 8).every(
      (t) => same(
        evaluate(u0 + (u1 - u0) * t, v0),
        evaluate(u0 + (u1 - u0) * t, v1),
      ),
    );
    final uPoles = [u0, u1].every(
      (u) => List.generate(
        9,
        (i) => i / 8,
      ).every((t) => same(evaluate(u, v0), evaluate(u, v0 + (v1 - v0) * t))),
    );
    final vPoles = [v0, v1].every(
      (v) => List.generate(
        9,
        (i) => i / 8,
      ).every((t) => same(evaluate(u0, v), evaluate(u0 + (u1 - u0) * t, v))),
    );
    return (uClosed && (vClosed || vPoles)) || (vClosed && uPoles);
  }

  @override
  Vector3 evaluate(double u, double v) {
    final a = _splineBasis(degreeU, knotsU, points.length, u),
        b = _splineBasis(degreeV, knotsV, points.first.length, v);
    var p = Vector3.zero(), w = 0.0;
    for (var i = 0; i < a.length; i++) {
      if (a[i] == 0) continue;
      for (var j = 0; j < b.length; j++) {
        final c = a[i] * b[j] * weights[i][j];
        p += points[i][j] * c;
        w += c;
      }
    }
    if (w <= 0) throw BrepUnsupported('Undefined rational surface');
    return p / w;
  }

  /// Position and exact rational first derivatives in the knot coordinates.
  (Vector3, Vector3, Vector3) jet(double u, double v) {
    final a = _splineBasis(degreeU, knotsU, points.length, u),
        b = _splineBasis(degreeV, knotsV, points.first.length, v),
        da = _splineBasisDerivative(degreeU, knotsU, points.length, u),
        db = _splineBasisDerivative(degreeV, knotsV, points.first.length, v);
    var h = Vector3.zero(), hu = Vector3.zero(), hv = Vector3.zero();
    var w = 0.0, wu = 0.0, wv = 0.0;
    for (var i = 0; i < a.length; i++) {
      if (a[i] == 0 && da[i] == 0) continue;
      for (var j = 0; j < b.length; j++) {
        if (b[j] == 0 && db[j] == 0) continue;
        final c = a[i] * b[j] * weights[i][j],
            cu = da[i] * b[j] * weights[i][j],
            cv = a[i] * db[j] * weights[i][j];
        h += points[i][j] * c;
        hu += points[i][j] * cu;
        hv += points[i][j] * cv;
        w += c;
        wu += cu;
        wv += cv;
      }
    }
    final p = h / w;
    return (p, (hu - p * wu) / w, (hv - p * wv) / w);
  }

  Vector3 du(double u, double v) => jet(u, v).$2;
  Vector3 dv(double u, double v) => jet(u, v).$3;

  @override
  Vector3 normalAt(double u, double v) {
    var n = du(u, v).cross(dv(u, v));
    if (n.length2 < 1e-24) {
      final a = u.clamp(u0 + (u1 - u0) * 1e-7, u1 - (u1 - u0) * 1e-7),
          b = v.clamp(v0 + (v1 - v0) * 1e-7, v1 - (v1 - v0) * 1e-7);
      n = du(a, b).cross(dv(a, b));
    }
    return n.normalized();
  }

  @override
  (double, double) parameters(Vector3 p) => _parameters(p, null);

  (double, double) _parameters(
    Vector3 p,
    List<(double, double, Vector3)>? samples,
  ) {
    var best = (u0, v0), bestD = double.infinity;
    final seeds = <(double, double, double)>[];
    final grid = samples ?? _parameterSamples();
    for (final sample in grid) {
      seeds.add((sample.$1, sample.$2, (sample.$3 - p).length2));
    }
    seeds.sort((a, b) => a.$3.compareTo(b.$3));
    for (final seed in seeds.take(4)) {
      var u = seed.$1, v = seed.$2;
      for (var k = 0; k < 32; k++) {
        final (q, x, y) = jet(u, v);
        final r = q - p,
            xx = x.length2,
            xy = x.dot(y),
            yy = y.length2,
            det = xx * yy - xy * xy;
        if (det <= 1e-20 * xx * yy) break;
        final su = (yy * x.dot(r) - xy * y.dot(r)) / det,
            sv = (xx * y.dot(r) - xy * x.dot(r)) / det;
        var step = 1.0, moved = false;
        while (step >= 1 / 1024) {
          final nu = (u - su * step).clamp(u0, u1),
              nv = (v - sv * step).clamp(v0, v1);
          if ((evaluate(nu, nv) - p).length2 < r.length2) {
            moved = (nu - u).abs() + (nv - v).abs() > 1e-12;
            u = nu;
            v = nv;
            break;
          }
          step *= .5;
        }
        if (!moved) break;
      }
      final d = (evaluate(u, v) - p).length2;
      if (d < bestD) {
        best = (u, v);
        bestD = d;
      }
    }
    return best;
  }

  // A uniform grid can alias a multi-turn thread to the wrong revolution.
  // Include every knot span so each polynomial patch can seed the inverse.
  List<(double, double, Vector3)> _parameterSamples() {
    List<double> stations(List<double> knots, double start, double end) {
      final values = <double>{
        for (var i = 0; i <= 6; i++) start + (end - start) * i / 6,
      };
      for (var i = 0; i + 1 < knots.length; i++) {
        final a = knots[i], b = knots[i + 1];
        if (a < start || b > end || b <= a) continue;
        values.addAll([a, (a + b) * .5, b]);
      }
      return values.toList()..sort();
    }

    return [
      for (final u in stations(knotsU, u0, u1))
        for (final v in stations(knotsV, v0, v1)) (u, v, evaluate(u, v)),
    ];
  }

  @override
  Vector3 closestPoint(Vector3 p) {
    final uv = parameters(p);
    return evaluate(uv.$1, uv.$2);
  }

  BrepSplineSurface transformed(BrepAffineTransform m) => BrepSplineSurface(
    degreeU: degreeU,
    degreeV: degreeV,
    points: [for (final row in points) row.map(m.point).toList()],
    knotsU: knotsU,
    knotsV: knotsV,
    weights: weights,
  );
  @override
  Map<String, dynamic> toJson() => {
    't': 'spline',
    'du': degreeU,
    'dv': degreeV,
    'points': [for (final row in points) row.map(_v).toList()],
    'ku': knotsU,
    'kv': knotsV,
    'weights': weights,
  };
  factory BrepSplineSurface.fromJson(Map<String, dynamic> j) =>
      BrepSplineSurface(
        degreeU: j['du'] as int,
        degreeV: j['dv'] as int,
        points: [
          for (final row in j['points'] as List) (row as List).map(_p).toList(),
        ],
        knotsU: _splineDoubles(j['ku']),
        knotsV: _splineDoubles(j['kv']),
        weights: [for (final row in j['weights'] as List) _splineDoubles(row)],
      );
}

/// Reuses seed positions within one tessellation/operation. The surface itself
/// retains no cache, so caller edits to exposed control points cannot stale it.
class BrepSplineSurfaceLocator {
  BrepSplineSurfaceLocator(this.surface)
    : samples = surface._parameterSamples();
  final BrepSplineSurface surface;
  final List<(double, double, Vector3)> samples;
  (double, double) parameters(Vector3 point) =>
      surface._parameters(point, samples);
}
