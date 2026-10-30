part of 'geom.dart';

/// Exact linear interpolation between two oriented, bounded curves.
/// Both parameters range from zero to one. The representation stores the
/// curves, not the triangles used to display the patch.
class BrepRuledSurface implements BrepSurface {
  BrepRuledSurface(this.lower, this.upper);

  final BrepCurve lower;
  final BrepCurve upper;

  Vector3 _point(BrepCurve c, double u) => c.evaluate(c.t0 + (c.t1 - c.t0) * u);

  Vector3 _derivative(BrepCurve c, double u) {
    final dt = c.t1 - c.t0;
    final t = c.t0 + dt * u;
    if (c is BrepLine) return c.b - c.a;
    if (c is BrepSplineCurve) return c.derivative(t) * dt;
    if (c is BrepAffineCurve) return c.transform.vector(_derivative(c.base, u));
    if (c is BrepHelix) return c.derivative(t) * dt;
    if (c is BrepCircle) return c.tangent(t) * (c.radius * dt);
    if (c is BrepEllipse) {
      return (-c.xDir * (c.rMajor * math.sin(t)) +
              c.yDir * (c.rMinor * math.cos(t))) *
          dt;
    }
    const h = 1e-5;
    return (_point(c, u + h) - _point(c, u - h)) / (2 * h);
  }

  Vector3 du(double u, double v) =>
      _derivative(lower, u) * (1 - v) + _derivative(upper, u) * v;

  Vector3 dv(double u, double v) => _point(upper, u) - _point(lower, u);

  @override
  Vector3 evaluate(double u, double v) =>
      _point(lower, u) * (1 - v) + _point(upper, u) * v;

  @override
  Vector3 normalAt(double u, double v) => du(u, v).cross(dv(u, v)).normalized();

  @override
  (double, double) parameters(Vector3 p) {
    final a = lower, b = upper;
    if (a is BrepHelix && b is BrepHelix &&
        (a.axis-b.axis).length2 < 1e-24 && (a.xDir-b.xDir).length2 < 1e-24 &&
        (a.lead-b.lead).abs() < 1e-12 && (a.t0-b.t0).abs() < 1e-12 &&
        (a.t1-b.t1).abs() < 1e-12 && (a.t1-a.t0).abs() <= math.pi &&
        (b.origin-a.origin).cross(a.axis).length2 < 1e-20) {
      final d = p-a.origin;
      var angle = math.atan2(d.dot(a.yDir), d.dot(a.xDir));
      angle += 2*math.pi * (((a.t0+a.t1)/2-angle)/(2*math.pi)).round();
      final u = (angle-a.t0)/(a.t1-a.t0);
      final start = evaluate(u, 0), direction = dv(u, 0);
      final v = (p-start).dot(direction)/direction.length2;
      if (u >= -1e-9 && u <= 1+1e-9 && v >= -1e-9 && v <= 1+1e-9 &&
          (evaluate(u,v)-p).length2 < 1e-18) return (u,v);
    }
    var best = (0.0, 0.0);
    var bestD = double.infinity;
    // Multiple starts matter at a concave section. Each step is accepted only
    // when the true Euclidean residual decreases, including at patch borders.
    for (var seed = 0; seed <= 8; seed++) {
      var u = seed / 8;
      final a = _point(lower, u), d = _point(upper, u) - _point(lower, u);
      var v = d.length2 < 1e-24
          ? 0.5
          : ((p - a).dot(d) / d.length2).clamp(0.0, 1.0);
      for (var iter = 0; iter < 24; iter++) {
        final q = evaluate(u, v), r = q - p;
        final x = du(u, v), y = dv(u, v);
        final xx = x.length2, xy = x.dot(y), yy = y.length2;
        final det = xx * yy - xy * xy;
        if (det <= 1e-20 * xx * yy) break;
        final su = (yy * x.dot(r) - xy * y.dot(r)) / det;
        final sv = (xx * y.dot(r) - xy * x.dot(r)) / det;
        var step = 1.0;
        var moved = false;
        while (step >= 1 / 1024) {
          final nu = (u - su * step).clamp(0.0, 1.0);
          final nv = (v - sv * step).clamp(0.0, 1.0);
          if ((evaluate(nu, nv) - p).length2 < r.length2) {
            moved = (nu - u).abs() + (nv - v).abs() > 1e-12;
            u = nu;
            v = nv;
            break;
          }
          step *= 0.5;
        }
        if (!moved) break;
      }
      final distance = (evaluate(u, v) - p).length2;
      if (distance < bestD) {
        bestD = distance;
        best = (u, v);
      }
    }
    return best;
  }

  @override
  Vector3 closestPoint(Vector3 p) {
    final uv = parameters(p);
    return evaluate(uv.$1, uv.$2);
  }

  @override
  Map<String, dynamic> toJson() => {
    't': 'ruled',
    'a': lower.toJson(),
    'b': upper.toJson(),
  };

  factory BrepRuledSurface.fromJson(Map<String, dynamic> json) =>
      BrepRuledSurface(
        BrepCurve.fromJson(json['a'] as Map<String, dynamic>),
        BrepCurve.fromJson(json['b'] as Map<String, dynamic>),
      );
}
