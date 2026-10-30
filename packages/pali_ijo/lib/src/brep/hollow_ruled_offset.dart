part of 'hollow.dart';

/// Offsets a translational sweep in its normal plane, retaining a native
/// ruled/spline surface. This covers curved text/SVG walls and rational sweeps,
/// including affine placements; it does not replace the body by display mesh.
BrepSurface? _offsetExtrusion(BrepFace face, double thickness) {
  BrepRuledSurface? ruled(BrepSurface s) {
    if (s is BrepRuledSurface) return s;
    if (s is BrepSplineSurface) return s.ruled;
    if (s is BrepAffineSurface) {
      final base = ruled(s.base);
      if (base != null) {
        return BrepRuledSurface(
          BrepAffineCurve(base.lower, s.transform),
          BrepAffineCurve(base.upper, s.transform),
        );
      }
    }
    return null;
  }

  final s = ruled(face.surface);
  if (s == null) return null;
  // Reflections reverse a ruled parameterization's handedness, while an
  // affine surface deliberately preserves its original outward normal.
  double parity(BrepSurface surface) => surface is BrepAffineSurface
      ? surface.transform.determinant.sign * parity(surface.base)
      : 1;
  final normalParity = parity(face.surface);
  final outwardSign = (face.reversed ? -1.0 : 1.0) * normalParity;
  List<Vector3>? hull(BrepCurve c) {
    if (c is BrepLine) return [c.a, c.b];
    if (c is BrepSplineCurve) return c.points;
    if (c is BrepAffineCurve) {
      return hull(c.base)?.map(c.transform.point).toList();
    }
    return null;
  }

  final lowerHull = hull(s.lower), upperHull = hull(s.upper);
  if (lowerHull != null && upperHull != null) {
    final origin = s.evaluate(.5, .5), n = s.normalAt(.5, .5);
    if (n.length2 > .99 &&
        [
          ...lowerHull,
          ...upperHull,
        ].every((p) => (p - origin).dot(n).abs() < 1e-9)) {
      return BrepPlane(
        origin: origin + n * (-outwardSign * thickness),
        normal: n * normalParity,
      );
    }
  }
  final delta = s.evaluate(0, 1) - s.evaluate(0, 0);
  if (delta.length < 1e-8) {
    throw BrepUnsupported('Hollow sweep has a collapsed ruling');
  }
  // Compare rational control nets, not just samples: a twisted rail can
  // agree at sample points while departing between them.
  if (!sameRationalCurve(
    s.lower,
    s.upper,
    translation: delta,
    tolerance: 1e-9,
  )) {
    return null;
  }
  final c = s.lower, span = c.t1 - c.t0;
  Vector3 derivative(BrepCurve c, double t) => switch (c) {
    BrepSplineCurve() => c.derivative(t),
    BrepAffineCurve() => c.transform.vector(derivative(c.base, t)),
    BrepLine() => c.b - c.a,
    BrepCircle() => c.tangent(t) * c.radius,
    BrepEllipse() =>
      -c.xDir * (c.rMajor * math.sin(t)) + c.yDir * (c.rMinor * math.cos(t)),
    BrepHelix() => c.derivative(t),
    _ =>
      (c.evaluate(t + span * 1e-6) - c.evaluate(t - span * 1e-6)) /
          (span * 2e-6),
  };
  Vector3 tangent(double u) => derivative(c, c.t0 + span * u) * span;
  Vector3 normal(double u) {
    final n = tangent(u).cross(delta);
    if (n.length < 1e-12) {
      throw BrepUnsupported('Hollow sweep has a singular normal');
    }
    return n.normalized();
  }

  final sign = -outwardSign;
  Vector3 offset(double u) => s.evaluate(u, 0) + normal(u) * (sign * thickness);
  final start = offset(0), end = offset(1);
  final d0 = tangent(0), d1 = tangent(1);
  // Extend the supporting wall tangentially at sharp trimmed corners, so
  // adjacent offset walls can form a miter instead of being clamped apart.
  var railLength = 0.0;
  for (var i = 1; i <= 32; i++) {
    railLength += (s.evaluate(i / 32, 0) - s.evaluate((i - 1) / 32, 0)).length;
  }
  final extension = math.max(16 * thickness, railLength);
  final periodic =
      (s.evaluate(0, 0) - s.evaluate(1, 0)).length < 1e-8 &&
      normal(0).dot(normal(1)) > 1 - 1e-8;
  final pad0 = periodic ? 0.0 : extension / math.max(d0.length, 1e-12);
  final pad1 = periodic ? 0.0 : extension / math.max(d1.length, 1e-12);
  Vector3 point(double u) => u < 0
      ? start + d0 * u
      : u > 1
      ? end + d1 * (u - 1)
      : offset(u);
  Vector3 slope(double u, double lo, double hi) {
    final h = math.min(1e-4, (hi - lo) * 1e-3);
    if (u - h < lo) {
      return (point(u) * -3 + point(u + h) * 4 - point(u + 2 * h)) / (2 * h);
    }
    if (u + h > hi) {
      return (point(u) * 3 - point(u - h) * 4 + point(u - 2 * h)) / (2 * h);
    }
    return (point(u + h) - point(u - h)) / (2 * h);
  }

  final controls = <Vector3>[], knots = <double>[];
  const tolerance = 1e-7;
  void fit(double a, double b, double lo, double hi, int depth) {
    final p0 = point(a), p3 = point(b);
    final p1 = p0 + slope(a, lo, hi) * ((b - a) / 3);
    final p2 = p3 - slope(b, lo, hi) * ((b - a) / 3);
    var error = 0.0;
    for (final f in [.125, .25, .375, .5, .625, .75, .875]) {
      final q =
          p0 * math.pow(1 - f, 3).toDouble() +
          p1 * (3 * f * (1 - f) * (1 - f)) +
          p2 * (3 * f * f * (1 - f)) +
          p3 * (f * f * f);
      final u = a + (b - a) * f;
      error = math.max(error, (q - point(u)).length);
      if (u > 0 && u < 1 && slope(u, 0, 1).dot(tangent(u)) <= 0) {
        throw BrepUnsupported('Hollow thickness folds a curved inner wall');
      }
    }
    if (error > tolerance) {
      if (depth >= 18) {
        throw BrepUnsupported('Hollow curved-wall offset did not converge');
      }
      final m = (a + b) / 2;
      fit(a, m, lo, hi, depth + 1);
      fit(m, b, lo, hi, depth + 1);
      return;
    }
    if (controls.isEmpty) {
      controls.addAll([p0, p1, p2, p3]);
      knots.addAll(List.filled(4, a));
    } else {
      controls.addAll([p1, p2, p3]);
    }
    knots.addAll(List.filled(3, b));
  }

  if (pad0 > 0) fit(-pad0, 0, -pad0, 0, 0);
  Iterable<double> curveKnots(BrepCurve c) {
    if (c is BrepSplineCurve) return c.knots;
    if (c is BrepAffineCurve) return curveKnots(c.base);
    return const [];
  }

  final cuts = <double>{
    0,
    1,
    for (final knot in curveKnots(c))
      if ((knot - c.t0) / span > 0 && (knot - c.t0) / span < 1)
        (knot - c.t0) / span,
  }.toList()..sort();
  for (var i = 1; i < cuts.length; i++) {
    fit(cuts[i - 1], cuts[i], cuts[i - 1], cuts[i], 0);
  }
  if (pad1 > 0) fit(1, 1 + pad1, 1, 1 + pad1, 0);
  knots.add(knots.last);
  final lower = BrepSplineCurve(degree: 3, points: controls, knots: knots);
  final upper = BrepSplineCurve(
    degree: 3,
    points: [for (final p in controls) p + delta],
    knots: knots,
  );
  // Keep the surface normal convention used by the source face.
  final result = BrepRuledSurface(lower, upper);
  if (normalParity < 0) {
    return BrepRuledSurface(upper, lower);
  }
  return result;
}
