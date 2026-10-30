import 'geom.dart';
import 'solid.dart';
import 'spline_convert.dart';
import 'rational_ruled.dart';
import 'package:vector_math/vector_math_64.dart';

/// Compare exact rational spans, including reversed/affine parameterizations.
/// Optional translation/projection is used to recognize coincident sweeps.
bool sameRationalCurve(
  BrepCurve a,
  BrepCurve b, {
  Vector3? translation,
  Vector3? ignoreAxis,
  bool reversed = false,
  double tolerance = 1e-7,
}) {
  bool supported(BrepCurve c) =>
      c is BrepLine ||
      c is BrepCircle ||
      c is BrepEllipse ||
      c is BrepSplineCurve ||
      c is BrepAffineCurve && supported(c.base);
  if (!supported(a) || !supported(b)) return false;
  final x = BrepSplineConvert.curve(a), y = BrepSplineConvert.curve(b);
  if (x.degree != y.degree) return false;
  final cuts = <double>{0, 1};
  for (final c in [x, y]) {
    for (final knot in c.knots) {
      var t = (knot - c.t0) / (c.t1 - c.t0);
      if (identical(c, y) && reversed) t = 1 - t;
      if (t > 0 && t < 1) cuts.add(t);
    }
  }
  final sorted = cuts.toList()..sort();
  try {
    for (var i = 1; i < sorted.length; i++) {
      final p = BrepRationalRuled.homogeneousSpan(x, sorted[i - 1], sorted[i]);
      final q = BrepRationalRuled.homogeneousSpan(
        y,
        reversed ? 1 - sorted[i - 1] : sorted[i - 1],
        reversed ? 1 - sorted[i] : sorted[i],
      );
      if (p.length != q.length) return false;
      final ratio = p.first.w / q.first.w;
      for (var j = 0; j < p.length; j++) {
        if ((p[j].w - q[j].w * ratio).abs() > 1e-10 * p[j].w.abs()) {
          return false;
        }
        var delta =
            Vector3(p[j].x, p[j].y, p[j].z) / p[j].w +
            (translation ?? Vector3.zero()) -
            Vector3(q[j].x, q[j].y, q[j].z) / q[j].w;
        if (ignoreAxis != null) delta -= ignoreAxis * delta.dot(ignoreAxis);
        if (delta.length > tolerance) return false;
      }
    }
  } on BrepUnsupported {
    return false;
  }
  return true;
}
