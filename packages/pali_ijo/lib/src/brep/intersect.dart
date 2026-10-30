import 'dart:math' as math;

import 'package:vector_math/vector_math_64.dart';

import 'geom.dart';
import 'solid.dart';
import 'spline_intersect.dart';

/// Surface–surface intersection curves clipped to neither face yet.
class BrepIntersect {
  BrepIntersect._();

  static List<BrepCurve> surfaces(BrepSurface a, BrepSurface b) {
    if (a is BrepSplineSurface) a = a.ruled ?? a;
    if (b is BrepSplineSurface) b = b.ruled ?? b;
    if (a is BrepRuledSurface && BrepQuadric.supports(b))
      return ruledQuadricBranches(a, b);
    if (b is BrepRuledSurface && BrepQuadric.supports(a))
      return ruledQuadricBranches(b, a);
    if (a is BrepPlane && b is BrepPlane) {
      return _planePlane(a, b);
    }
    if (a is BrepPlane && b is BrepSphere) {
      return _planeSphere(a, b);
    }
    if (a is BrepSphere && b is BrepPlane) {
      return _planeSphere(b, a);
    }
    if (a is BrepPlane && b is BrepCylinder) {
      return _planeCylinder(a, b);
    }
    if (a is BrepCylinder && b is BrepPlane) {
      return _planeCylinder(b, a);
    }
    if (a is BrepSphere && b is BrepSphere) {
      return _sphereSphere(a, b);
    }
    if (a is BrepPlane && b is BrepTorus) {
      return _planeTorusPrincipal(a, b);
    }
    if (a is BrepTorus && b is BrepPlane) {
      return _planeTorusPrincipal(b, a);
    }
    if (a is BrepCylinder && b is BrepSphere) {
      return _cylinderSphere(a, b);
    }
    if (a is BrepSphere && b is BrepCylinder) {
      return _cylinderSphere(b, a);
    }
    if (a is BrepCylinder && b is BrepCylinder) {
      return _cylinderCylinder(a, b);
    }
    if (a is BrepCylinder && b is BrepCone) return _coaxialCylinderCone(a, b);
    if (b is BrepCylinder && a is BrepCone) return _coaxialCylinderCone(b, a);
    if (a is BrepPlane && b is BrepCone) {
      return _planeCone(a, b);
    }
    if (a is BrepCone && b is BrepPlane) {
      return _planeCone(b, a);
    }
    if (a is BrepSplineSurface ||
        b is BrepSplineSurface ||
        a is BrepRuledSurface && b is BrepRuledSurface)
      return BrepSplineIntersect.intersect(a, b);
    throw BrepUnsupported(
      'Unsupported surface intersection (${a.runtimeType} × ${b.runtimeType})',
    );
  }

  /// Real branches for cylinder/quadric intersections, with all angular
  /// discriminant roots isolated before evaluating a branch. A quartic in
  /// tan(u/2) finds small loops that a fixed angular sampling can miss.
  static List<BrepIntersectionCurve> quadricBranches(
    BrepCylinder carrier,
    BrepSurface target,
  ) {
    final probe = BrepIntersectionCurve(
      carrier: carrier,
      target: target,
      branch: 1,
    );
    final f0 = probe.discriminant(0), fp = probe.discriminant(math.pi);
    final fh = probe.discriminant(math.pi / 2),
        f3h = probe.discriminant(3 * math.pi / 2);
    final c0 = (f0 + fp + fh + f3h) / 4;
    final c1 = (f0 - fp) / 2, c2 = (fh - f3h) / 2;
    final c3 = (f0 + fp) / 2 - c0;
    final c4 = probe.discriminant(math.pi / 4) - c0 - (c1 + c2) / math.sqrt2;
    final roots = _realPolynomialRoots([
      c0 + c1 + c3,
      2 * c2 + 4 * c4,
      2 * c0 - 6 * c3,
      2 * c2 - 4 * c4,
      c0 - c1 + c3,
    ]);
    final cuts = <double>[0, math.pi * 2];
    for (final root in roots) {
      var u = 2 * math.atan(root);
      if (u < 0) u += math.pi * 2;
      if (u > 1e-10 && u < math.pi * 2 - 1e-10) cuts.add(u);
    }
    final scale = math.max(
      1.0,
      [f0, fp, fh, f3h].map((v) => v.abs()).reduce(math.max),
    );
    if (fp.abs() < scale * 1e-10) cuts.add(math.pi);
    cuts.sort();
    final unique = <double>[];
    for (final u in cuts) {
      if (unique.isEmpty || u - unique.last > 1e-9) unique.add(u);
    }
    final out = <BrepIntersectionCurve>[];
    for (var i = 1; i < unique.length; i++) {
      final lo = unique[i - 1], hi = unique[i];
      if (probe.discriminant((lo + hi) / 2) <= scale * 1e-12) continue;
      for (final branch in [-1, 1]) {
        out.add(
          BrepIntersectionCurve(
            carrier: carrier,
            target: target,
            branch: branch,
            t0: lo,
            t1: hi,
          ),
        );
      }
    }
    return out;
  }

  /// Ray hits of a bounded ruled wall. Eliminate its ruling coordinate using
  /// two planes through the ray, then solve only the remaining scalar equation.
  /// Null means a coincident ruling requires a different ray direction.
  static List<double>? rayRuled(
    BrepRuledSurface s,
    Vector3 origin,
    Vector3 ray,
  ) {
    final n1 = brepFrameU(ray), n2 = ray.cross(n1);
    double determinant(double u) {
      final a = s.evaluate(u, 0) - origin, d = s.dv(u, 0);
      return a.dot(n1) * d.dot(n2) - a.dot(n2) * d.dot(n1);
    }

    final probes = [for (var i = 0; i <= 8; i++) determinant(i / 8)];
    final scale = math.max(1.0, probes.map((x) => x.abs()).reduce(math.max));
    if (probes.every((v) => v.abs() < scale * 1e-13)) return null;
    final out = <double>[];
    for (final u in unitRoots(determinant)) {
      final a = s.evaluate(u, 0), d = s.dv(u, 0), h = a - origin;
      final x = d.dot(n1), y = d.dot(n2);
      if (math.max(x.abs(), y.abs()) < 1e-12 * d.length) {
        if ((h - ray * h.dot(ray)).length < 1e-7) return null;
        continue;
      }
      final v = x.abs() > y.abs() ? -h.dot(n1) / x : -h.dot(n2) / y;
      if (v < -1e-9 || v > 1 + 1e-9) continue;
      final p = s.evaluate(u, v), t = (p - origin).dot(ray);
      if ((p - origin - ray * t).length > 1e-6) return null;
      out.add(t);
    }
    return out;
  }

  /// Ruled/quadric branches split at discriminant roots and patch rims.
  /// Each retained segment evaluates the original implicit equation.
  static List<BrepCurve> ruledQuadricBranches(
    BrepRuledSurface s,
    BrepSurface target,
  ) {
    final probe = BrepRuledQuadricCurve(surface: s, target: target, branch: 1);
    final cuts = <double>[0, 1];
    cuts.addAll(unitRoots(probe.discriminant));
    cuts.addAll(unitRoots((u) => probe.coefficients(u).$1));
    cuts.addAll(unitRoots((u) => BrepQuadric.value(target, s.evaluate(u, 0))));
    cuts.addAll(unitRoots((u) => BrepQuadric.value(target, s.evaluate(u, 1))));
    cuts.sort();
    final ts = <double>[];
    for (final t in cuts) {
      if (ts.isEmpty || t - ts.last > 1e-10) ts.add(t);
    }
    final out = <BrepCurve>[];
    for (var i = 1; i < ts.length; i++) {
      final lo = ts[i - 1], hi = ts[i], mid = (lo + hi) / 2;
      final (a, b, c) = probe.coefficients(mid);
      final scale = math.max(
        1.0,
        math.max(a.abs(), math.max(b.abs(), c.abs())),
      );
      if (probe.discriminant(mid) < -1e-12 * scale * scale) continue;
      if (a.abs() < 1e-14 * scale && b.abs() < 1e-14 * scale) {
        if (c.abs() < 1e-10 * scale)
          throw BrepUnsupported('Coincident ruled/quadric patch');
        continue;
      }
      for (final branch in a.abs() < 1e-14 * scale ? [1] : [-1, 1]) {
        final curve = BrepRuledQuadricCurve(
          surface: s,
          target: target,
          branch: branch,
          t0: lo,
          t1: hi,
        );
        final v = curve.ruling(mid);
        if (v < -1e-8 || v > 1 + 1e-8) continue;
        final endpoints = [curve.ruling(lo), curve.ruling(hi)];
        if (endpoints.any((v) => !v.isFinite || v < -1e-5 || v > 1 + 1e-5)) {
          throw BrepUnsupported(
            'Ruled intersection branch reaches an unresolved singularity',
          );
        }
        out.add(curve);
      }
    }
    return out;
  }

  /// Brackets both crossings and stationary zeros. Stationary-point refinement
  /// finds narrow real branches that endpoint sign tests alone can miss.
  static List<double> unitRoots(double Function(double) f) {
    const n = 256;
    final values = [for (var i = 0; i <= n; i++) f(i / n)];
    if (values.any((v) => !v.isFinite))
      throw BrepUnsupported(
        'Intersection equation is not finite on its domain',
      );
    final scale = math.max(1.0, values.map((v) => v.abs()).reduce(math.max));
    if (values.every((v) => v.abs() < 1e-13 * scale)) return [];
    final knots = <double>[0, 1];
    for (var i = 1; i < n; i++) {
      if ((values[i] - values[i - 1]) * (values[i + 1] - values[i]) > 0)
        continue;
      final minimum = values[i] <= values[i - 1] && values[i] <= values[i + 1];
      var a = (i - 1) / n, b = (i + 1) / n;
      for (var k = 0; k < 48; k++) {
        final x = a + (b - a) / 3, y = b - (b - a) / 3;
        if ((f(x) < f(y)) == minimum) {
          b = y;
        } else {
          a = x;
        }
      }
      knots.add((a + b) / 2);
    }
    knots.addAll([for (var i = 1; i < n; i++) i / n]);
    knots.sort();
    final roots = <double>[];
    for (var i = 0; i < knots.length; i++) {
      final x = knots[i], fx = f(x);
      if (fx.abs() < 1e-12 * scale) roots.add(x);
      if (i == 0) continue;
      var a = knots[i - 1], b = x, fa = f(a);
      if (fa * fx >= 0) continue;
      for (var k = 0; k < 52; k++) {
        final m = (a + b) / 2, fm = f(m);
        if (fa * fm <= 0) {
          b = m;
        } else {
          a = m;
          fa = fm;
        }
      }
      roots.add((a + b) / 2);
    }
    roots.sort();
    final unique = <double>[];
    for (final r in roots) {
      if (unique.isEmpty || r - unique.last > 1e-9) unique.add(r);
    }
    return unique;
  }

  /// Intersections of a line with the analytic torus. Work near the torus
  /// center in radius-sized units to keep the quartic well conditioned.
  static List<double> rayTorus(
    BrepTorus torus,
    Vector3 origin,
    Vector3 direction,
  ) {
    final length = direction.length;
    if (length < 1e-15) return [];
    final axis = direction / length;
    final shift = -(origin - torus.center).dot(axis);
    final scale = torus.majorRadius + torus.minorRadius;
    final p = (origin + axis * shift - torus.center) / scale;
    if (p.length2 > 1 + 1e-12) return [];
    final minor = torus.minorRadius / scale, major = torus.majorRadius / scale;
    final z = p.dot(torus.axis), dz = axis.dot(torus.axis);
    final g0 = p.length2 + major * major - minor * minor, g1 = 2 * p.dot(axis);
    return [
      for (final root in _realPolynomialRoots([
        g0 * g0 - 4 * major * major * (p.length2 - z * z),
        2 * g0 * g1 - 8 * major * major * (p.dot(axis) - z * dz),
        g1 * g1 + 2 * g0 - 4 * major * major * (1 - dz * dz),
        2 * g1,
        1,
      ]))
        (shift + root * scale) / length,
    ];
  }

  static List<double> _realPolynomialRoots(List<double> coefficients) {
    final c = List<double>.from(coefficients);
    final scale = c.map((v) => v.abs()).fold(0.0, math.max);
    if (scale == 0) return [];
    for (var i = 0; i < c.length; i++) {
      c[i] /= scale;
    }
    while (c.length > 1 && c.last.abs() < 1e-13) {
      c.removeLast();
    }
    if (c.length == 1) return [];
    if (c.length == 2) return [-c[0] / c[1]];
    double value(double x) {
      var result = c.last;
      for (var i = c.length - 2; i >= 0; i--) {
        result = result * x + c[i];
      }
      return result;
    }

    final bound =
        1 +
        c.take(c.length - 1).map((v) => (v / c.last).abs()).reduce(math.max);
    final stationary = _realPolynomialRoots([
      for (var i = 1; i < c.length; i++) i * c[i],
    ]).where((x) => x > -bound && x < bound).toList()..sort();
    final intervals = [-bound, ...stationary, bound];
    final roots = <double>[];
    for (final x in stationary) {
      final error = c.asMap().entries.fold(
        0.0,
        (double s, e) => s + e.value.abs() * math.pow(x.abs(), e.key),
      );
      if (value(x).abs() < math.max(1.0, error) * 1e-12) roots.add(x);
    }
    for (var i = 1; i < intervals.length; i++) {
      var lo = intervals[i - 1], hi = intervals[i];
      var flo = value(lo), fhi = value(hi);
      if (flo * fhi >= 0) continue;
      for (var k = 0; k < 100; k++) {
        final mid = (lo + hi) / 2, fm = value(mid);
        if (flo * fm <= 0) {
          hi = mid;
          fhi = fm;
        } else {
          lo = mid;
          flo = fm;
        }
      }
      roots.add((lo + hi) / 2);
    }
    return roots..sort();
  }

  static List<BrepCurve> _planePlane(BrepPlane a, BrepPlane b) {
    var dir = a.normal.cross(b.normal);
    if (dir.length < 1e-10) return const [];
    dir.normalize();
    // Point on both planes: solve in the plane spanned by nA, nB.
    final n1 = a.normal, n2 = b.normal;
    final d1 = a.origin.dot(n1), d2 = b.origin.dot(n2);
    final n1n1 = n1.dot(n1), n2n2 = n2.dot(n2), n1n2 = n1.dot(n2);
    final det = n1n1 * n2n2 - n1n2 * n1n2;
    if (det.abs() < 1e-14) return const [];
    final c1 = (d1 * n2n2 - d2 * n1n2) / det;
    final c2 = (d2 * n1n1 - d1 * n1n2) / det;
    final p = n1 * c1 + n2 * c2;
    return [BrepLine(p - dir * 1e4, p + dir * 1e4)];
  }

  static List<BrepCurve> _planeSphere(BrepPlane a, BrepSphere b) {
    final h = a.height(b.center);
    if (h.abs() > b.radius + 1e-9) return const [];
    final r = math.sqrt(math.max(0, b.radius * b.radius - h * h));
    if (r < 1e-9) return const [];
    final c = b.center - a.normal * h;
    return [BrepCircle(center: c, normal: a.normal, radius: r, xDir: a.xDir)];
  }

  static List<BrepCurve> _planeCylinder(BrepPlane a, BrepCylinder b) {
    final an = a.normal.dot(b.axis).abs();
    if (an > 0.999) {
      final h = a.height(b.origin);
      final along = -h / a.normal.dot(b.axis);
      final c = b.origin + b.axis * along;
      if ((c - a.closestPoint(c)).length > 1e-6) {
        // Plane ⟂ axis: circle of radius r.
      }
      return [
        BrepCircle(
          center: a.closestPoint(
            b.origin + b.axis * ((a.origin - b.origin).dot(b.axis)),
          ),
          normal: a.normal,
          radius: b.radius,
          xDir: b.xDir,
        ),
      ];
    }
    if (an < 1e-3) {
      final h = a.height(b.origin);
      final d = h.abs();
      if (d > b.radius + 1e-9) return const [];
      if ((b.radius - d).abs() < 1e-6) {
        // Tangent generator: project the axis onto the plane (not the
        // closest point to the plane origin, which misses the contact).
        final proj = b.origin - a.normal * h;
        final radial = proj - b.origin;
        final x = radial.length2 < 1e-20
            ? b.origin + a.normal * (h > 0 ? -b.radius : b.radius)
            : b.origin + radial.normalized() * b.radius;
        return [BrepLine(x - b.axis * 1e4, x + b.axis * 1e4)];
      }
      final t = math.sqrt(b.radius * b.radius - d * d);
      final inward = a.normal * (h > 0 ? -1.0 : 1.0);
      var lat = inward.cross(b.axis);
      if (lat.length < 1e-8) lat = b.xDir.clone();
      lat.normalize();
      final mid = b.origin + a.normal * (-h);
      // The signed projection already chooses the correct side of the axis.
      // Normalizing it and multiplying by -h again reflects positive-height
      // intersections off the cutting plane.
      final p1 = mid + lat * t;
      final p2 = mid - lat * t;
      return [
        BrepLine(p1 - b.axis * 1e4, p1 + b.axis * 1e4),
        BrepLine(p2 - b.axis * 1e4, p2 + b.axis * 1e4),
      ];
    }
    return [_planeCylinderEllipse(a, b)];
  }

  static BrepEllipse _planeCylinderEllipse(BrepPlane plane, BrepCylinder cyl) {
    final n = plane.normal;
    final denom = n.dot(cyl.axis);
    final along = (plane.origin - cyl.origin).dot(n) / denom;
    final center = cyl.origin + cyl.axis * along;
    final cosA = denom.abs().clamp(1e-9, 1.0);
    var majorDir = cyl.axis - n * cyl.axis.dot(n);
    if (majorDir.length < 1e-8) {
      majorDir = plane.xDir;
    }
    majorDir.normalize();
    return BrepEllipse(
      center: center,
      normal: n,
      xDir: majorDir,
      rMajor: cyl.radius / cosA,
      rMinor: cyl.radius,
    );
  }

  static List<BrepCurve> _sphereSphere(BrepSphere a, BrepSphere b) {
    final d = b.center - a.center;
    final dist = d.length;
    if (dist < 1e-9 || dist > a.radius + b.radius + 1e-9) return const [];
    if (dist + math.min(a.radius, b.radius) <
        math.max(a.radius, b.radius) - 1e-9) {
      return const [];
    }
    final n = d.normalized();
    final x =
        (dist * dist + a.radius * a.radius - b.radius * b.radius) / (2 * dist);
    final r = math.sqrt(math.max(0, a.radius * a.radius - x * x));
    if (r < 1e-9) return const [];
    final c = a.center + n * x;
    return [BrepCircle(center: c, normal: n, radius: r)];
  }

  static List<BrepCurve> _planeTorusPrincipal(BrepPlane a, BrepTorus b) {
    final an = a.normal.dot(b.axis).abs();
    final h = a.height(b.center);
    if (an > 0.999) {
      final rho2 = b.minorRadius * b.minorRadius - h * h;
      if (rho2 < -1e-9) return const [];
      final rho = math.sqrt(math.max(0, rho2));
      final c = b.center - a.normal * h;
      var xDir = a.xDir.clone();
      if (xDir.cross(a.normal).length < 1e-8) {
        xDir = b.xDir.clone();
      }
      if (rho < 1e-9) {
        return [
          BrepCircle(
            center: c,
            normal: a.normal,
            radius: b.majorRadius,
            xDir: xDir,
          ),
        ];
      }
      return [
        BrepCircle(
          center: c,
          normal: a.normal,
          radius: b.majorRadius + rho,
          xDir: xDir,
        ),
        BrepCircle(
          center: c,
          normal: a.normal,
          radius: (b.majorRadius - rho).abs(),
          xDir: xDir,
        ),
      ];
    }
    // Plane contains the torus axis: two tube cross-section circles.
    if (an < 1e-3 && h.abs() <= 1e-6) {
      var dir = a.normal.cross(b.axis);
      if (dir.length < 1e-12) return const [];
      dir.normalize();
      var xDir = b.axis.clone();
      if (xDir.cross(a.normal).length < 1e-8) {
        xDir = dir.clone();
      }
      return [
        BrepCircle(
          center: b.center + dir * b.majorRadius,
          normal: a.normal,
          radius: b.minorRadius,
          xDir: xDir,
        ),
        BrepCircle(
          center: b.center - dir * b.majorRadius,
          normal: a.normal,
          radius: b.minorRadius,
          xDir: xDir,
        ),
      ];
    }
    if (!_planeMeetsTorus(a, b)) return const [];
    return BrepSplineIntersect.intersect(a, b);
  }

  static bool _planeMeetsTorus(BrepPlane pl, BrepTorus tor) {
    var pos = false;
    var neg = false;
    const n = 16;
    for (var i = 0; i < n; i++) {
      for (var j = 0; j < n; j++) {
        final h = pl.height(
          tor.evaluate(i / n * math.pi * 2, j / n * math.pi * 2),
        );
        if (h > 1e-4) pos = true;
        if (h < -1e-4) neg = true;
        if (pos && neg) return true;
      }
    }
    return false;
  }

  /// Coaxial intersections are circles; offset intersections retain exact
  /// quadric branches. A narrow cylinder passing through a larger sphere still
  /// has two rims even when its entire circular footprint is inside the sphere.
  static List<BrepCurve> _cylinderSphere(BrepCylinder cyl, BrepSphere sph) {
    final d = sph.center - cyl.origin;
    final along = d.dot(cyl.axis);
    final radial = d - cyl.axis * along;
    if (radial.length > 1e-6) return quadricBranches(cyl, sph);
    final h2 = sph.radius * sph.radius - cyl.radius * cyl.radius;
    if (h2 < -1e-9) return const [];
    final mid = cyl.origin + cyl.axis * along;
    if (h2 < 1e-12) {
      return [
        BrepCircle(
          center: mid,
          normal: cyl.axis,
          radius: cyl.radius,
          xDir: cyl.xDir,
        ),
      ];
    }
    final h = math.sqrt(h2);
    return [
      BrepCircle(
        center: mid + cyl.axis * h,
        normal: cyl.axis,
        radius: cyl.radius,
        xDir: cyl.xDir,
      ),
      BrepCircle(
        center: mid - cyl.axis * h,
        normal: cyl.axis,
        radius: cyl.radius,
        xDir: cyl.xDir,
      ),
    ];
  }

  static List<BrepCurve> _coaxialCylinderCone(
    BrepCylinder cylinder,
    BrepCone cone,
  ) {
    final d = cylinder.origin - cone.origin;
    if (cylinder.axis.cross(cone.axis).length <= 1e-8) {
      final distance = (d - cone.axis * d.dot(cone.axis)).length;
      final rmin = math.min(cone.radius, cone.topRadius),
          rmax = math.max(cone.radius, cone.topRadius);
      if (distance > cylinder.radius + rmax + 1e-8 ||
          cylinder.radius > distance + rmax + 1e-8 ||
          rmin > distance + cylinder.radius + 1e-8)
        return [];
    }
    if (cylinder.axis.cross(cone.axis).length > 1e-8 ||
        (d - cone.axis * d.dot(cone.axis)).length > 1e-7) {
      throw BrepUnsupported(
        'Non-coaxial cylinder-cone intersection needs a general quadric solver',
      );
    }
    final slope = (cone.topRadius - cone.radius) / cone.height;
    // Equal-radius coaxial surfaces are coincident or disjoint, with no
    // isolated intersection curve.
    if (slope.abs() < 1e-14) return [];
    final z = (cylinder.radius - cone.radius) / slope;
    if (z < -1e-8 || z > cone.height + 1e-8) return [];
    return [
      BrepCircle(
        center: cone.origin + cone.axis * z,
        normal: cone.axis,
        radius: cylinder.radius,
        xDir: cone.xDir,
      ),
    ];
  }

  static List<BrepCurve> _planeCone(BrepPlane a, BrepCone b) {
    final an = a.normal.dot(b.axis).abs();
    if (an > 0.999) {
      final den = a.normal.dot(b.axis);
      if (den.abs() < 1e-12) return const [];
      final z = -a.height(b.origin) / den;
      if (z < -1e-6 || z > b.height + 1e-6) return const [];
      final r = b.radiusAt(z);
      if (r < 1e-9) return const [];
      return [
        BrepCircle(
          center: b.origin + b.axis * z,
          normal: a.normal,
          radius: r,
          xDir: b.xDir,
        ),
      ];
    }
    // Finite-wall cut: every generator meets the plane once → an ellipse
    // (or a circle, already handled). A hyperbola misses some generators.
    final apex = b.origin + b.axis * b.height;
    const n = 32;
    final hits = <Vector3>[];
    var missed = false;
    var parallel = false;
    for (var i = 0; i < n; i++) {
      final base = b.evaluate(i / n * math.pi * 2, 0);
      final d = base - apex;
      final den = a.normal.dot(d);
      if (den.abs() < 1e-14) {
        parallel = true;
        continue;
      }
      final t = -a.height(apex) / den;
      if (t < -1e-4 || t > 1 + 1e-3) {
        missed = true;
        continue;
      }
      hits.add(apex + d * t);
    }
    if (hits.isEmpty) return const [];
    if (missed || parallel) {
      throw BrepUnsupported(
        'Cone–plane intersection (oblique) is not implemented yet',
      );
    }
    return [_ellipseThrough(hits, a.normal)];
  }

  /// Best-fit ellipse of coplanar [pts] (a true plane–cone cut is exact).
  static BrepEllipse _ellipseThrough(List<Vector3> pts, Vector3 normal) {
    var c = Vector3.zero();
    for (final p in pts) {
      c += p;
    }
    c /= pts.length.toDouble();
    var x = normal.cross(Vector3(0, 0, 1));
    if (x.length2 < 1e-12) x = normal.cross(Vector3(0, 1, 0));
    x.normalize();
    final y = normal.cross(x).normalized();
    var xx = 0.0, xy = 0.0, yy = 0.0;
    for (final p in pts) {
      final u = (p - c).dot(x);
      final v = (p - c).dot(y);
      xx += u * u;
      xy += u * v;
      yy += v * v;
    }
    final inv = 1.0 / pts.length;
    xx *= inv;
    xy *= inv;
    yy *= inv;
    final tr = xx + yy;
    final disc = math.sqrt(math.max(0, tr * tr / 4 - (xx * yy - xy * xy)));
    final l1 = tr / 2 + disc;
    Vector3 major;
    if (xy.abs() > 1e-12) {
      major = (x * (l1 - yy) + y * xy);
    } else {
      major = xx >= yy ? x : y;
    }
    if (major.length < 1e-12) major = x.clone();
    major.normalize();
    final minor = normal.cross(major).normalized();
    var rMaj = 0.0, rMin = 0.0;
    for (final p in pts) {
      rMaj = math.max(rMaj, (p - c).dot(major).abs());
      rMin = math.max(rMin, (p - c).dot(minor).abs());
    }
    return BrepEllipse(
      center: c,
      normal: normal,
      xDir: major,
      rMajor: math.max(rMaj, 1e-8),
      rMinor: math.max(rMin, 1e-8),
    );
  }

  /// Parallel cylinders → generator lines. Crossing axes are refused.
  static List<BrepCurve> _cylinderCylinder(BrepCylinder a, BrepCylinder b) {
    if (a.axis.cross(b.axis).length > 1e-6) {
      throw BrepUnsupported(
        'Skew or crossing cylinder–cylinder intersection is not implemented yet',
      );
    }
    final delta = b.origin - a.origin;
    final radial = delta - a.axis * delta.dot(a.axis);
    final dist = radial.length;
    final rSum = a.radius + b.radius;
    final rDiff = (a.radius - b.radius).abs();
    if (dist > rSum + 1e-9 || dist + 1e-9 < rDiff) return const [];
    if (dist < 1e-9) return const [];
    final n = radial.normalized();
    var lat = a.axis.cross(n);
    if (lat.length < 1e-8) lat = a.xDir;
    lat.normalize();
    final x =
        (dist * dist + a.radius * a.radius - b.radius * b.radius) / (2 * dist);
    final h2 = a.radius * a.radius - x * x;
    if (h2 < -1e-9) return const [];
    final h = math.sqrt(math.max(0, h2));
    final mid = a.origin + n * x;
    if (h < 1e-8) {
      return [BrepLine(mid - a.axis * 1e4, mid + a.axis * 1e4)];
    }
    return [
      BrepLine(mid + lat * h - a.axis * 1e4, mid + lat * h + a.axis * 1e4),
      BrepLine(mid - lat * h - a.axis * 1e4, mid - lat * h + a.axis * 1e4),
    ];
  }
}
