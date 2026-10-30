import 'dart:math' as math;

import 'package:vector_math/vector_math_64.dart';
import 'solid.dart' show BrepUnsupported;

part 'intersection_curve.dart';
part 'helix.dart';
part 'ruled_surface.dart';
part 'offset_surface.dart';
part 'affine_geometry.dart';
part 'spline_geometry.dart';
part 'surface_intersection_curve.dart';
part 'plane_section_curve.dart';
part 'ruled_quadric_curve.dart';

const brepEps = 1e-8;

Vector3 brepFrameU(Vector3 axis) {
  var u = axis.cross(Vector3(0, 0, 1));
  if (u.length < 1e-8) u = axis.cross(Vector3(0, 1, 0));
  return u.normalized();
}

/// Analytic surface. UV domains: plane unbounded; cylinder θ∈[0,2π), z along axis;
/// sphere θ∈[0,2π), φ∈[0,π]; torus u,v ∈[0,2π); cone like cylinder.
sealed class BrepSurface {
  Vector3 evaluate(double u, double v);
  Vector3 normalAt(double u, double v);
  (double, double)? parameters(Vector3 p);
  Vector3 closestPoint(Vector3 p);
  Map<String, dynamic> toJson();

  static BrepSurface fromJson(Map<String, dynamic> json) {
    switch (json['t'] as String) {
      case 'offset':
        return BrepOffsetSurface.fromJson(json);
      case 'spline':
        return BrepSplineSurface.fromJson(json);
      case 'plane':
        return BrepPlane.fromJson(json);
      case 'cyl':
        return BrepCylinder.fromJson(json);
      case 'sphere':
        return BrepSphere.fromJson(json);
      case 'torus':
        return BrepTorus.fromJson(json);
      case 'cone':
        return BrepCone.fromJson(json);
      case 'ruled':
        return BrepRuledSurface.fromJson(json);
      case 'affine':
        return BrepAffineSurface.fromJson(json);
      default:
        throw FormatException('Unknown surface ${json['t']}');
    }
  }
}

class BrepPlane implements BrepSurface {
  BrepPlane({required this.origin, required Vector3 normal, Vector3? xDir})
    : normal = normal.normalized(),
      xDir = _planeX(normal.normalized(), xDir),
      yDir = normal
          .normalized()
          .cross(_planeX(normal.normalized(), xDir))
          .normalized();

  final Vector3 origin;
  final Vector3 normal;
  final Vector3 xDir;
  final Vector3 yDir;

  static Vector3 _planeX(Vector3 n, Vector3? x) {
    if (x != null) {
      final p = x - n * x.dot(n);
      if (p.length > 1e-10) return p.normalized();
    }
    return brepFrameU(n);
  }

  double height(Vector3 p) => (p - origin).dot(normal);

  @override
  Vector3 evaluate(double u, double v) => origin + xDir * u + yDir * v;

  @override
  Vector3 normalAt(double u, double v) => normal.clone();

  @override
  (double, double)? parameters(Vector3 p) {
    final d = p - origin;
    return (d.dot(xDir), d.dot(yDir));
  }

  @override
  Vector3 closestPoint(Vector3 p) => p - normal * height(p);

  BrepPlane offset(double d) =>
      BrepPlane(origin: origin + normal * d, normal: normal, xDir: xDir);

  @override
  Map<String, dynamic> toJson() => {
    't': 'plane',
    'o': _v(origin),
    'n': _v(normal),
    'x': _v(xDir),
  };

  factory BrepPlane.fromJson(Map<String, dynamic> json) => BrepPlane(
    origin: _p(json['o']),
    normal: _p(json['n']),
    xDir: _p(json['x']),
  );
}

class BrepCylinder implements BrepSurface {
  BrepCylinder({
    required this.origin,
    required Vector3 axis,
    required this.radius,
    Vector3? xDir,
  }) : axis = axis.normalized(),
       xDir = _planeX(axis.normalized(), xDir) {
    yDir = this.axis.cross(this.xDir).normalized();
  }

  final Vector3 origin;
  final Vector3 axis;
  final double radius;
  final Vector3 xDir;
  late final Vector3 yDir;

  static Vector3 _planeX(Vector3 n, Vector3? x) {
    if (x != null) {
      final p = x - n * x.dot(n);
      if (p.length > 1e-10) return p.normalized();
    }
    return brepFrameU(n);
  }

  @override
  Vector3 evaluate(double u, double v) =>
      origin +
      axis * v +
      xDir * (radius * math.cos(u)) +
      yDir * (radius * math.sin(u));

  @override
  Vector3 normalAt(double u, double v) =>
      (xDir * math.cos(u) + yDir * math.sin(u)).normalized();

  @override
  (double, double)? parameters(Vector3 p) {
    final d = p - origin;
    final v = d.dot(axis);
    final rad = d - axis * v;
    if (rad.length < 1e-12) return (0.0, v);
    return (math.atan2(rad.dot(yDir), rad.dot(xDir)), v);
  }

  @override
  Vector3 closestPoint(Vector3 p) {
    final d = p - origin;
    final v = d.dot(axis);
    var rad = d - axis * v;
    if (rad.length < 1e-12) rad = xDir.clone();
    rad.normalize();
    return origin + axis * v + rad * radius;
  }

  BrepCylinder offset(double d) =>
      BrepCylinder(origin: origin, axis: axis, radius: radius + d, xDir: xDir);

  @override
  Map<String, dynamic> toJson() => {
    't': 'cyl',
    'o': _v(origin),
    'a': _v(axis),
    'r': radius,
    'x': _v(xDir),
  };

  factory BrepCylinder.fromJson(Map<String, dynamic> json) => BrepCylinder(
    origin: _p(json['o']),
    axis: _p(json['a']),
    radius: (json['r'] as num).toDouble(),
    xDir: _p(json['x']),
  );
}

class BrepSphere implements BrepSurface {
  BrepSphere({required this.center, required this.radius});

  final Vector3 center;
  final double radius;

  @override
  Vector3 evaluate(double u, double v) {
    final s = math.sin(v);
    return center +
        Vector3(
          radius * s * math.cos(u),
          radius * s * math.sin(u),
          radius * math.cos(v),
        );
  }

  @override
  Vector3 normalAt(double u, double v) =>
      (evaluate(u, v) - center).normalized();

  @override
  (double, double)? parameters(Vector3 p) {
    final d = p - center;
    if (d.length < 1e-12) return (0.0, 0.0);
    return (math.atan2(d.y, d.x), math.acos((d.z / d.length).clamp(-1.0, 1.0)));
  }

  @override
  Vector3 closestPoint(Vector3 p) {
    var d = p - center;
    if (d.length < 1e-12) return center + Vector3(0, 0, radius);
    return center + d.normalized() * radius;
  }

  /// Pole of the spherical cap bounded by a closed [curve] when that
  /// circle is walked CCW around ±[BrepCircle.normal]. [coedgeReversed]
  /// flips the sense. The filled cap is the one on the +n side of the
  /// circle's plane.
  Vector3 circularCapPole(BrepCircle curve, {required bool coedgeReversed}) {
    final small = closestPoint(curve.center);
    final large = center * 2.0 - small;
    final n = coedgeReversed ? -curve.normal : curve.normal;
    return n.dot(large - curve.center) >= n.dot(small - curve.center)
        ? large
        : small;
  }

  /// Which spherical cap a loop fills, from the first coedge's sense.
  Vector3 loopFillPole(
    Vector3 p0,
    Vector3 tangent, {
    required bool coedgeReversed,
    required bool faceReversed,
    required Vector3 small,
    required Vector3 large,
  }) {
    var tang = Vector3.copy(tangent);
    if (coedgeReversed) tang = -tang;
    var outN = p0 - center;
    if (outN.length < 1e-12) outN = Vector3(0, 0, 1);
    outN = Vector3.copy(outN).normalized();
    if (faceReversed) outN = -outN;
    var inward = outN.cross(tang);
    if (inward.length < 1e-12) return small;
    final interior = closestPoint(p0 + inward.normalized() * radius * 0.25);
    return (interior - small).length <= (interior - large).length
        ? small
        : large;
  }

  BrepSphere offset(double d) => BrepSphere(center: center, radius: radius + d);

  @override
  Map<String, dynamic> toJson() => {
    't': 'sphere',
    'c': _v(center),
    'r': radius,
  };

  factory BrepSphere.fromJson(Map<String, dynamic> json) =>
      BrepSphere(center: _p(json['c']), radius: (json['r'] as num).toDouble());
}

class BrepTorus implements BrepSurface {
  BrepTorus({
    required this.center,
    required Vector3 axis,
    required this.majorRadius,
    required this.minorRadius,
    Vector3? xDir,
  }) : axis = axis.normalized(),
       xDir = BrepCylinder._planeX(axis.normalized(), xDir) {
    yDir = this.axis.cross(this.xDir).normalized();
  }

  final Vector3 center;
  final Vector3 axis;
  final double majorRadius;
  final double minorRadius;
  final Vector3 xDir;
  late final Vector3 yDir;

  @override
  Vector3 evaluate(double u, double v) {
    final cr = xDir * math.cos(u) + yDir * math.sin(u);
    return center +
        cr * (majorRadius + minorRadius * math.cos(v)) +
        axis * (minorRadius * math.sin(v));
  }

  @override
  Vector3 normalAt(double u, double v) {
    final cr = xDir * math.cos(u) + yDir * math.sin(u);
    return (cr * math.cos(v) + axis * math.sin(v)).normalized();
  }

  @override
  (double, double)? parameters(Vector3 p) {
    final d = p - center;
    final z = d.dot(axis);
    final rad = d - axis * z;
    final u = math.atan2(rad.dot(yDir), rad.dot(xDir));
    final rho = rad.length;
    final v = math.atan2(z, rho - majorRadius);
    return (u, v);
  }

  @override
  Vector3 closestPoint(Vector3 p) {
    final uv = parameters(p);
    if (uv == null) return evaluate(0, 0);
    return evaluate(uv.$1, uv.$2);
  }

  BrepTorus offset(double d) => BrepTorus(
    center: center,
    axis: axis,
    majorRadius: majorRadius,
    minorRadius: minorRadius + d,
    xDir: xDir,
  );

  @override
  Map<String, dynamic> toJson() => {
    't': 'torus',
    'c': _v(center),
    'a': _v(axis),
    'R': majorRadius,
    'r': minorRadius,
    'x': _v(xDir),
  };

  factory BrepTorus.fromJson(Map<String, dynamic> json) => BrepTorus(
    center: _p(json['c']),
    axis: _p(json['a']),
    majorRadius: (json['R'] as num).toDouble(),
    minorRadius: (json['r'] as num).toDouble(),
    xDir: _p(json['x']),
  );
}

class BrepCone implements BrepSurface {
  BrepCone({
    required this.origin,
    required Vector3 axis,
    required this.radius,
    required this.height,
    this.topRadius = 0,
    Vector3? xDir,
  }) : axis = axis.normalized(),
       xDir = BrepCylinder._planeX(axis.normalized(), xDir) {
    yDir = this.axis.cross(this.xDir).normalized();
  }

  final Vector3 origin;
  final Vector3 axis;
  final double radius;
  final double topRadius;
  final double height;
  final Vector3 xDir;
  late final Vector3 yDir;

  double radiusAt(double z) {
    if (height.abs() < 1e-18) return radius;
    return radius + (topRadius - radius) * (z / height);
  }

  @override
  Vector3 evaluate(double u, double v) {
    final r = radiusAt(v);
    return origin +
        axis * v +
        xDir * (r * math.cos(u)) +
        yDir * (r * math.sin(u));
  }

  @override
  Vector3 normalAt(double u, double v) {
    final sl = math.atan((radius - topRadius) / height);
    final cr = xDir * math.cos(u) + yDir * math.sin(u);
    return (cr * math.cos(sl) + axis * math.sin(sl)).normalized();
  }

  @override
  (double, double)? parameters(Vector3 p) {
    final d = p - origin;
    final v = d.dot(axis);
    final rad = d - axis * v;
    if (rad.length < 1e-12) return (0.0, v);
    return (math.atan2(rad.dot(yDir), rad.dot(xDir)), v);
  }

  @override
  Vector3 closestPoint(Vector3 p) {
    final d = p - origin;
    final z = d.dot(axis);
    final radial = d - axis * z;
    final rho = radial.length;
    final k = (topRadius - radius) / height;
    var v = (z + k * (rho - radius)) / (1 + k * k);
    if (k > 1e-15) v = math.max(v, -radius / k);
    if (k < -1e-15) v = math.min(v, -radius / k);
    final direction = rho > 1e-12 ? radial / rho : xDir;
    return origin + axis * v + direction * radiusAt(v);
  }

  @override
  Map<String, dynamic> toJson() => {
    't': 'cone',
    'o': _v(origin),
    'a': _v(axis),
    'r': radius,
    'rt': topRadius,
    'h': height,
    'x': _v(xDir),
  };

  factory BrepCone.fromJson(Map<String, dynamic> json) => BrepCone(
    origin: _p(json['o']),
    axis: _p(json['a']),
    radius: (json['r'] as num).toDouble(),
    topRadius: (json['rt'] as num?)?.toDouble() ?? 0,
    height: (json['h'] as num).toDouble(),
    xDir: _p(json['x']),
  );
}

sealed class BrepCurve {
  Vector3 evaluate(double t);
  Vector3 closestPoint(Vector3 p);
  double closestParam(Vector3 p);
  Vector3 tangent(double t);
  double get t0;
  double get t1;
  Map<String, dynamic> toJson();

  static BrepCurve fromJson(Map<String, dynamic> json) {
    switch (json['t'] as String) {
      case 'helix':
        return BrepHelix.fromJson(json);
      case 'spline':
        return BrepSplineCurve.fromJson(json);
      case 'surface_intersection':
        return BrepSurfaceIntersectionCurve.fromJson(json);
      case 'ruled_quadric':
        return BrepRuledQuadricCurve.fromJson(json);
      case 'plane_section':
        return BrepPlaneSectionCurve.fromJson(json);
      case 'line':
        return BrepLine.fromJson(json);
      case 'circle':
        return BrepCircle.fromJson(json);
      case 'ellipse':
        return BrepEllipse.fromJson(json);
      case 'intersection':
        return BrepIntersectionCurve.fromJson(json);
      case 'affine':
        return BrepAffineCurve.fromJson(json);
      default:
        throw FormatException('Unknown curve ${json['t']}');
    }
  }
}

class BrepLine implements BrepCurve {
  BrepLine(this.a, this.b);

  final Vector3 a;
  final Vector3 b;

  @override
  double get t0 => 0;
  @override
  double get t1 => 1;

  @override
  Vector3 evaluate(double t) => a + (b - a) * t;

  @override
  Vector3 tangent(double t) {
    final d = b - a;
    if (d.length < 1e-12) return Vector3(0, 0, 1);
    return d.normalized();
  }

  @override
  double closestParam(Vector3 p) {
    final d = b - a;
    final len2 = d.length2;
    if (len2 < 1e-20) return 0;
    return ((p - a).dot(d) / len2).clamp(0.0, 1.0);
  }

  @override
  Vector3 closestPoint(Vector3 p) => evaluate(closestParam(p));

  @override
  Map<String, dynamic> toJson() => {'t': 'line', 'a': _v(a), 'b': _v(b)};

  factory BrepLine.fromJson(Map<String, dynamic> json) =>
      BrepLine(_p(json['a']), _p(json['b']));
}

class BrepCircle implements BrepCurve {
  BrepCircle({
    required this.center,
    required Vector3 normal,
    required this.radius,
    Vector3? xDir,
    this.t0 = 0,
    this.t1 = math.pi * 2,
  }) : normal = normal.normalized(),
       xDir = BrepPlane._planeX(normal.normalized(), xDir) {
    yDir = this.normal.cross(this.xDir).normalized();
  }

  final Vector3 center;
  final Vector3 normal;
  final Vector3 xDir;
  late final Vector3 yDir;
  final double radius;
  @override
  final double t0;
  @override
  final double t1;

  bool get isClosed => (t1 - t0).abs() >= math.pi * 2 - 1e-6;

  @override
  Vector3 evaluate(double t) =>
      center + xDir * (radius * math.cos(t)) + yDir * (radius * math.sin(t));

  @override
  Vector3 tangent(double t) =>
      (-xDir * math.sin(t) + yDir * math.cos(t)).normalized();

  @override
  double closestParam(Vector3 p) {
    final d = p - center;
    final angle = math.atan2(d.dot(yDir), d.dot(xDir));
    final lo = math.min(t0, t1), hi = math.max(t0, t1);
    final middle = (lo + hi) / 2;
    final t = angle + ((middle - angle) / (2 * math.pi)).round() * 2 * math.pi;
    if (t >= lo && t <= hi) return t;
    // Compare actual endpoint distances; a rounding error near the start must
    // not wrap one full turn and select the opposite end of a short arc.
    return (evaluate(lo) - p).length2 <= (evaluate(hi) - p).length2 ? lo : hi;
  }

  @override
  Vector3 closestPoint(Vector3 p) => evaluate(closestParam(p));

  @override
  Map<String, dynamic> toJson() => {
    't': 'circle',
    'c': _v(center),
    'n': _v(normal),
    'x': _v(xDir),
    'r': radius,
    'a': t0,
    'b': t1,
  };

  factory BrepCircle.fromJson(Map<String, dynamic> json) => BrepCircle(
    center: _p(json['c']),
    normal: _p(json['n']),
    xDir: _p(json['x']),
    radius: (json['r'] as num).toDouble(),
    t0: (json['a'] as num?)?.toDouble() ?? 0,
    t1: (json['b'] as num?)?.toDouble() ?? math.pi * 2,
  );
}

class BrepEllipse implements BrepCurve {
  BrepEllipse({
    required this.center,
    required Vector3 normal,
    required this.xDir,
    required this.rMajor,
    required this.rMinor,
    this.t0 = 0,
    this.t1 = math.pi * 2,
  }) : normal = normal.normalized() {
    yDir = normal.cross(xDir.normalized()).normalized();
  }

  final Vector3 center;
  final Vector3 normal;
  final Vector3 xDir;
  late final Vector3 yDir;
  final double rMajor;
  final double rMinor;
  @override
  final double t0;
  @override
  final double t1;

  bool get isClosed => (t1 - t0).abs() >= math.pi * 2 - 1e-6;

  @override
  Vector3 evaluate(double t) =>
      center + xDir * (rMajor * math.cos(t)) + yDir * (rMinor * math.sin(t));

  @override
  Vector3 tangent(double t) {
    final d = -xDir * (rMajor * math.sin(t)) + yDir * (rMinor * math.cos(t));
    if (d.length < 1e-12) return xDir.clone();
    return d.normalized();
  }

  @override
  double closestParam(Vector3 p) {
    final d = p - center;
    final angle = math.atan2(
      d.dot(yDir) / rMinor,
      d.dot(xDir.normalized()) / (rMajor * xDir.length),
    );
    final lo = math.min(t0, t1), hi = math.max(t0, t1);
    final middle = (lo + hi) / 2;
    final t = angle + ((middle - angle) / (2 * math.pi)).round() * 2 * math.pi;
    if (t >= lo && t <= hi) return t;
    return (evaluate(lo) - p).length2 <= (evaluate(hi) - p).length2 ? lo : hi;
  }

  @override
  Vector3 closestPoint(Vector3 p) => evaluate(closestParam(p));

  @override
  Map<String, dynamic> toJson() => {
    't': 'ellipse',
    'c': _v(center),
    'n': _v(normal),
    'x': _v(xDir),
    'A': rMajor,
    'B': rMinor,
    'a': t0,
    'b': t1,
  };

  factory BrepEllipse.fromJson(Map<String, dynamic> json) => BrepEllipse(
    center: _p(json['c']),
    normal: _p(json['n']),
    xDir: _p(json['x']),
    rMajor: (json['A'] as num).toDouble(),
    rMinor: (json['B'] as num).toDouble(),
    t0: (json['a'] as num?)?.toDouble() ?? 0,
    t1: (json['b'] as num?)?.toDouble() ?? math.pi * 2,
  );
}

List<double> _v(Vector3 p) => [p.x, p.y, p.z];

Vector3 _p(dynamic json) {
  final l = json as List<dynamic>;
  return Vector3(
    (l[0] as num).toDouble(),
    (l[1] as num).toDouble(),
    (l[2] as num).toDouble(),
  );
}
