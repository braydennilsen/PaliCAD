import 'dart:math' as math;

import 'package:vector_math/vector_math_64.dart';

import '../brep/brep.dart';
import '../mesh/mesh_data.dart';
import 'face_domain.dart';

/// An axis-aligned box.
class IjoBounds {
  IjoBounds(this.min, this.max);

  final Vector3 min;
  final Vector3 max;

  static IjoBounds? fromPoints(Iterable<Vector3> points) {
    final acc = BoundsBuilder();
    for (final p in points) {
      acc.add(p);
    }
    return acc.result;
  }

  Vector3 get center => (min + max) * 0.5;
  Vector3 get size => max - min;
  double get diagonal => size.length;

  double get boxVolume {
    final s = size;
    return s.x * s.y * s.z;
  }

  bool containsPoint(Vector3 p, {double tolerance = 0}) =>
      p.x >= min.x - tolerance &&
      p.x <= max.x + tolerance &&
      p.y >= min.y - tolerance &&
      p.y <= max.y + tolerance &&
      p.z >= min.z - tolerance &&
      p.z <= max.z + tolerance;

  bool intersects(IjoBounds other, {double tolerance = 0}) =>
      min.x <= other.max.x + tolerance &&
      max.x >= other.min.x - tolerance &&
      min.y <= other.max.y + tolerance &&
      max.y >= other.min.y - tolerance &&
      min.z <= other.max.z + tolerance &&
      max.z >= other.min.z - tolerance;

  IjoBounds merge(IjoBounds other) => IjoBounds(
    Vector3(
      math.min(min.x, other.min.x),
      math.min(min.y, other.min.y),
      math.min(min.z, other.min.z),
    ),
    Vector3(
      math.max(max.x, other.max.x),
      math.max(max.y, other.max.y),
      math.max(max.z, other.max.z),
    ),
  );

  IjoBounds expanded(double d) =>
      IjoBounds(min - Vector3.all(d), max + Vector3.all(d));

  Map<String, dynamic> toJson() => {
    'min': [min.x, min.y, min.z],
    'max': [max.x, max.y, max.z],
  };

  @override
  String toString() => 'IjoBounds(${_f(min)} .. ${_f(max)})';

  static String _f(Vector3 v) =>
      '${v.x.toStringAsFixed(4)}, ${v.y.toStringAsFixed(4)}, '
      '${v.z.toStringAsFixed(4)}';
}

class BoundsBuilder {
  double _minX = double.infinity, _minY = double.infinity;
  double _minZ = double.infinity, _maxX = double.negativeInfinity;
  double _maxY = double.negativeInfinity, _maxZ = double.negativeInfinity;
  var _any = false;

  void add(Vector3 p) {
    if (!p.x.isFinite || !p.y.isFinite || !p.z.isFinite) return;
    _any = true;
    if (p.x < _minX) _minX = p.x;
    if (p.y < _minY) _minY = p.y;
    if (p.z < _minZ) _minZ = p.z;
    if (p.x > _maxX) _maxX = p.x;
    if (p.y > _maxY) _maxY = p.y;
    if (p.z > _maxZ) _maxZ = p.z;
  }

  IjoBounds? get result => _any
      ? IjoBounds(Vector3(_minX, _minY, _minZ), Vector3(_maxX, _maxY, _maxZ))
      : null;
}

IjoBounds? meshBounds(MeshData mesh) => IjoBounds.fromPoints(mesh.vertices);

/// Conservative sphere-face envelope. Circular cap constraints shrink the
/// supporting sphere without relying on sparse boundary samples.
IjoBounds conservativeSphereFaceBounds(BrepSolid solid, BrepFace face) {
  final sphere = face.surface as BrepSphere;
  final r = sphere.radius;
  final lo = sphere.center - Vector3.all(r);
  final hi = sphere.center + Vector3.all(r);
  final constraints = <(Vector3, double)>[];
  final closedCaps =
      face.loops.isNotEmpty &&
      face.loops.every(
        (loop) =>
            loop.coedges.length == 1 &&
            solid.edges[loop.coedges.single.edge].curve is BrepCircle &&
            (solid.edges[loop.coedges.single.edge].curve as BrepCircle)
                .isClosed,
      );
  if (closedCaps) {
    for (final loop in face.loops) {
      final co = loop.coedges.single;
      final circle = solid.edges[co.edge].curve as BrepCircle;
      var n = co.reversed ? -circle.normal : circle.normal;
      if (face.loops.length == 1) n = -n;
      constraints.add((n, (circle.center - sphere.center).dot(n)));
    }
  } else if (face.loops.length == 1) {
    // Convex great-circle wires (including spherical lunes with antipodal
    // corners) need no finite UV chart. Each oriented arc bounds an exact
    // hemisphere. Verify every complete arc, not just sampled vertices, before
    // using these constraints; large complements/concave trims stay unbounded.
    final uses = face.loops.single.coedges;
    final arcs = <BrepCircle>[];
    for (final co in uses) {
      final c = solid.edges[co.edge].curve;
      if (c is! BrepCircle ||
          c.t1 <= c.t0 ||
          c.t1 - c.t0 > math.pi + 1e-10 ||
          c.center.distanceTo(sphere.center) > r * 1e-10 ||
          (c.radius - r).abs() > r * 1e-10) {
        arcs.clear();
        break;
      }
      arcs.add(c);
    }
    if (arcs.length == uses.length && arcs.isNotEmpty) {
      final normals = [
        for (var i = 0; i < arcs.length; i++)
          arcs[i].normal * ((uses[i].reversed != face.reversed) ? -1.0 : 1.0),
      ];
      if (normals.every(
        (n) => arcs.every(
          (c) => _circleSupport(c, sphere.center, -n) <= r * 1e-10,
        ),
      )) {
        constraints.addAll(normals.map((n) => (-n, 0.0)));
      }
    }
  }
  for (final constraint in constraints) {
    final (n, h) = constraint;
    if (h.abs() > r) continue;
    double support(Vector3 axis) {
      final d = axis.dot(n);
      if (r * d <= h) return r;
      return h * d +
          math.sqrt(math.max(0, r * r - h * h)) *
              math.sqrt(math.max(0, 1 - d * d));
    }

    for (var i = 0; i < 3; i++) {
      final axis = Vector3.zero()..[i] = 1;
      lo[i] = math.max(lo[i], sphere.center[i] - support(-axis));
      hi[i] = math.min(hi[i], sphere.center[i] + support(axis));
    }
  }
  return IjoBounds(lo, hi);
}

// Exact support of an increasing circular arc along a unit direction.
double _circleSupport(BrepCircle c, Vector3 origin, Vector3 n) {
  var result = math.max(
    (c.evaluate(c.t0) - origin).dot(n),
    (c.evaluate(c.t1) - origin).dot(n),
  );
  final phase = math.atan2(c.yDir.dot(n), c.xDir.dot(n));
  final t = phase + ((c.t0 - phase) / (2 * math.pi)).ceil() * 2 * math.pi;
  if (t <= c.t1) result = math.max(result, (c.evaluate(t) - origin).dot(n));
  return result;
}

IjoBounds? _rectangularCylinderBounds(
  BrepSolid s,
  BrepFace f,
  BrepCylinder cylinder,
) {
  if (f.loops.length != 1 || f.loops.single.coedges.length != 4) return null;
  final curves = f.loops.single.coedges
      .map((c) => s.edges[c.edge].curve)
      .toList();
  final arcs = curves.whereType<BrepCircle>().toList();
  final lines = curves.whereType<BrepLine>().toList();
  if (arcs.length != 2 || lines.length != 2) return null;
  final tol = math.max(1e-9, cylinder.radius * 1e-10);
  if (arcs.any(
        (c) =>
            c.t1 <= c.t0 ||
            c.normal.cross(cylinder.axis).length > 1e-10 ||
            (c.radius - cylinder.radius).abs() > tol ||
            (c.center - cylinder.origin).cross(cylinder.axis).length > tol,
      ) ||
      lines.any((c) => (c.b - c.a).cross(cylinder.axis).length > tol)) {
    return null;
  }
  // The two translated rims must cover the same angular interval.
  Vector3 radial(Vector3 p) {
    final d = p - cylinder.origin;
    return d - cylinder.axis * d.dot(cylinder.axis);
  }

  final a = arcs[0], b = arcs[1];
  bool match(bool reverse) => [0.0, .5, 1.0].every(
    (t) =>
        radial(a.evaluate(a.t0 + (a.t1 - a.t0) * t)).distanceTo(
          radial(b.evaluate(b.t0 + (b.t1 - b.t0) * (reverse ? 1 - t : t))),
        ) <=
        tol,
  );
  if ((a.t1 - a.t0 - (b.t1 - b.t0)).abs() > 1e-10 ||
      (!match(false) && !match(true))) {
    return null;
  }
  final acc = BoundsBuilder();
  for (final c in curves) {
    _addCurveExtrema(acc, c);
  }
  return acc.result;
}

/// Conservative envelope for rejecting disjoint boolean operands. Unlike a
/// display/measurement bound, this may overestimate a trimmed surface, but
/// must never discard a possible intersection. Returns null if a cylindrical
/// face has no boundary from which to bound its axial extent.
IjoBounds? conservativeSolidBounds(BrepSolid solid) {
  final transformed = BrepXform.unwrapAffine(solid);
  if (transformed != null) {
    final box = conservativeSolidBounds(transformed.solid);
    return box == null ? null : _affineBounds(box, transformed.transform);
  }
  final acc = BoundsBuilder();
  void ball(Vector3 center, double radius, BoundsBuilder into) {
    into.add(center - Vector3.all(radius));
    into.add(center + Vector3.all(radius));
  }

  for (final face in solid.faces) {
    final surface = face.surface;
    if (surface is BrepRuledSurface &&
        surface.lower is BrepHelix &&
        surface.upper is BrepHelix) {
      // Every ruling is a convex combination of its two exact rail points.
      for (final p in (surface.lower as BrepHelix).boundsPoints) acc.add(p);
      for (final p in (surface.upper as BrepHelix).boundsPoints) acc.add(p);
      continue;
    }
    if (surface is BrepSplineSurface) {
      for (final p in surface.points.expand((r) => r)) {
        acc.add(p);
      }
      continue;
    }
    if (surface is BrepSphere) {
      final box = conservativeSphereFaceBounds(solid, face);
      acc.add(box.min);
      acc.add(box.max);
      continue;
    }
    if (surface is BrepCylinder) {
      final box = _rectangularCylinderBounds(solid, face, surface);
      if (box != null) {
        acc.add(box.min);
        acc.add(box.max);
        continue;
      }
    }
    if (surface is BrepTorus) {
      ball(surface.center, surface.majorRadius + surface.minorRadius, acc);
      continue;
    }
    final boundary = BoundsBuilder();
    for (final loop in face.loops) {
      for (final c in loop.coedges) {
        final edge = solid.edges[c.edge];
        boundary.add(solid.vertices[edge.v0]);
        boundary.add(solid.vertices[edge.v1]);
        switch (edge.curve) {
          case BrepHelix curve:
            _addCurveExtrema(boundary, curve);
          case BrepSurfaceIntersectionCurve curve:
            _addCurveExtrema(boundary, curve);
          case BrepRuledQuadricCurve curve:
            _addCurveExtrema(boundary, curve);
          case BrepPlaneSectionCurve curve:
            _addCurveExtrema(boundary, curve);
          case BrepSplineCurve curve:
            for (final p in curve.points) {
              boundary.add(p);
            }
          case BrepAffineCurve curve:
            _addCurveExtrema(boundary, curve);
          case BrepIntersectionCurve curve:
            final (lo, hi) = curve.axialBounds;
            ball(
              curve.carrier.origin + curve.carrier.axis * lo,
              curve.carrier.radius,
              boundary,
            );
            ball(
              curve.carrier.origin + curve.carrier.axis * hi,
              curve.carrier.radius,
              boundary,
            );
          case BrepLine curve:
            boundary.add(curve.a);
            boundary.add(curve.b);
          case BrepCircle curve:
            if (surface is BrepPlane && curve.t1 >= curve.t0) {
              _addCurveExtrema(boundary, curve);
              break;
            }
            final extent = Vector3.zero();
            for (var i = 0; i < 3; i++) {
              extent[i] =
                  curve.radius *
                  math.sqrt(
                    curve.xDir[i] * curve.xDir[i] +
                        curve.yDir[i] * curve.yDir[i],
                  );
            }
            boundary.add(curve.center - extent);
            boundary.add(curve.center + extent);
          case BrepEllipse curve:
            final extent = Vector3.zero();
            for (var i = 0; i < 3; i++) {
              extent[i] = math.sqrt(
                math.pow(curve.xDir[i] * curve.rMajor, 2) +
                    math.pow(curve.yDir[i] * curve.rMinor, 2),
              );
            }
            boundary.add(curve.center - extent);
            boundary.add(curve.center + extent);
        }
      }
    }
    final box = boundary.result;
    if (box == null) return null;
    acc.add(box.min);
    acc.add(box.max);
    if (surface is BrepCylinder || surface is BrepCone) {
      final axis = surface is BrepCylinder
          ? surface.axis
          : (surface as BrepCone).axis;
      final origin = surface is BrepCylinder
          ? surface.origin
          : (surface as BrepCone).origin;
      var z0 = double.infinity, z1 = double.negativeInfinity;
      for (final x in [box.min.x, box.max.x]) {
        for (final y in [box.min.y, box.max.y]) {
          for (final z in [box.min.z, box.max.z]) {
            final along = (Vector3(x, y, z) - origin).dot(axis);
            z0 = math.min(z0, along);
            z1 = math.max(z1, along);
          }
        }
      }
      final radius = surface is BrepCylinder
          ? surface.radius
          : math.max(
              (surface as BrepCone).radiusAt(z0).abs(),
              surface.radiusAt(z1).abs(),
            );
      final extent = Vector3.zero();
      for (var i = 0; i < 3; i++) {
        extent[i] = radius * math.sqrt(math.max(0, 1 - axis[i] * axis[i]));
      }
      for (final z in [z0, z1]) {
        acc.add(origin + axis * z - extent);
        acc.add(origin + axis * z + extent);
      }
    }
  }
  return acc.result;
}

/// Axis-aligned bounds from curves and surface extrema. Affine images of
/// general solids and implicit intersection curves may have conservative
/// bounds; full spheres and ellipsoids use exact support extents.
IjoBounds? solidBounds(BrepSolid solid) {
  final transformed = BrepXform.unwrapAffine(solid);
  if (transformed != null) {
    final base = transformed.solid, m = transformed.transform;
    if (base.faces.length == 1 &&
        base.faces.single.surface is BrepSphere &&
        base.faces.single.loops.isEmpty) {
      final sphere = base.faces.single.surface as BrepSphere;
      final c = m.point(sphere.center);
      final extent =
          Vector3(
            m.transposeVector(Vector3(1, 0, 0)).length,
            m.transposeVector(Vector3(0, 1, 0)).length,
            m.transposeVector(Vector3(0, 0, 1)).length,
          ) *
          sphere.radius;
      return IjoBounds(c - extent, c + extent);
    }
    final box = solidBounds(base);
    return box == null ? null : _affineBounds(box, m);
  }
  final acc = BoundsBuilder();
  for (final v in solid.vertices) {
    acc.add(v);
  }
  for (final e in solid.edges) {
    _addCurveExtrema(acc, e.curve);
  }
  for (final face in solid.faces) {
    _addSurfaceBulge(acc, solid, face);
  }
  return acc.result;
}

void _addCurveExtrema(BoundsBuilder acc, BrepCurve curve) {
  switch (curve) {
    case BrepHelix():
      for (final p in curve.boundsPoints) acc.add(p);
    case BrepSurfaceIntersectionCurve():
      for (final row in curve.first.points) for (final p in row) acc.add(p);
    case BrepRuledQuadricCurve():
      _addCurveExtrema(acc, curve.surface.lower);
      _addCurveExtrema(acc, curve.surface.upper);
    case BrepPlaneSectionCurve():
      _addCurveExtrema(acc, curve.surface.lower);
      _addCurveExtrema(acc, curve.surface.upper);
    case BrepSplineCurve():
      for (final p in curve.points) {
        acc.add(p);
      }
    case BrepAffineCurve():
      final bounds = BoundsBuilder();
      _addCurveExtrema(bounds, curve.base);
      final box = bounds.result;
      if (box != null) {
        final out = _affineBounds(box, curve.transform);
        acc.add(out.min);
        acc.add(out.max);
      }
    case BrepIntersectionCurve():
      final (lo, hi) = curve.axialBounds;
      for (final z in [lo, hi]) {
        final p = curve.carrier.origin + curve.carrier.axis * z;
        acc.add(p - Vector3.all(curve.carrier.radius));
        acc.add(p + Vector3.all(curve.carrier.radius));
      }
    case BrepLine():
      acc.add(curve.a);
      acc.add(curve.b);
    case BrepCircle():
      _addArcExtrema(
        acc,
        curve.xDir * curve.radius,
        curve.yDir * curve.radius,
        curve.t0,
        curve.isClosed ? curve.t0 + math.pi * 2 : curve.t1,
        curve.evaluate,
      );
    case BrepEllipse():
      _addArcExtrema(
        acc,
        curve.xDir.normalized() * curve.rMajor,
        curve.yDir * curve.rMinor,
        curve.t0,
        curve.t1,
        curve.evaluate,
      );
  }
}

void _addArcExtrema(
  BoundsBuilder acc,
  Vector3 a,
  Vector3 b,
  double t0,
  double t1,
  Vector3 Function(double) eval,
) {
  acc.add(eval(t0));
  acc.add(eval(t1));

  const twoPi = math.pi * 2;
  final full = (t1 - t0).abs() >= twoPi - 1e-9;
  for (var i = 0; i < 3; i++) {
    final phi = math.atan2(b[i], a[i]);
    for (final candidate in [phi, phi + math.pi]) {
      if (full) {
        acc.add(eval(candidate));
        continue;
      }
      var t = candidate;
      while (t < t0) {
        t += twoPi;
      }
      while (t > t0 + twoPi) {
        t -= twoPi;
      }
      if (t <= t1 + 1e-12) acc.add(eval(t));
    }
  }
}

void _addSurfaceBulge(BoundsBuilder acc, BrepSolid solid, BrepFace face) {
  final surface = face.surface;
  if (surface is BrepSplineSurface) {
    for (final p in surface.points.expand((r) => r)) {
      acc.add(p);
    }
    return;
  }
  if (surface is! BrepSphere && surface is! BrepTorus) return;

  final FaceDomain domain;
  try {
    domain = faceDomain(solid, face);
  } on BrepUnsupported {
    return;
  }

  bool inside(double u, double v) => switch (domain) {
    FaceDomainRect() => domain.contains(u, v, uPeriodic: true),
    FaceDomainLoops() => domain.contains(u, v),
  };

  for (final (u, v) in _bulgeCandidates(surface)) {
    if (inside(u, v)) acc.add(surface.evaluate(u, v));
  }
}

List<(double, double)> _bulgeCandidates(BrepSurface surface) {
  const half = math.pi / 2;
  switch (surface) {
    case BrepSphere():
      return const [
        (0, half),
        (math.pi, half),
        (half, half),
        (-half, half),
        (0, 0),
        (0, math.pi),
      ];
    case BrepTorus():
      final out = <(double, double)>[];
      for (var i = 0; i < 3; i++) {
        final alpha = surface.xDir[i];
        final beta = surface.yDir[i];
        final gamma = surface.axis[i];
        final phi = math.atan2(beta, alpha);
        final rho = math.sqrt(alpha * alpha + beta * beta);
        for (final (u, sign) in [(phi, 1.0), (phi + math.pi, -1.0)]) {
          final v = math.atan2(gamma, sign * rho);
          out.add((u, v));
          out.add((u, v + math.pi));
        }
      }
      return out;
    default:
      return const [];
  }
}

IjoBounds _affineBounds(IjoBounds box, BrepAffineTransform transform) =>
    IjoBounds.fromPoints([
      for (final x in [box.min.x, box.max.x])
        for (final y in [box.min.y, box.max.y])
          for (final z in [box.min.z, box.max.z])
            transform.point(Vector3(x, y, z)),
    ])!;
