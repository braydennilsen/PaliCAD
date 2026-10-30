import 'dart:math' as math;

import 'package:vector_math/vector_math_64.dart';

import 'geom.dart';
import 'solid.dart';
import 'tessellate.dart';
import 'xform.dart';
import 'intersect.dart';
import 'spline_intersect.dart';
import '../measure/bounds.dart';

enum BrepSide { outside, on, inside }

/// Inside / outside tests against an analytic solid.
class BrepClassify {
  BrepClassify._();

  static Expando<Map<BrepFace, List<List<Vector3>>>>? _trimCache;
  static Expando<IjoBounds>? _boundsCache;
  static Expando<IjoBounds>? _helicalBounds;
  static IjoBounds? _helixBounds(BrepSurface surface) {
    if (surface is! BrepRuledSurface ||
        surface.lower is! BrepHelix ||
        surface.upper is! BrepHelix)
      return null;
    final cached = _helicalBounds?[surface];
    if (cached != null) return cached;
    final box = IjoBounds.fromPoints([
      ...(surface.lower as BrepHelix).boundsPoints,
      ...(surface.upper as BrepHelix).boundsPoints,
    ])!;
    if (_helicalBounds != null) _helicalBounds![surface] = box;
    return box;
  }

  static double _boxDistance(IjoBounds box, Vector3 p) {
    var d = 0.0;
    for (var i = 0; i < 3; i++) {
      final delta = p[i] - p[i].clamp(box.min[i], box.max[i]);
      d += delta * delta;
    }
    return math.sqrt(d);
  }

  static bool _rayBox(IjoBounds box, Vector3 p, Vector3 ray, double tol) {
    var near = 0.0, far = double.infinity;
    for (var i = 0; i < 3; i++) {
      if (ray[i].abs() < 1e-20) {
        if (p[i] < box.min[i] - tol || p[i] > box.max[i] + tol) return false;
        continue;
      }
      final a = (box.min[i] - tol - p[i]) / ray[i],
          b = (box.max[i] + tol - p[i]) / ray[i];
      near = math.max(near, math.min(a, b));
      far = math.min(far, math.max(a, b));
      if (far < near) return false;
    }
    return true;
  }

  /// Reuse parameterized trims only during a synchronous modelling operation.
  /// Restore the outer scope on reentry; do not retain data across model edits.
  static T withTrimCache<T>(T Function() action) {
    final previous = _trimCache,
        previousBounds = _boundsCache,
        previousHelical = _helicalBounds;
    _trimCache = Expando();
    _boundsCache = Expando();
    _helicalBounds = Expando();
    try {
      return action();
    } finally {
      _trimCache = previous;
      _boundsCache = previousBounds;
      _helicalBounds = previousHelical;
    }
  }

  static List<List<Vector3>> _trimPolygons(BrepSolid solid, BrepFace face) {
    List<List<Vector3>> build() {
      final period = _closedRuledUPeriod(face.surface);
      return [
        for (final loop in face.loops)
          () {
            final points = <Vector3>[];
            for (final p in BrepTessellate.sampleLoop(solid, loop, 64)) {
              final uv = face.surface.parameters(p)!;
              var u = uv.$1;
              // The two sides of a closed sweep's seam occupy the same world
              // points but different UVs. Keep the wire continuous in its chart;
              // otherwise it cuts diagonally across the wall and loses material.
              if (period != null && points.isNotEmpty) {
                u += ((points.last.x - u) / period).round() * period;
              }
              points.add(Vector3(u, uv.$2, 0));
            }
            return points;
          }(),
      ];
    }

    final cache = _trimCache;
    if (cache == null) return build();
    return (cache[solid] ??= {}).putIfAbsent(face, build);
  }

  static double? _closedRuledUPeriod(BrepSurface surface) {
    final ruled = surface is BrepRuledSurface
        ? surface
        : surface is BrepSplineSurface
        ? surface.ruled
        : null;
    if (ruled == null ||
        (ruled.evaluate(0, 0) - ruled.evaluate(1, 0)).length > 1e-8 ||
        (ruled.evaluate(0, 1) - ruled.evaluate(1, 1)).length > 1e-8) {
      return null;
    }
    return surface is BrepSplineSurface ? surface.u1 - surface.u0 : 1;
  }

  static bool isInside(BrepSolid solid, Vector3 p, {double onTol = 1e-6}) {
    final s = side(solid, p, onTol: onTol);
    return s == BrepSide.inside || s == BrepSide.on;
  }

  static BrepSide side(BrepSolid solid, Vector3 p, {double onTol = 1e-6}) {
    if (solid.faces.isEmpty) return BrepSide.outside;
    final transformed = BrepXform.unwrapAffine(solid);
    if (transformed != null) {
      final m = transformed.transform;
      final bound = math.sqrt(
        m.inverseVector(Vector3(1, 0, 0)).length2 +
            m.inverseVector(Vector3(0, 1, 0)).length2 +
            m.inverseVector(Vector3(0, 0, 1)).length2,
      );
      return side(transformed.solid, m.inversePoint(p), onTol: onTol * bound);
    }
    final analytic = _analyticPrimitiveSide(solid, p, onTol);
    if (analytic != null) return analytic;
    var bounds = _boundsCache?[solid];
    if (bounds == null) {
      bounds = conservativeSolidBounds(solid);
      if (bounds != null && _boundsCache != null) _boundsCache![solid] = bounds;
    }
    if (bounds != null && !bounds.containsPoint(p, tolerance: onTol))
      return BrepSide.outside;
    var best = double.infinity;
    Vector3? qBest;
    Vector3? nBest;
    for (final face in solid.faces) {
      final box = _helixBounds(face.surface);
      if (box != null && _boxDistance(box, p) > onTol) continue;
      final q = closestOnFace(solid, face, p);
      final d = (p - q).length;
      if (d < best) {
        best = d;
        qBest = q;
        final uv = face.surface.parameters(q) ?? (0.0, 0.0);
        nBest = face.outwardNormal(uv.$1, uv.$2);
      }
    }
    if (best < onTol) return BrepSide.on;
    final raySide = _analyticRaySide(solid, p, onTol);
    if (raySide != null) return raySide;
    // Only an unresolved ray needs a nearest-normal fallback. Distant screw
    // patches cannot make a point lie on the boundary, so defer those costly
    // closest-point searches until this fallback actually needs them.
    for (final face in solid.faces) {
      final box = _helixBounds(face.surface);
      if (box == null ||
          _boxDistance(box, p) <= onTol ||
          _boxDistance(box, p) > best)
        continue;
      final q = closestOnFace(solid, face, p), d = (p - q).length;
      if (d < best) {
        best = d;
        qBest = q;
        final uv = face.surface.parameters(q)!;
        nBest = face.outwardNormal(uv.$1, uv.$2);
      }
    }
    if (qBest == null || nBest == null) return BrepSide.outside;
    if (solid.faces.any((f) => f.surface is BrepSplineSurface))
      throw BrepUnsupported(
        'Spline material classification has unresolved ray intersections',
      );
    return (p - qBest).dot(nBest) < 0 ? BrepSide.inside : BrepSide.outside;
  }

  // A nearest face normal is not a membership test at concave joins or at
  // prism edges/vertices: several bounded faces can share the closest point.
  // Count crossings for planar/cylindrical solids too, including triangular
  // cutters. Circularly trimmed spheres are supported in the same path.
  // Unsupported trim families retain the existing path.
  static BrepSide? _analyticRaySide(BrepSolid solid, Vector3 p, double tol) {
    for (final f in solid.faces) {
      if (f.surface is BrepSphere) {
        if (!_circularSphereFace(solid, f) &&
            BrepTessellate.sphereFaceFrame(solid, f) == null)
          return null;
      } else if (f.surface is! BrepPlane &&
          f.surface is! BrepTorus &&
          f.surface is! BrepRuledSurface &&
          f.surface is! BrepSplineSurface &&
          !BrepQuadric.supports(f.surface)) {
        return null;
      }
    }
    for (final direction in [
      Vector3(1, .371, .529),
      Vector3(.293, 1, .617),
      Vector3(.431, .239, 1),
    ]) {
      final ray = direction.normalized();
      final hits = <double>[];
      var ambiguous = false;
      for (final face in solid.faces) {
        final surface = face.surface;
        final box = _helixBounds(surface);
        if (box != null && !_rayBox(box, p, ray, tol)) continue;
        final roots = <double>[];
        if (surface is BrepPlane) {
          final den = surface.normal.dot(ray);
          if (den.abs() < 1e-12) continue;
          roots.add(-surface.height(p) / den);
        } else if (surface is BrepTorus) {
          roots.addAll(BrepIntersect.rayTorus(surface, p, ray));
        } else if (surface is BrepSplineSurface) {
          final hits = BrepSplineIntersect.ray(surface, p, ray);
          if (hits == null) {
            ambiguous = true;
            break;
          }
          roots.addAll(hits);
        } else if (surface is BrepRuledSurface &&
            !BrepQuadric.supports(surface)) {
          final hits = BrepIntersect.rayRuled(surface, p, ray);
          if (hits == null) {
            ambiguous = true;
            break;
          }
          roots.addAll(hits);
        } else {
          final coefficients = BrepQuadric.line(surface, p, ray);
          final a = coefficients.$1,
              b = coefficients.$2 / 2,
              c = coefficients.$3;
          if (a.abs() < 1e-20) {
            if (b.abs() > 1e-20) roots.add(-c / (2 * b));
          } else {
            final disc = b * b - a * c;
            if (disc < 0) continue;
            final root = math.sqrt(disc);
            // Stable quadratic roots, including a ray starting near a surface.
            final q = -b - (b >= 0 ? root : -root);
            if (q.abs() < 1e-20) {
              roots.add(-b / a);
            } else {
              roots.addAll([q / a, c / q]);
            }
          }
        }
        for (final t in roots) {
          if (t <= tol) continue;
          final hit = p + ray * t;
          if (surface is BrepSphere) {
            if (!_inSphereFace(solid, face, hit, tol)) continue;
          } else if ((closestOnFace(solid, face, hit) - hit).length > tol) {
            continue;
          }
          final uv = surface.parameters(hit)!;
          if (surface.normalAt(uv.$1, uv.$2).dot(ray).abs() < 1e-8) {
            ambiguous = true;
            break;
          }
          // Do not count a ray through an edge twice (or guess at a vertex).
          for (final loop in face.loops) {
            for (final co in loop.coedges) {
              if ((solid.edges[co.edge].curve.closestPoint(hit) - hit).length <=
                  tol) {
                ambiguous = true;
              }
            }
          }
          hits.add(t);
        }
        if (ambiguous) break;
      }
      if (ambiguous) continue;
      return hits.length.isOdd ? BrepSide.inside : BrepSide.outside;
    }
    return null;
  }

  static bool _circularSphereFace(
    BrepSolid solid,
    BrepFace face,
  ) => face.loops.every(
    (loop) =>
        loop.coedges.length == 1 &&
        solid.edges[loop.coedges.single.edge].curve is BrepCircle &&
        (solid.edges[loop.coedges.single.edge].curve as BrepCircle).isClosed,
  );

  static bool _inSphereFace(
    BrepSolid solid,
    BrepFace face,
    Vector3 q,
    double tol,
  ) {
    final frame = BrepTessellate.sphereFaceFrame(solid, face);
    if (frame != null) {
      if ((q - frame.origin).dot(frame.normal) <= 0) return false;
      // Use exact boundary sides close to a trim; chordal polygons otherwise
      // misclassify the small probes used to decide which boolean side to keep.
      var closest = double.infinity;
      bool? nearSide;
      for (final loop in face.loops) {
        for (final co in loop.coedges) {
          final c = solid.edges[co.edge].curve,
              t = c.closestParam(q),
              p = c.evaluate(t);
          final distance = (p - q).length;
          if (distance < tol) return true;
          if (distance >= 1e-3 ||
              distance >= closest ||
              t <= c.t0 + 1e-8 ||
              t >= c.t1 - 1e-8)
            continue;
          var tangent = c.tangent(t);
          if (co.reversed != face.reversed) tangent = -tangent;
          final inward = (p - frame.origin).normalized().cross(tangent);
          nearSide = (q - p).dot(inward) >= 0;
          closest = distance;
        }
      }
      if (nearSide != null) return nearSide;
      final uv = BrepTessellate.sphereProject(frame, q);
      var inside = false;
      for (final loop in face.loops) {
        final poly = [
          for (final p in BrepTessellate.sampleLoop(solid, loop, 96))
            BrepTessellate.sphereProject(frame, p),
        ];
        if (_inPoly(
          [for (final p in poly) Vector3(p.x, p.y, 0)],
          Vector3(uv.x, uv.y, 0),
          Vector3(0, 0, 1),
        ))
          inside = !inside;
      }
      return inside;
    }
    for (final loop in face.loops) {
      final co = loop.coedges.single;
      final circle = solid.edges[co.edge].curve as BrepCircle;
      final n = co.reversed ? -circle.normal : circle.normal;
      final side = (q - circle.center).dot(n);
      // One loop encloses a cap; multiple closed loops exclude caps from the
      // remainder (including the band between two cube attachment circles).
      if (face.loops.length == 1 ? side < -tol : side > tol) return false;
    }
    return true;
  }

  /// Closed cone / cylinder / sphere primitives: closest-face + normal is
  /// wrong near a bounded cap rim (the rim wins, then the cap half-space
  /// swallows points that sit outside the wall).
  static BrepSide? _analyticPrimitiveSide(
    BrepSolid solid,
    Vector3 p,
    double onTol,
  ) {
    BrepCone? cone;
    BrepSphere? sph;
    var mixed = false;
    for (final face in solid.faces) {
      final s = face.surface;
      if (s is BrepCone) {
        if (cone != null) return null;
        cone = s;
      } else if (s is BrepSphere) {
        if (sph != null) return null;
        sph = s;
      } else if (s is BrepPlane) {
        // Caps of a cone primitive.
      } else {
        mixed = true;
      }
    }
    if (mixed || (cone != null && sph != null)) return null;
    if (cone != null && solid.faces.length <= 2) {
      return _coneSolidSide(cone, p, onTol);
    }
    if (sph != null &&
        solid.faces.length == 1 &&
        solid.faces.single.loops.isEmpty) {
      final d = (p - sph.center).length - sph.radius;
      if (d.abs() <= onTol) return BrepSide.on;
      return d < 0 ? BrepSide.inside : BrepSide.outside;
    }
    return null;
  }

  static BrepSide _coneSolidSide(BrepCone cone, Vector3 p, double onTol) {
    final d = p - cone.origin;
    final z = d.dot(cone.axis);
    final rad = (d - cone.axis * z).length;
    if (z < -onTol || z > cone.height + onTol) return BrepSide.outside;
    final r = cone.radiusAt(z.clamp(0.0, cone.height));
    if (z <= onTol) {
      if (rad > cone.radius + onTol) return BrepSide.outside;
      return BrepSide.on;
    }
    if ((cone.height - z).abs() <= onTol) {
      if (rad > cone.topRadius + onTol) return BrepSide.outside;
      return BrepSide.on;
    }
    if ((rad - r).abs() <= onTol) return BrepSide.on;
    return rad < r ? BrepSide.inside : BrepSide.outside;
  }

  /// Exact even-odd containment for planar line/conic wires. Split conics at
  /// ray-coordinate extrema so tangencies never become spurious crossings.
  /// Null asks the caller to use the general trim path for other curve types.
  static bool? containsPlaneFace(BrepSolid solid, BrepFace face, Vector3 p) {
    final plane = face.surface;
    if (plane is! BrepPlane || face.loops.isEmpty) return null;
    final curves = [
      for (final l in face.loops)
        for (final co in l.coedges) solid.edges[co.edge].curve,
    ];
    if (curves.any(
      (c) => c is! BrepLine && c is! BrepCircle && c is! BrepEllipse,
    ))
      return null;
    // Work relative to the query, and give shared trim endpoints the same
    // half-open ray convention despite tiny trigonometric roundoff at seams.
    (double, double) rayCoordinates(Vector3 point) {
      final d = point - p, y = d.dot(plane.yDir);
      return (d.dot(plane.xDir), y.abs() < 1e-10 ? 0.0 : y);
    }

    var inside = false;
    for (final curve in curves) {
      if ((curve.closestPoint(p) - p).length < 1e-8) return true;
      final lo = math.min(curve.t0, curve.t1),
          hi = math.max(curve.t0, curve.t1);
      final cuts = <double>[lo, hi];
      if (curve is BrepCircle || curve is BrepEllipse) {
        final x = curve is BrepCircle
            ? curve.xDir * curve.radius
            : (curve as BrepEllipse).xDir * curve.rMajor;
        final y = curve is BrepCircle
            ? curve.yDir * curve.radius
            : (curve as BrepEllipse).yDir * curve.rMinor;
        final extremum = math.atan2(y.dot(plane.yDir), x.dot(plane.yDir));
        for (
          var k = ((lo - extremum) / math.pi).ceil();
          extremum + k * math.pi < hi;
          k++
        ) {
          final t = extremum + k * math.pi;
          if (t > lo) cuts.add(t);
        }
      }
      cuts.sort();
      for (var i = 1; i < cuts.length; i++) {
        var a = cuts[i - 1], b = cuts[i];
        final pa = rayCoordinates(curve.evaluate(a)),
            pb = rayCoordinates(curve.evaluate(b));
        if ((pa.$2 > 0.0) == (pb.$2 > 0.0)) continue;
        double t;
        if (curve is BrepLine) {
          t = a + (b - a) * (0.0 - pa.$2) / (pb.$2 - pa.$2);
        } else {
          final increasing = pb.$2 > pa.$2;
          for (var j = 0; j < 45; j++) {
            final m = (a + b) / 2, y = (curve.evaluate(m) - p).dot(plane.yDir);
            if ((y < 0.0) == increasing) {
              a = m;
            } else {
              b = m;
            }
          }
          t = (a + b) / 2;
        }
        if ((curve.evaluate(t) - p).dot(plane.xDir) > 0) inside = !inside;
      }
    }
    return inside;
  }

  static Vector3 closestOnFace(BrepSolid solid, BrepFace face, Vector3 p) {
    final s = face.surface;
    final q = s.closestPoint(p);
    if (s is BrepSphere &&
        face.loops.isNotEmpty &&
        (_circularSphereFace(solid, face) ||
            BrepTessellate.sphereFaceFrame(solid, face) != null)) {
      if (_inSphereFace(solid, face, q, 1e-9)) return q;
      Vector3? best;
      var distance = double.infinity;
      for (final co in face.loops.expand((l) => l.coedges)) {
        final curve = solid.edges[co.edge].curve;
        final candidate = curve.closestPoint(p);
        if (!_inSphereFace(solid, face, candidate, 1e-8)) continue;
        final d = (candidate - p).length2;
        if (d < distance) {
          best = candidate;
          distance = d;
        }
      }
      if (best != null) return best;
    }
    if ((s is BrepRuledSurface || s is BrepSplineSurface) &&
        face.loops.isNotEmpty) {
      final uv = s.parameters(q)!;
      var inside = false;
      final period = _closedRuledUPeriod(s);
      for (final poly in _trimPolygons(solid, face)) {
        var u = uv.$1;
        if (period != null && poly.isNotEmpty) {
          final center = poly.fold(0.0, (double a, p) => a + p.x) / poly.length;
          u += ((center - u) / period).round() * period;
        }
        bool contains(double x) =>
            _inPoly(poly, Vector3(x, uv.$2, 0), Vector3(0, 0, 1));
        // Unequal sampling can bias a full-period polygon's mean away from
        // its midpoint. Check adjacent charts too, so neither seam side is lost.
        if (contains(u) ||
            period != null && (contains(u - period) || contains(u + period))) {
          inside = !inside;
        }
      }
      if (!inside) {
        Vector3? best;
        var distance = double.infinity;
        for (final loop in face.loops) {
          for (final co in loop.coedges) {
            final c = solid.edges[co.edge].curve.closestPoint(p),
                d = (c - p).length2;
            if (d < distance) {
              best = c;
              distance = d;
            }
          }
        }
        if (best != null) return best;
      }
    }
    if (s is BrepPlane && face.loops.isNotEmpty) {
      final inside = containsPlaneFace(solid, face, q);
      if (inside == true) return q;
      if (inside == false) {
        Vector3? best;
        var distance = double.infinity;
        for (final loop in face.loops)
          for (final co in loop.coedges) {
            final candidate = solid.edges[co.edge].curve.closestPoint(p);
            final d = (candidate - p).length2;
            if (d < distance) {
              best = candidate;
              distance = d;
            }
          }
        if (best != null) return best;
      }
      final poly = BrepTessellate.sampleLoop(solid, face.loops.first, 48);
      if (poly.length >= 3 && !_inPoly(poly, q, s.normal)) {
        return _closestOnPoly(poly, p);
      }
      for (var i = 1; i < face.loops.length; i++) {
        final hole = BrepTessellate.sampleLoop(solid, face.loops[i], 48);
        if (hole.length >= 3 && _inPoly(hole, q, s.normal)) {
          return _closestOnPoly(hole, p);
        }
      }
    }
    if (s is BrepCylinder && face.loops.isNotEmpty) {
      return _closestOnCylinderPatch(solid, face, p, q);
    }
    if (s is BrepTorus && face.loops.isNotEmpty) {
      final domain = BrepTessellate.torusFaceDomain(solid, s, face);
      if (domain != null) {
        final uv = s.parameters(q)!;
        double clampAngle(double angle, double lo, double hi) {
          angle +=
              (((lo + hi) / 2 - angle) / (2 * math.pi)).round() * 2 * math.pi;
          return angle.clamp(lo, hi);
        }

        return s.evaluate(
          clampAngle(uv.$1, domain.$1, domain.$2),
          clampAngle(uv.$2, domain.$3, domain.$4),
        );
      }
    }
    if (s is BrepCone) {
      final cone = s;
      final uv = cone.parameters(q) ?? cone.parameters(p);
      if (uv != null) {
        final domain = BrepTessellate.coneFaceDomain(solid, face);
        var u = uv.$1;
        u +=
            (((domain.$1 + domain.$2) / 2 - u) / (2 * math.pi)).round() *
            2 *
            math.pi;
        return cone.evaluate(
          u.clamp(domain.$1, domain.$2),
          uv.$2.clamp(domain.$3, domain.$4),
        );
      }
    }
    return q;
  }

  /// Infinite-cylinder [closestPoint] is wrong for a wall *patch* (D-notch
  /// half-cylinder). Snap onto the face's U/V domain.
  static Vector3 _closestOnCylinderPatch(
    BrepSolid solid,
    BrepFace face,
    Vector3 p,
    Vector3 infinite,
  ) {
    final cyl = face.surface as BrepCylinder;
    if (face.loops.any(
      (l) => l.coedges.any(
        (c) =>
            solid.edges[c.edge].curve is BrepSurfaceIntersectionCurve ||
            solid.edges[c.edge].curve is BrepRuledQuadricCurve ||
            solid.edges[c.edge].curve is BrepPlaneSectionCurve ||
            solid.edges[c.edge].curve is BrepSplineCurve ||
            solid.edges[c.edge].curve is BrepAffineCurve ||
            solid.edges[c.edge].curve is BrepIntersectionCurve ||
            solid.edges[c.edge].curve is BrepEllipse,
      ),
    )) {
      final q = cyl.parameters(infinite)!;
      var inside = false;
      double nearestTrim = 1e-2;
      bool? nearInside;
      for (var li = 0; li < face.loops.length; li++) {
        final loop = face.loops[li];
        final points = BrepTessellate.sampleLoop(solid, loop, 96);
        final poly = <Vector3>[];
        for (final point in points) {
          final uv = cyl.parameters(point)!;
          var u = uv.$1;
          if (poly.isNotEmpty) {
            while (u - poly.last.x > math.pi) {
              u -= 2 * math.pi;
            }
            while (u - poly.last.x < -math.pi) {
              u += 2 * math.pi;
            }
          }
          poly.add(Vector3(u, uv.$2, 0));
        }
        if (poly.isEmpty) continue;
        final center =
            poly.fold(0.0, (double sum, p) => sum + p.x) / poly.length;
        var u = q.$1;
        while (u - center > math.pi) {
          u -= 2 * math.pi;
        }
        while (u - center < -math.pi) {
          u += 2 * math.pi;
        }
        if (_inPoly(poly, Vector3(u, q.$2, 0), Vector3(0, 0, 1)))
          inside = !inside;
        // A sampled trim can miss a ray hit just inside an elliptical seam.
        // Near a smooth trim use its exact local side; retain the polygon
        // test away from the edge and at ambiguous trim vertices.
        var area = 0.0;
        for (var i = 0; i < poly.length; i++) {
          final a = poly[i], b = poly[(i + 1) % poly.length];
          area += a.x * b.y - b.x * a.y;
        }
        if (area.abs() < 1e-12) continue;
        for (final co in loop.coedges) {
          final c = solid.edges[co.edge].curve;
          final t = c.closestParam(infinite), point = c.evaluate(t);
          final distance = (point - infinite).length;
          if (distance >= nearestTrim) continue;
          final lo = math.min(c.t0, c.t1), hi = math.max(c.t0, c.t1);
          if (t <= lo + (hi - lo) * 1e-7 || t >= hi - (hi - lo) * 1e-7)
            continue;
          var tangent = c.tangent(t);
          if ((c.t1 < c.t0) != co.reversed) tangent = -tangent;
          if ((area < 0) != (li > 0)) tangent = -tangent;
          final uv = cyl.parameters(point)!;
          final inward = cyl.normalAt(uv.$1, uv.$2).cross(tangent);
          nearestTrim = distance;
          nearInside = (infinite - point).dot(inward) >= -1e-9;
        }
      }
      if (nearInside ?? inside) return infinite;
      Vector3? best;
      var distance = double.infinity;
      for (final loop in face.loops) {
        for (final co in loop.coedges) {
          final point = solid.edges[co.edge].curve.closestPoint(p);
          final d = (point - p).length2;
          if (d < distance) {
            best = point;
            distance = d;
          }
        }
      }
      if (best != null) return best;
    }
    final domain = BrepTessellate.cylinderFaceDomain(solid, face);
    var u0 = domain.$1;
    var u1 = domain.$2;
    final z0 = domain.$3;
    final z1 = domain.$4;
    if (u1 <= u0 + 1e-9) u1 += math.pi * 2;
    final fullU = (u1 - u0) >= math.pi * 2 - 1e-3;
    final uv = cyl.parameters(infinite) ?? cyl.parameters(p);
    if (uv == null) return infinite;

    bool inU(double u) {
      if (fullU) return true;
      var x = u;
      while (x < u0 - 1e-9) {
        x += math.pi * 2;
      }
      while (x > u0 + math.pi * 2 + 1e-9) {
        x -= math.pi * 2;
      }
      return x >= u0 - 1e-9 && x <= u1 + 1e-9;
    }

    if (inU(uv.$1)) {
      return cyl.evaluate(uv.$1, uv.$2.clamp(z0, z1));
    }

    var best = cyl.evaluate(u0, z0);
    var bestD = (p - best).length2;
    void consider(Vector3 c) {
      final d = (p - c).length2;
      if (d < bestD) {
        bestD = d;
        best = c;
      }
    }

    consider(_closestOnSeg(p, cyl.evaluate(u0, z0), cyl.evaluate(u0, z1)));
    consider(_closestOnSeg(p, cyl.evaluate(u1, z0), cyl.evaluate(u1, z1)));
    consider(_closestOnCylArc(p, cyl, u0, u1, z0));
    consider(_closestOnCylArc(p, cyl, u0, u1, z1));
    return best;
  }

  static Vector3 _closestOnSeg(Vector3 p, Vector3 a, Vector3 b) {
    final d = b - a;
    final len2 = d.length2;
    if (len2 < 1e-20) return a.clone();
    final t = ((p - a).dot(d) / len2).clamp(0.0, 1.0);
    return a + d * t;
  }

  static Vector3 _closestOnCylArc(
    Vector3 p,
    BrepCylinder cyl,
    double u0,
    double u1,
    double z,
  ) {
    final q = cyl.closestPoint(p);
    final uv = cyl.parameters(q);
    if (uv != null) {
      var u = uv.$1;
      while (u < u0) {
        u += math.pi * 2;
      }
      while (u > u0 + math.pi * 2) {
        u -= math.pi * 2;
      }
      if (u >= u0 && u <= u1) return cyl.evaluate(u, z);
    }
    final a = cyl.evaluate(u0, z);
    final b = cyl.evaluate(u1, z);
    return (p - a).length2 <= (p - b).length2 ? a : b;
  }

  /// Even-odd (concave profiles). On-edge counts as inside.
  static bool _inPoly(List<Vector3> poly, Vector3 p, Vector3 n) {
    if (poly.length < 3 || n.length2 < 1e-20) return false;
    final nn = n.normalized();
    if (_distToEdges(p, poly) < 1e-6) return true;
    var xAxis = nn.cross(Vector3(0, 0, 1));
    if (xAxis.length2 < 1e-12) xAxis = nn.cross(Vector3(0, 1, 0));
    xAxis.normalize();
    final yAxis = nn.cross(xAxis);
    final o = poly.first;
    double xd(Vector3 q) => (q - o).dot(xAxis);
    double yd(Vector3 q) => (q - o).dot(yAxis);
    final px = xd(p);
    final py = yd(p);
    var inside = false;
    for (var i = 0, j = poly.length - 1; i < poly.length; j = i++) {
      final yi = yd(poly[i]);
      final yj = yd(poly[j]);
      if ((yi > py) == (yj > py)) continue;
      final xi = xd(poly[i]);
      final xj = xd(poly[j]);
      final xHit = (xj - xi) * (py - yi) / (yj - yi) + xi;
      if (px < xHit) inside = !inside;
    }
    return inside;
  }

  static double _distToEdges(Vector3 p, List<Vector3> poly) {
    var best = double.infinity;
    for (var i = 0; i < poly.length; i++) {
      final a = poly[i];
      final b = poly[(i + 1) % poly.length];
      final d = b - a;
      final len2 = d.length2;
      final t = len2 < 1e-20 ? 0.0 : ((p - a).dot(d) / len2).clamp(0.0, 1.0);
      final dist = (p - (a + d * t)).length;
      if (dist < best) best = dist;
    }
    return best;
  }

  static Vector3 _closestOnPoly(List<Vector3> poly, Vector3 p) {
    var best = poly.first;
    var bestD = (p - best).length2;
    for (var i = 0; i < poly.length; i++) {
      final a = poly[i];
      final b = poly[(i + 1) % poly.length];
      final d = b - a;
      final len2 = d.length2;
      final t = len2 < 1e-20 ? 0.0 : ((p - a).dot(d) / len2).clamp(0.0, 1.0);
      final q = a + d * t;
      final dist = (p - q).length2;
      if (dist < bestD) {
        bestD = dist;
        best = q;
      }
    }
    return best;
  }

  static Vector3 faceCentroid(BrepSolid solid, BrepFace face) {
    if (face.loops.isEmpty) {
      if (face.surface is BrepSphere) {
        return (face.surface as BrepSphere).center;
      }
      return Vector3.zero();
    }
    final pts = BrepTessellate.sampleLoop(solid, face.loops.first, 24);
    var c = Vector3.zero();
    for (final p in pts) {
      c += p;
    }
    if (pts.isEmpty) return Vector3.zero();
    return c / pts.length.toDouble();
  }
}
