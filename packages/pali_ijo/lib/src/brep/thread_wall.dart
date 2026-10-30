import 'dart:math' as math;
import 'package:vector_math/vector_math_64.dart';
import '../thread/form.dart';
import 'geom.dart';
import 'boolean.dart';
import 'solid.dart';
import 'heal.dart';
import 'primitives.dart';

/// Direct boundary construction for complete coaxial cylindrical walls.
/// Work is proportional to the number of turns, without repeated booleans.
/// Parts with interrupted walls continue through the general intersection path.
class BrepThreadWall {
  static double? _angleInInterval(double angle, double a, double b) {
    final lo = math.min(a, b), hi = math.max(a, b);
    final t = angle + ((lo - angle) / (2 * math.pi)).ceil() * 2 * math.pi;
    return t <= hi ? t : null;
  }

  // Conservative support bounds, never tessellation/sample-based clearance.
  // Unknown curves return -infinity and retain the general trimming path.
  static double _support(BrepCurve c, Vector3 origin, Vector3 direction) {
    if (c is BrepLine) {
      return math.min(
        (c.a - origin).dot(direction),
        (c.b - origin).dot(direction),
      );
    }
    if (c is BrepCircle) {
      final x = c.xDir.dot(direction), y = c.yDir.dot(direction);
      final t = _angleInInterval(math.atan2(y, x) + math.pi, c.t0, c.t1);
      final minimum = math.min(
        x * math.cos(c.t0) + y * math.sin(c.t0),
        x * math.cos(c.t1) + y * math.sin(c.t1),
      );
      return (c.center - origin).dot(direction) +
          c.radius * (t == null ? minimum : -math.sqrt(x * x + y * y));
    }
    if (c is BrepHelix) {
      final x = c.xDir.dot(direction), y = c.yDir.dot(direction);
      final axial = c.axis.dot(direction) * c.lead / (2 * math.pi);
      return (c.origin - origin).dot(direction) +
          math.min(axial * c.t0, axial * c.t1) -
          math.max(c.radiusAt(c.t0), c.radiusAt(c.t1)) *
              math.sqrt(x * x + y * y);
    }
    if (c is BrepPlaneSectionCurve &&
        c.t0 >= 0 &&
        c.t0 <= 1 &&
        c.t1 >= 0 &&
        c.t1 <= 1) {
      final a = c.surface.lower, b = c.surface.upper;
      if (a is BrepHelix &&
          b is BrepHelix &&
          c.plane.normal.cross(a.axis).length < 1e-10 &&
          c.plane.normal.cross(b.axis).length < 1e-10) {
        // With an axial cutting plane the ruling is linear-fractional in u.
        // Same-sign denominators bound v by its endpoint values. Bound the
        // linear interpolation too, including any extrapolation beyond rails.
        final d0 = c.surface.dv(c.t0, 0).dot(c.plane.normal);
        final d1 = c.surface.dv(c.t1, 0).dot(c.plane.normal);
        if (d0 * d1 > 0 && d0.abs() > 1e-12 && d1.abs() > 1e-12) {
          final v0 = c.ruling(c.t0), v1 = c.ruling(c.t1);
          final amin = _support(a, origin, direction);
          final amax = -_support(a, origin, -direction);
          final bmin = _support(b, origin, direction);
          final bmax = -_support(b, origin, -direction);
          double bound(double v) =>
              (1 - v) * (v <= 1 ? amin : amax) + v * (v >= 0 ? bmin : bmax);
          return math.min(bound(v0), bound(v1));
        }
      }
    }
    return double.negativeInfinity;
  }

  static Vector3 _radialDirection(Vector3 p, Vector3 origin, Vector3 axis) {
    final delta = p - origin;
    return (delta - axis * delta.dot(axis)).normalized();
  }

  static Vector3 _curveAnchor(BrepCurve c) => switch (c) {
    BrepHelix() => c.origin,
    BrepCircle() => c.center,
    BrepPlaneSectionCurve() => _curveAnchor(c.surface.lower),
    _ => c.evaluate((c.t0 + c.t1) / 2),
  };

  static bool _curveOutsideBore(
    BrepCurve c,
    Vector3 origin,
    Vector3 axis,
    double radius,
  ) {
    if (c is BrepCircle && c.normal.cross(axis).length < 1e-10) {
      // The closest point on a trimmed arc is an endpoint or the radial
      // stationary point. Its unused supporting circle is not an obstruction.
      final delta = origin - c.center;
      final t = _angleInInterval(
        math.atan2(delta.dot(c.yDir), delta.dot(c.xDir)),
        c.t0,
        c.t1,
      );
      return [c.t0, c.t1, if (t != null) t].every(
        (u) => (c.evaluate(u) - origin).cross(axis).length > radius + 1e-7,
      );
    }
    final direction = _radialDirection(_curveAnchor(c), origin, axis);
    return _support(c, origin, direction) > radius + 1e-7;
  }

  static bool _outsideBore(
    BrepSolid body,
    BrepFace face,
    Vector3 origin,
    Vector3 axis,
    double radius,
  ) {
    final surface = face.surface;
    if (surface is BrepCylinder && surface.axis.cross(axis).length < 1e-10) {
      // For a bounded parallel-axis cylindrical patch, every angular station reaches
      // its boundary. Circular rims and axial seams therefore bound the entire
      // wall, including split/partial cylinders, without sampling its surface.
      final edges = face.loops.expand((l) => l.coedges).toList();
      return edges.any((co) => body.edges[co.edge].curve is BrepCircle) &&
          edges.every((co) {
            final c = body.edges[co.edge].curve;
            if (c is BrepCircle &&
                c.normal.cross(axis).length < 1e-10 &&
                (c.center - surface.origin).cross(axis).length < 1e-7 &&
                (c.radius - surface.radius).abs() < 1e-7) {
              return _curveOutsideBore(c, origin, axis, radius);
            }
            if (c is BrepLine && (c.b - c.a).cross(axis).length < 1e-10) {
              return (c.a - origin).cross(axis).length > radius + 1e-7;
            }
            return false;
          });
    }
    if (surface is BrepRuledSurface) {
      final direction = _radialDirection(
        (_curveAnchor(surface.lower) + _curveAnchor(surface.upper)) / 2,
        origin,
        axis,
      );
      return math.min(
            _support(surface.lower, origin, direction),
            _support(surface.upper, origin, direction),
          ) >
          radius + 1e-7;
    }
    if (surface is BrepPlane &&
        face.loops.isNotEmpty &&
        face.loops.first.coedges.isNotEmpty) {
      final first = body.edges[face.loops.first.coedges.first.edge].curve;
      final direction = _radialDirection(_curveAnchor(first), origin, axis);
      // A single separating half-plane for the whole planar boundary also
      // separates its interior. Testing each edge independently would not.
      return face.loops
          .expand((l) => l.coedges)
          .every(
            (co) =>
                _support(body.edges[co.edge].curve, origin, direction) >
                radius + 1e-7,
          );
    }
    return false;
  }

  /// Replace a complete external cylindrical wall on a larger planar part.
  /// Adjacent cap/shoulder loops inherit the new helical rim; the head and its
  /// other faces are retained exactly. Interrupted walls use the general CSG.
  static BrepSolid? tryReplaceCylinder(
    BrepSolid body, {
    required Vector3 origin,
    required Vector3 axis,
    required double radius,
    required double pitch,
    required double lo,
    required double hi,
    required List<(double, double)> regions,
    required bool rightHand,
    required int starts,
  }) {
    final walls = body.faces
        .where(
          (f) =>
              f.surface is BrepCylinder &&
              ((f.surface as BrepCylinder).radius - radius).abs() < 1e-7 &&
              (f.surface as BrepCylinder).axis.cross(axis).length < 1e-8 &&
              ((f.surface as BrepCylinder).origin - origin).cross(axis).length <
                  1e-7,
        )
        .toList();
    if (walls.length != 1 ||
        body.faces.any(
          (f) => f.surface is! BrepPlane && f.surface is! BrepCylinder,
        ))
      return null;
    final wall = walls.single, cylinder = wall.surface as BrepCylinder;
    final n = axis.normalized();
    final rootRadius = radius - ThreadForm(pitch: pitch, internal: false).depth;
    bool insideRoot(BrepCurve c) {
      if (c is BrepCircle)
        return c.normal.cross(n).length < 1e-8 &&
            (c.center - origin).cross(n).length + c.radius < rootRadius - 1e-7;
      if (c is BrepLine)
        return [
          c.a,
          c.b,
        ].every((p) => (p - origin).cross(n).length < rootRadius - 1e-7);
      return false;
    }

    if (wall.reversed ||
        (cylinder.radius - radius).abs() > 1e-7 ||
        cylinder.axis.cross(n).length > 1e-8 ||
        (cylinder.origin - origin).cross(n).length > 1e-7)
      return null;
    final wallEdges = {
      for (final l in wall.loops)
        for (final c in l.coedges) c.edge,
    };
    final rims = <double, Set<int>>{lo: {}, hi: {}};
    for (final id in wallEdges) {
      final c = body.edges[id].curve;
      if (c is BrepLine) {
        if ((c.evaluate(c.t1) - c.evaluate(c.t0)).cross(n).length > 1e-7)
          return null;
        continue;
      }
      if (c is! BrepCircle ||
          (c.radius - radius).abs() > 1e-7 ||
          c.normal.cross(n).length > 1e-8)
        return null;
      final z = (c.center - origin).dot(n);
      final end = (z - lo).abs() < 1e-7
          ? lo
          : (z - hi).abs() < 1e-7
          ? hi
          : null;
      if (end == null) return null;
      rims[end]!.add(id);
    }
    for (final rim in rims.values) {
      final angle = rim.fold(0.0, (v, id) {
        final c = body.edges[id].curve as BrepCircle;
        return v + (c.t1 - c.t0).abs();
      });
      if ((angle - 2 * math.pi).abs() > 1e-7) return null;
    }
    final replacements = <BrepFace, (int, double)>{};
    for (final f in body.faces) {
      if (identical(f, wall)) continue;
      if (f.surface is BrepCylinder) {
        final c = f.surface as BrepCylinder;
        if (c.axis.cross(n).length > 1e-8 ||
            (c.origin - origin).cross(n).length + c.radius >= rootRadius - 1e-7)
          return null;
        continue;
      }
      final plane = f.surface as BrepPlane;
      final z = (plane.origin - origin).dot(n);
      if (plane.normal.cross(n).length < 1e-8) {
        final end = (z - lo).abs() < 1e-7
            ? lo
            : (z - hi).abs() < 1e-7
            ? hi
            : null;
        if (end == null) {
          if (z > lo && z < hi) return null;
          continue;
        }
        final ids = rims[end]!;
        final index = f.loops.indexWhere(
          (l) =>
              l.coedges.length == ids.length &&
              l.coedges.every((c) => ids.contains(c.edge)),
        );
        if (index < 0) {
          if (f.loops.every(
            (l) =>
                l.coedges.every((co) => insideRoot(body.edges[co.edge].curve)),
          ))
            continue;
          return null;
        }
        // Additional boundaries must stay outside the cylinder's envelope.
        // This deliberately excludes holes or pockets intersecting the thread.
        for (var i = 0; i < f.loops.length; i++) {
          if (i == index) continue;
          for (final co in f.loops[i].coedges) {
            final c = body.edges[co.edge].curve;
            if (insideRoot(c)) continue;
            if (c is! BrepLine) return null;
            final a = c.evaluate(c.t0) - origin;
            final b = c.evaluate(c.t1) - origin;
            final p = a - n * a.dot(n), q = b - n * b.dot(n), d = q - p;
            final t = d.length2 == 0
                ? 0.0
                : (-p.dot(d) / d.length2).clamp(0.0, 1.0);
            if ((p + d * t).length < radius - 1e-7) return null;
          }
        }
        replacements[f] = (index, end);
      } else if (plane.normal.dot(n).abs() > 1e-8 ||
          (origin - plane.origin).dot(plane.normal).abs() < radius - 1e-7) {
        return null;
      }
    }
    if (replacements.length != 2 ||
        replacements.values.map((v) => v.$2).toSet().length != 2)
      return null;
    final patch = tryBuild(
      BrepPrimitives.cylinder(
        origin: origin + n * lo,
        axis: n,
        radius: radius,
        height: hi - lo,
      ),
      origin: origin,
      axis: n,
      radius: radius,
      pitch: pitch,
      lo: lo,
      hi: hi,
      regions: regions,
      internal: false,
      rightHand: rightHand,
      starts: starts,
    );
    if (patch == null) return null;
    final caps = <double, BrepFace>{};
    for (final f in patch.faces) {
      final s = f.surface;
      if (s is! BrepPlane) continue;
      final z = (s.origin - origin).dot(n);
      if ((z - lo).abs() < 1e-7) caps[lo] = f;
      if ((z - hi).abs() < 1e-7) caps[hi] = f;
    }
    if (caps.length != 2 || caps.values.any((f) => f.loops.length != 1))
      return null;
    final joined = BrepSolid.compound([body, patch]);
    final offset = body.edges.length;
    BrepLoop remap(BrepLoop l, bool reverse) => BrepLoop([
      for (final c in reverse ? l.coedges.reversed : l.coedges)
        BrepCoedge(c.edge + offset, reversed: c.reversed != reverse),
    ]);
    final faces = <BrepFace>[];
    for (final f in body.faces) {
      if (identical(f, wall)) continue;
      final replacement = replacements[f];
      final loops = List<BrepLoop>.from(f.loops);
      if (replacement != null) {
        final cap = caps[replacement.$2]!;
        final reverse = f.outwardNormal(0, 0).dot(cap.outwardNormal(0, 0)) < 0;
        loops[replacement.$1] = remap(cap.loops.single, reverse);
      }
      faces.add(
        BrepFace(
          surface: f.surface,
          loops: loops,
          reversed: f.reversed,
          origin: f.origin,
        ),
      );
    }
    for (var i = 0; i < patch.faces.length; i++) {
      if (!caps.values.contains(patch.faces[i]))
        faces.add(joined.faces[body.faces.length + i]);
    }
    final result = BrepSolid(
      vertices: joined.vertices,
      edges: joined.edges,
      faces: faces,
    );
    return BrepHeal.topologyIsClosed(result) ? result : null;
  }

  /// Replace a complete bore in a planar part without intersecting every
  /// helical face against the body. Clearance checks reject interrupted bores.
  static BrepSolid? tryReplaceBore(
    BrepSolid body, {
    required Vector3 origin,
    required Vector3 axis,
    required double radius,
    required double pitch,
    required double lo,
    required double hi,
    required List<(double, double)> regions,
    required bool rightHand,
    required int starts,
  }) {
    final majorRadius = radius + ThreadForm(pitch: pitch, internal: true).depth;
    final walls = body.faces
        .where(
          (f) =>
              f.surface is BrepCylinder &&
              f.reversed &&
              (f.surface as BrepCylinder).radius >= radius - 1e-7 &&
              (f.surface as BrepCylinder).radius <= majorRadius + 1e-7 &&
              (f.surface as BrepCylinder).axis.cross(axis).length < 1e-8 &&
              ((f.surface as BrepCylinder).origin - origin).cross(axis).length <
                  1e-7,
        )
        .toList();
    if (walls.isEmpty) return null;
    final cylinder = walls.first.surface as BrepCylinder;
    final n = axis.normalized();
    final boreRadius = cylinder.radius;
    if (walls.any(
      (f) => ((f.surface as BrepCylinder).radius - boreRadius).abs() > 1e-7,
    )) {
      return null;
    }
    // A boolean may split a complete bore at its periodic seam. Treat all
    // matching patches as one wall, requiring paired longitudinal seams and
    // one complete circular boundary at each end.
    final wallEdges = <int, int>{};
    for (final wall in walls) {
      for (final co in wall.loops.expand((l) => l.coedges)) {
        wallEdges.update(co.edge, (n) => n + 1, ifAbsent: () => 1);
      }
    }
    final rims = <double, Set<int>>{lo: {}, hi: {}};
    for (final id in wallEdges.keys) {
      final c = body.edges[id].curve;
      if (c is BrepLine) {
        if (wallEdges[id] != 2 ||
            (c.evaluate(c.t1) - c.evaluate(c.t0)).cross(n).length > 1e-7)
          return null;
        continue;
      }
      if (wallEdges[id] != 1 ||
          c is! BrepCircle ||
          (c.radius - boreRadius).abs() > 1e-7 ||
          c.normal.cross(n).length > 1e-8)
        return null;
      final z = (c.center - origin).dot(n);
      final end = (z - lo).abs() < 1e-7
          ? lo
          : (z - hi).abs() < 1e-7
          ? hi
          : null;
      if (end == null) return null;
      rims[end]!.add(id);
    }
    for (final rim in rims.values) {
      final angle = rim.fold(0.0, (v, id) {
        final c = body.edges[id].curve as BrepCircle;
        return v + (c.t1 - c.t0).abs();
      });
      if ((angle - 2 * math.pi).abs() > 1e-7) return null;
    }
    final replacements = <BrepFace, (int, double)>{};
    for (final f in body.faces) {
      if (walls.contains(f)) continue;
      // Existing threads and other separated patches must not force a global
      // helical boolean. Only skip surfaces with a proven separating bound.
      if (_outsideBore(body, f, origin, n, majorRadius)) continue;
      if (f.surface is! BrepPlane && f.surface is! BrepCylinder) return null;
      if (f.surface is BrepCylinder) {
        final c = f.surface as BrepCylinder;
        // Only accept other cylinders entirely outside the swept envelope.
        // Intersecting or coaxial stepped bores need general trimming.
        final distance = (c.origin - origin).cross(n).length;
        if (c.axis.cross(n).length > 1e-8 ||
            (distance - c.radius <= majorRadius + 1e-7 &&
                c.radius - distance <= majorRadius + 1e-7))
          return null;
        continue;
      }
      final plane = f.surface as BrepPlane;
      final z = (plane.origin - origin).dot(n);
      if (plane.normal.cross(n).length < 1e-8) {
        final end = (z - lo).abs() < 1e-7
            ? lo
            : (z - hi).abs() < 1e-7
            ? hi
            : null;
        if (end == null) {
          if (z > lo && z < hi) return null;
          continue;
        }
        final ids = rims[end]!;
        final index = f.loops.indexWhere(
          (l) =>
              l.coedges.length == ids.length &&
              l.coedges.every((c) => ids.contains(c.edge)),
        );
        if (index < 0) return null;
        // Additional boundaries must stay outside the cylinder's envelope.
        // This deliberately excludes holes or pockets intersecting the thread.
        for (var i = 0; i < f.loops.length; i++) {
          if (i == index) continue;
          for (final co in f.loops[i].coedges) {
            final c = body.edges[co.edge].curve;
            if (c is BrepCircle) {
              if (c.normal.cross(n).length > 1e-8 ||
                  !_curveOutsideBore(c, origin, n, majorRadius))
                return null;
              continue;
            }
            if (c is! BrepLine) {
              if (!_curveOutsideBore(c, origin, n, majorRadius)) return null;
              continue;
            }
            final a = c.evaluate(c.t0) - origin;
            final b = c.evaluate(c.t1) - origin;
            final p = a - n * a.dot(n), q = b - n * b.dot(n), d = q - p;
            final t = d.length2 == 0
                ? 0.0
                : (-p.dot(d) / d.length2).clamp(0.0, 1.0);
            if ((p + d * t).length <= majorRadius + 1e-7) return null;
          }
        }
        replacements[f] = (index, end);
      } else if (plane.normal.dot(n).abs() > 1e-8 ||
          (origin - plane.origin).dot(plane.normal).abs() <=
              majorRadius + 1e-7) {
        return null;
      }
    }
    if (replacements.length != 2 ||
        replacements.values.map((v) => v.$2).toSet().length != 2)
      return null;
    // A temporary coaxial tube supplies a closed helical bore and its rims.
    // Its outer wall is discarded; the actual body's exterior is retained.
    final tube = BrepBoolean.run(
      BrepPrimitives.cylinder(
        origin: origin + n * lo,
        axis: n,
        radius: majorRadius + pitch,
        height: hi - lo,
      ),
      BrepPrimitives.cylinder(
        origin: origin + n * (lo - pitch),
        axis: n,
        radius: boreRadius,
        height: hi - lo + 2 * pitch,
      ),
      BrepBoolean.subtractOp,
    );
    final patch = tryBuild(
      tube,
      origin: origin,
      axis: n,
      radius: radius,
      pitch: pitch,
      lo: lo,
      hi: hi,
      regions: regions,
      internal: true,
      rightHand: rightHand,
      starts: starts,
    );
    if (patch == null) return null;
    final caps = <double, BrepFace>{};
    for (final f in patch.faces) {
      final s = f.surface;
      if (s is! BrepPlane) continue;
      final z = (s.origin - origin).dot(n);
      if ((z - lo).abs() < 1e-7) caps[lo] = f;
      if ((z - hi).abs() < 1e-7) caps[hi] = f;
    }
    if (caps.length != 2 || caps.values.any((f) => f.loops.length != 2))
      return null;
    final joined = BrepSolid.compound([body, patch]);
    final offset = body.edges.length;
    BrepLoop remap(BrepLoop l, bool reverse) => BrepLoop([
      for (final c in reverse ? l.coedges.reversed : l.coedges)
        BrepCoedge(c.edge + offset, reversed: c.reversed != reverse),
    ]);
    final faces = <BrepFace>[];
    for (final f in body.faces) {
      if (walls.contains(f)) continue;
      final replacement = replacements[f];
      final loops = List<BrepLoop>.from(f.loops);
      if (replacement != null) {
        final cap = caps[replacement.$2]!;
        final reverse = f.outwardNormal(0, 0).dot(cap.outwardNormal(0, 0)) < 0;
        loops[replacement.$1] = remap(cap.loops[1], reverse);
      }
      faces.add(
        BrepFace(
          surface: f.surface,
          loops: loops,
          reversed: f.reversed,
          origin: f.origin,
        ),
      );
    }
    for (var i = 0; i < patch.faces.length; i++) {
      final f = patch.faces[i];
      final surface = f.surface;
      final outerWall =
          surface is BrepCylinder &&
          (surface.radius - (majorRadius + pitch)).abs() < 1e-7;
      if (!caps.values.contains(f) && !outerWall)
        faces.add(joined.faces[body.faces.length + i]);
    }
    final result = BrepSolid(
      vertices: joined.vertices,
      edges: joined.edges,
      faces: faces,
    );
    return BrepHeal.topologyIsClosed(result) ? result : null;
  }

  static BrepSolid? tryBuild(
    BrepSolid body, {
    required Vector3 origin,
    required Vector3 axis,
    required double radius,
    required double pitch,
    required double lo,
    required double hi,
    required List<(double, double)> regions,
    required bool internal,
    required bool rightHand,
    required int starts,
  }) {
    final cylinders = body.faces
        .where((f) => f.surface is BrepCylinder)
        .toList();
    final planes = body.faces.where((f) => f.surface is BrepPlane).toList();
    if (planes.length != 2 ||
        cylinders.isEmpty ||
        cylinders.length > 2 ||
        planes.length + cylinders.length != body.faces.length ||
        body.edges.any((e) => e.curve is! BrepLine && e.curve is! BrepCircle))
      return null;
    final n = axis.normalized(), x = brepFrameU(n), y = n.cross(x);
    final caps = <double>[];
    for (final face in planes) {
      final p = face.surface as BrepPlane;
      if (p.normal.cross(n).length > 1e-8 ||
          face.loops.any(
            (l) =>
                l.coedges.any((c) => body.edges[c.edge].curve is! BrepCircle),
          ))
        return null;
      caps.add((p.origin - origin).dot(n));
    }
    caps.sort();
    if ((caps.first - lo).abs() > 1e-7 || (caps.last - hi).abs() > 1e-7)
      return null;
    final radii = <double>[];
    for (final face in cylinders) {
      final c = face.surface as BrepCylinder;
      if (c.axis.cross(n).length > 1e-8 ||
          (c.origin - origin).cross(n).length > 1e-7)
        return null;
      radii.add(c.radius);
    }
    radii.sort();
    final outer = radii.last, inner = radii.length == 2 ? radii.first : 0.0;
    if ((!internal && (radius - outer).abs() > 1e-7) ||
        internal && (inner == 0 || radius > inner + 1e-7))
      return null;
    final form = ThreadForm(pitch: pitch, internal: internal);
    if (internal && inner > radius + form.depth + 1e-7) return null;
    if ((hi - lo) / pitch + starts > 800)
      throw BrepUnsupported('Thread sweep exceeds the current face budget');
    // A nominal bore is a modeling convention: form the female crests inward
    // from its major diameter. A smaller, predrilled bore still cuts outward.
    final nominalBore = internal && (inner - radius - form.depth).abs() <= 1e-7;
    final wallRadius = internal && !nominalBore ? inner : radius;
    final plainRadius = internal ? inner : radius;
    final depthAmount = internal && !nominalBore
        ? radius + form.depth - inner
        : form.depth;
    final mouth = form.rootHalfWidth + depthAmount / math.sqrt(3);
    if (depthAmount <= 1e-7) return null;
    if (internal && radius + form.depth >= outer - 1e-7 ||
        !internal && outer - form.depth <= inner + 1e-7)
      return null;
    final b = BrepBuilder();
    final uses = <int, List<BrepCoedge>>{}, byPair = <(int, int), List<int>>{};
    int vertex(Vector3 p) {
      for (var i = 0; i < b.vertices.length; i++) {
        if ((b.vertices[i] - p).length < 1e-7) return i;
      }
      return b.vertex(p);
    }

    BrepCoedge edge(BrepCurve c) {
      if (c is BrepPlaneSectionCurve) {
        final a = c.surface.lower, z = c.surface.upper;
        if (a is BrepHelix &&
            z is BrepHelix &&
            a.radiusSlope == 0 &&
            z.radiusSlope == 0 &&
            (a.radius - z.radius).abs() < 1e-10 &&
            c.plane.normal.cross(n).length < 1e-10) {
          // Constant-radius axial sections are exact circular arcs. Use the
          // same representation as the adjacent unthreaded wall.
          c = BrepCircle(
            center: origin + n * (c.plane.origin - origin).dot(n),
            normal: n,
            xDir: a.xDir,
            radius: a.radius,
            t0: a.t0 + (a.t1 - a.t0) * c.t0,
            t1: a.t0 + (a.t1 - a.t0) * c.t1,
          );
        }
      }
      final backward = c.t1 < c.t0;
      if (backward) {
        c = switch (c) {
          BrepHelix() => c.trimmed(c.t1, c.t0),
          BrepPlaneSectionCurve() => c.trimmed(c.t1, c.t0),
          BrepCircle() => BrepCircle(
            center: c.center,
            normal: c.normal,
            xDir: c.xDir,
            radius: c.radius,
            t0: c.t1,
            t1: c.t0,
          ),
          _ => throw BrepUnsupported(
            'Thread edge has an unsupported parameter direction',
          ),
        };
      }
      final i = vertex(c.evaluate(c.t0)), j = vertex(c.evaluate(c.t1));
      final candidates = byPair.putIfAbsent((
        math.min(i, j),
        math.max(i, j),
      ), () => []);
      for (final id in candidates) {
        final old = b.edges[id], reverse = old.v0 != i;
        if ([.25, .5, .75].every(
          (t) =>
              (old.curve.evaluate(
                        old.curve.t0 +
                            (old.curve.t1 - old.curve.t0) *
                                (reverse ? 1 - t : t),
                      ) -
                      c.evaluate(c.t0 + (c.t1 - c.t0) * t))
                  .length <
              1e-7,
        ))
          return BrepCoedge(id, reversed: reverse != backward);
      }
      final id = b.edge(i, j, c);
      candidates.add(id);
      return BrepCoedge(id, reversed: backward);
    }

    void face(BrepSurface s, List<BrepCoedge> loop, {bool reversed = false}) {
      final oriented = reversed
          ? [
              for (final c in loop.reversed)
                BrepCoedge(c.edge, reversed: !c.reversed),
            ]
          : loop;
      for (final c in oriented) (uses[c.edge] ??= []).add(c);
      b.face(s, [BrepLoop(oriented)], reversed: reversed);
    }

    List<Vector2> clip(List<Vector2> polygon, double Function(Vector2) signed) {
      final out = <Vector2>[];
      for (var i = 0; i < polygon.length; i++) {
        final a = polygon[i],
            d = polygon[(i + 1) % polygon.length],
            fa = signed(a),
            fb = signed(d);
        if (fa >= -1e-12) out.add(a);
        if ((fa < -1e-12 && fb > 1e-12) || (fb < -1e-12 && fa > 1e-12))
          out.add(a + (d - a) * (fa / (fa - fb)));
      }
      return out;
    }

    final lead = pitch * starts * (rightHand ? 1 : -1),
        h = lead / (2 * math.pi);
    for (final region in regions) {
      for (var sector = 0; sector < 4; sector++) {
        final a = sector * math.pi / 2, d = (sector + 1) * math.pi / 2;
        final qlo = region.$1 - math.max(h * a, h * d) - lo,
            qhi = region.$2 - math.min(h * a, h * d) - lo;
        final first = (qlo / pitch - .5).floor(),
            last = (qhi / pitch + .5).ceil();
        for (var turn = first; turn <= last; turn++) {
          final q = [
            -pitch / 2,
            -mouth,
            -form.rootHalfWidth,
            form.rootHalfWidth,
            mouth,
            pitch / 2,
          ];
          final depth = [0.0, 0.0, depthAmount, depthAmount, 0.0, 0.0];
          for (var band = 0; band < 5; band++) {
            final za = lo + turn * pitch + q[band],
                zb = lo + turn * pitch + q[band + 1];
            final ra = wallRadius + (internal ? 1 : -1) * depth[band],
                rb = wallRadius + (internal ? 1 : -1) * depth[band + 1];
            final lower = BrepHelix(
              origin: origin + n * za,
              axis: n,
              xDir: x,
              radius: ra,
              lead: lead,
              t0: a,
              t1: d,
            );
            final upper = BrepHelix(
              origin: origin + n * zb,
              axis: n,
              xDir: x,
              radius: rb,
              lead: lead,
              t0: a,
              t1: d,
            );
            final s = BrepRuledSurface(lower, upper);
            double axial(Vector2 p) =>
                h * (a + (d - a) * p.x) + za + (zb - za) * p.y;
            var poly = [
              Vector2(0, 0),
              Vector2(1, 0),
              Vector2(1, 1),
              Vector2(0, 1),
            ];
            poly = clip(poly, (p) => axial(p) - region.$1);
            poly = clip(poly, (p) => region.$2 - axial(p));
            if (poly.length < 3) continue;
            final co = <BrepCoedge>[];
            for (var k = 0; k < poly.length; k++) {
              final p = poly[k], r = poly[(k + 1) % poly.length];
              if ((p - r).length < 1e-10) continue;
              BrepCurve c;
              if ((p.y - r.y).abs() < 1e-9 &&
                  (p.y.abs() < 1e-9 || (p.y - 1).abs() < 1e-9)) {
                c = (p.y < .5 ? lower : upper).trimmed(
                  a + (d - a) * p.x,
                  a + (d - a) * r.x,
                );
              } else if ((p.x - r.x).abs() < 1e-9) {
                c = BrepLine(s.evaluate(p.x, p.y), s.evaluate(r.x, r.y));
              } else {
                c = BrepPlaneSectionCurve(
                  surface: s,
                  plane: BrepPlane(origin: origin + n * axial(p), normal: n),
                  t0: p.x,
                  t1: r.x,
                );
              }
              if ((c.evaluate(c.t0) - c.evaluate(c.t1)).length > 1e-8)
                co.add(edge(c));
            }
            if (co.length >= 3) face(s, co, reversed: internal);
          }
        }
      }
    }
    // Split the remaining circular walls at every existing endpoint. This
    // makes the thread/plain transitions share real edges instead of T-joints.
    void plain(double r, double low, double high, bool inward) {
      final cuts = <double>[
        0,
        math.pi / 2,
        math.pi,
        3 * math.pi / 2,
        2 * math.pi,
      ];
      for (final p in b.vertices) {
        final delta = p - origin, z = delta.dot(n);
        if ((z - low).abs() > 1e-7 && (z - high).abs() > 1e-7) continue;
        // Include the angles of the opposite rim too. A partial thread has
        // a thin planar transition between different radii; a long display
        // chord on the plain rim can otherwise cross the finer inner rail.
        // Splitting the exact circle at those angles keeps the cap ordered
        // even at draft quality, without changing the native surface.
        var t = math.atan2(delta.dot(y), delta.dot(x));
        if (t < 0) t += 2 * math.pi;
        if (t > 1e-8 && t < 2 * math.pi - 1e-8) cuts.add(t);
      }
      cuts.sort();
      final angles = <double>[];
      for (final t in cuts)
        if (angles.isEmpty || t - angles.last > 1e-8) angles.add(t);
      final cylinder = BrepCylinder(
        origin: origin,
        axis: n,
        xDir: x,
        radius: r,
      );
      for (var i = 1; i < angles.length; i++) {
        final a = angles[i - 1], d = angles[i];
        face(cylinder, [
          edge(
            BrepCircle(
              center: origin + n * low,
              normal: n,
              xDir: x,
              radius: r,
              t0: a,
              t1: d,
            ),
          ),
          edge(BrepLine(cylinder.evaluate(d, low), cylinder.evaluate(d, high))),
          edge(
            BrepCircle(
              center: origin + n * high,
              normal: n,
              xDir: x,
              radius: r,
              t0: d,
              t1: a,
            ),
          ),
          edge(BrepLine(cylinder.evaluate(a, high), cylinder.evaluate(a, low))),
        ], reversed: inward);
      }
    }

    var at = lo;
    for (final region in regions) {
      if (region.$1 - at > 1e-8) plain(plainRadius, at, region.$1, internal);
      at = region.$2;
    }
    if (hi - at > 1e-8) plain(plainRadius, at, hi, internal);
    if (internal) plain(outer, lo, hi, false);
    if (!internal && inner > 0) plain(inner, lo, hi, true);
    // Opposite ends can introduce different angular splits on the same rim.
    // Sew coincident circular intervals before walking the remaining cap
    // boundaries; otherwise they create zero-width overlapping end faces.
    final circularSplits = <int, List<BrepCoedge>>{};
    final shared = <(int, int), List<int>>{};
    final oldEdges = List<BrepEdge>.from(b.edges);
    for (var id = 0; id < oldEdges.length; id++) {
      final old = oldEdges[id], c = old.curve;
      if (c is! BrepCircle) continue;
      final cuts = <(double, int)>[(c.t0, old.v0), (c.t1, old.v1)];
      for (var v = 0; v < b.vertices.length; v++) {
        if (v == old.v0 || v == old.v1) continue;
        final delta = b.vertices[v] - c.center;
        if (delta.dot(c.normal).abs() > 1e-7 ||
            (delta.length - c.radius).abs() > 1e-7)
          continue;
        final t = _angleInInterval(
          math.atan2(delta.dot(c.yDir), delta.dot(c.xDir)),
          c.t0,
          c.t1,
        );
        if (t != null && t > c.t0 + 1e-8 && t < c.t1 - 1e-8) cuts.add((t, v));
      }
      cuts.sort((a, z) => a.$1.compareTo(z.$1));
      final chain = <BrepCoedge>[];
      for (var i = 1; i < cuts.length; i++) {
        final a = cuts[i - 1], z = cuts[i];
        if (z.$1 - a.$1 < 1e-8) continue;
        final curve = BrepCircle(
          center: c.center,
          normal: c.normal,
          xDir: c.xDir,
          radius: c.radius,
          t0: a.$1,
          t1: z.$1,
        );
        final key = (math.min(a.$2, z.$2), math.max(a.$2, z.$2));
        final candidates = shared.putIfAbsent(key, () => []);
        final mid = curve.evaluate((curve.t0 + curve.t1) / 2);
        final found = candidates.where((e) {
          final other = b.edges[e].curve;
          return (other.evaluate((other.t0 + other.t1) / 2) - mid).length <
              1e-7;
        }).firstOrNull;
        final e = found ?? b.edge(a.$2, z.$2, curve);
        if (found == null) candidates.add(e);
        chain.add(BrepCoedge(e, reversed: b.edges[e].v0 != a.$2));
      }
      circularSplits[id] = chain;
    }
    for (var i = 0; i < b.faces.length; i++) {
      final f = b.faces[i];
      b.faces[i] = BrepFace(
        surface: f.surface,
        reversed: f.reversed,
        origin: f.origin,
        loops: [
          for (final l in f.loops)
            BrepLoop([
              for (final co in l.coedges)
                if (!circularSplits.containsKey(co.edge))
                  co
                else ...[
                  for (final split
                      in co.reversed
                          ? circularSplits[co.edge]!.reversed
                          : circularSplits[co.edge]!)
                    BrepCoedge(
                      split.edge,
                      reversed: split.reversed != co.reversed,
                    ),
                ],
            ]),
        ],
      );
    }
    uses.clear();
    for (final f in b.faces) {
      for (final co in f.loops.expand((l) => l.coedges))
        (uses[co.edge] ??= []).add(co);
    }
    final pending = <BrepCoedge>[
      for (final entry in uses.entries)
        if (entry.value.length == 1)
          BrepCoedge(entry.key, reversed: !entry.value.single.reversed),
    ];
    int start(BrepCoedge c) =>
        c.reversed ? b.edges[c.edge].v1 : b.edges[c.edge].v0;
    int end(BrepCoedge c) =>
        c.reversed ? b.edges[c.edge].v0 : b.edges[c.edge].v1;
    final capLoops = <(double, List<BrepLoop>)>[];
    while (pending.isNotEmpty) {
      final loop = <BrepCoedge>[pending.removeLast()],
          first = start(loop.first);
      while (end(loop.last) != first) {
        final i = pending.indexWhere((c) => start(c) == end(loop.last));
        if (i < 0)
          throw BrepUnsupported(
            'Thread boundary could not close an axial transition',
          );
        loop.add(pending.removeAt(i));
      }
      final z = (b.vertices[first] - origin).dot(n);
      if (loop.any(
        (c) => ((b.vertices[start(c)] - origin).dot(n) - z).abs() > 1e-7,
      )) {
        throw BrepUnsupported(
          'Thread boundary has an unpaired longitudinal seam',
        );
      }
      final i = capLoops.indexWhere((c) => (c.$1 - z).abs() < 1e-7);
      if (i < 0)
        capLoops.add((z, [BrepLoop(loop)]));
      else
        capLoops[i].$2.add(BrepLoop(loop));
    }
    for (final cap in capLoops) {
      final z = cap.$1;
      // The pending-edge walk can discover the bore before the outer rim.
      // Planar faces require the enclosing loop first, followed by holes.
      // These caps are coaxial, and the wall-depth guard above guarantees
      // that their radial ranges do not cross.
      double radialExtent(BrepLoop loop) => loop.coedges.fold(0.0, (r, co) {
        final delta = b.vertices[start(co)] - origin;
        return math.max(r, (delta - n * delta.dot(n)).length2);
      });
      cap.$2.sort((a, b) => radialExtent(b).compareTo(radialExtent(a)));
      // Cut grooves face into the removed material. Nominal female threads
      // add crests, so their internal axial transition faces point outward.
      final sign = (z - lo).abs() < 1e-7
          ? -1
          : (z - hi).abs() < 1e-7
          ? 1
          : regions.any((r) => (r.$1 - z).abs() < 1e-7)
          ? (nominalBore ? -1 : 1)
          : (nominalBore ? 1 : -1);
      b.face(
        BrepPlane(origin: origin + n * z, normal: n * sign.toDouble()),
        cap.$2,
      );
    }
    // Separate groove ends can form disjoint planar islands (multi-start
    // threads). They are separate faces, not holes in the first end patch.
    final result = BrepHeal.splitPlanarIslands(b.build());
    if (!BrepHeal.topologyIsClosed(result))
      throw BrepUnsupported('Thread boundary is not a closed manifold');
    return result;
  }
}
