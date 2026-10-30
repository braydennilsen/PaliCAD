import 'dart:math' as math;

import 'package:vector_math/vector_math_64.dart';

import 'geom.dart';
import 'solid.dart';

part 'xform_editing.dart';

class BrepXform {
  BrepXform._();

  /// Expose exact world-space primitives for dimensioned editing. This also
  /// accepts older saved affine wrappers, preserving topology and provenance.
  /// Nonuniform scaling must not be undone before applying millimetre offsets.
  /// Editors that already handle affine spline sweeps can preserve their
  /// boundary/display parameterization with [preserveAffineSweeps].
  static BrepSolid materializeForEditing(
    BrepSolid solid, {
    bool preserveAffineSweeps = false,
  }) => _materializeForEditing(solid, preserveAffineSweeps);

  static BrepSolid translate(BrepSolid s, Vector3 d) => map(s, (p) => p + d);

  /// Place a CAD instance without wrapping its analytic faces in a generic
  /// affine surface. Native plane/cylinder/cone editing remains available.
  static BrepSolid rigid(BrepSolid solid, BrepAffineTransform transform) {
    final x = transform.x, y = transform.y, z = transform.z;
    if ((x.length2 - 1).abs() > 1e-9 ||
        (y.length2 - 1).abs() > 1e-9 ||
        (z.length2 - 1).abs() > 1e-9 ||
        x.dot(y).abs() > 1e-9 ||
        x.dot(z).abs() > 1e-9 ||
        y.dot(z).abs() > 1e-9 ||
        transform.determinant < 0) {
      throw BrepUnsupported('Rigid placement requires a proper rotation');
    }
    return BrepSolid(
      vertices: [for (final p in solid.vertices) transform.point(p)],
      edges: [
        for (final e in solid.edges)
          BrepEdge(
            v0: e.v0,
            v1: e.v1,
            curve: _mapCurve(e.curve, transform.point),
            origin: e.origin,
          ),
      ],
      faces: [
        for (final f in solid.faces)
          BrepFace(
            surface: _mapSurface(f.surface, transform.point),
            loops: f.loops,
            reversed: f.reversed,
            origin: f.origin,
          ),
      ],
    );
  }

  /// Applies [f] to every vertex, curve and surface.
  ///
  /// [reverseWinding] reverses each loop *and* flips each face's `reversed`
  /// flag. [flipFaceNormals] splits the second effect out, because the two are
  /// not always wanted together: after a reflection the loops must reverse to
  /// stay counter-clockwise about the outward normal, but the outward normal
  /// itself simply reflects, so flipping `reversed` as well turns the solid
  /// inside out. See [mirror].
  static BrepSolid map(
    BrepSolid s,
    Vector3 Function(Vector3 p) f, {
    bool reverseWinding = false,
    bool? flipFaceNormals,
  }) {
    final flipNormals = flipFaceNormals ?? reverseWinding;
    if (s.faces.any(
          (face) =>
              face.surface is BrepOffsetSurface ||
              face.surface is BrepRuledSurface ||
              face.surface is BrepSplineSurface,
        ) &&
        !flipNormals) {
      return affine(s, BrepAffineTransform.fromFunction(f));
    }
    final b = BrepBuilder();
    for (final v in s.vertices) {
      b.vertex(f(v));
    }
    // A transform is the one operation whose output faces and edges correspond
    // one-to-one with its input, so provenance carries across by index and
    // needs no geometric matching. This is what lets a face colour survive a
    // move, a scale or a mirror.
    for (final e in s.edges) {
      b.edge(e.v0, e.v1, _mapCurve(e.curve, f), origin: e.origin);
    }
    for (final face in s.faces) {
      final loops = [
        for (final loop in face.loops)
          BrepLoop(
            reverseWinding
                ? [
                    for (final c in loop.coedges.reversed)
                      BrepCoedge(c.edge, reversed: !c.reversed),
                  ]
                : loop.coedges,
          ),
      ];
      b.face(
        _mapSurface(face.surface, f),
        loops,
        reversed: flipNormals ? !face.reversed : face.reversed,
        origin: face.origin,
      );
    }
    return b.build();
  }

  /// Reflects the solid across the plane through [origin] with [normal].
  ///
  /// A reflection reverses handedness, so every loop is reversed to keep it
  /// counter-clockwise about its face's outward normal. The outward normals are
  /// *not* flipped: `_mapSurface` already reflects each stored surface normal,
  /// and reflecting an outward-pointing vector leaves it outward-pointing.
  /// Flipping them as well would produce a solid of negative volume whose
  /// tessellation faces inward.
  static BrepSolid mirror(BrepSolid s, Vector3 origin, Vector3 normal) {
    if (!_finite(origin) || !_finite(normal) || normal.length2 < 1e-24) {
      throw BrepUnsupported('Mirror needs a finite plane and a nonzero normal');
    }
    final n = normal.normalized();
    return map(
      s,
      (p) => p - n * (2 * (p - origin).dot(n)),
      reverseWinding: true,
      flipFaceNormals: false,
    );
  }

  static BrepSolid scaleAbout(BrepSolid s, Vector3 center, Vector3 scale) {
    if (!_finite(center) ||
        !_finite(scale) ||
        scale.x.abs() < 1e-12 ||
        scale.y.abs() < 1e-12 ||
        scale.z.abs() < 1e-12) {
      throw BrepUnsupported('Scale factors must be finite and nonzero');
    }
    Vector3 stretch(Vector3 v) =>
        Vector3(v.x * scale.x, v.y * scale.y, v.z * scale.z);
    bool equal(double a, double b) =>
        (a - b).abs() <= 1e-10 * (a.abs() + b.abs());
    bool perpendicular(Vector3 a, Vector3 b) =>
        a.dot(b).abs() <= 1e-10 * a.length * b.length;
    final uniform =
        equal(scale.x.abs(), scale.y.abs()) &&
        equal(scale.y.abs(), scale.z.abs());
    for (final face in s.faces) {
      final surface = face.surface;
      if (surface is BrepPlane) continue;
      if (surface is BrepOffsetSurface ||
          surface is BrepRuledSurface ||
          surface is BrepAffineSurface) {
        return affine(
          s,
          BrepAffineTransform.fromFunction((p) => center + stretch(p - center)),
        );
      }
      if (surface is BrepCylinder || surface is BrepCone) {
        final axis = surface is BrepCylinder
            ? surface.axis
            : (surface as BrepCone).axis;
        final x = surface is BrepCylinder
            ? surface.xDir
            : (surface as BrepCone).xDir;
        final y = surface is BrepCylinder
            ? surface.yDir
            : (surface as BrepCone).yDir;
        final ax = stretch(axis), xx = stretch(x), yy = stretch(y);
        if (equal(xx.length, yy.length) &&
            perpendicular(xx, yy) &&
            perpendicular(xx, ax) &&
            perpendicular(yy, ax))
          continue;
      } else if (uniform) {
        continue;
      }
      return affine(
        s,
        BrepAffineTransform.fromFunction((p) => center + stretch(p - center)),
      );
    }
    for (final edge in s.edges) {
      final c = edge.curve;
      if (c is BrepLine) continue;
      if (c is BrepCircle &&
          equal(stretch(c.xDir).length, stretch(c.yDir).length) &&
          perpendicular(stretch(c.xDir), stretch(c.yDir)))
        continue;
      if (c is BrepEllipse && perpendicular(stretch(c.xDir), stretch(c.yDir)))
        continue;
      if (uniform && c is BrepIntersectionCurve) continue;
      return affine(
        s,
        BrepAffineTransform.fromFunction((p) => center + stretch(p - center)),
      );
    }
    return map(
      s,
      (p) => center + stretch(p - center),
      reverseWinding: scale.x * scale.y * scale.z < 0,
      flipFaceNormals: false,
    );
  }

  /// Applies an exact affine representation, composing existing maps rather
  /// than nesting them. All face and edge identities survive by index.
  static BrepSolid affine(BrepSolid solid, BrepAffineTransform transform) {
    final unwrapped = unwrapAffine(solid);
    final base = unwrapped?.solid ?? solid;
    final m = unwrapped == null
        ? transform
        : unwrapped.transform.then(transform);
    return BrepSolid(
      vertices: [for (final p in base.vertices) m.point(p)],
      edges: [
        for (final e in base.edges)
          BrepEdge(
            v0: e.v0,
            v1: e.v1,
            curve: BrepAffineCurve(e.curve, m),
            origin: e.origin,
          ),
      ],
      faces: [
        for (final f in base.faces)
          BrepFace(
            surface: BrepAffineSurface(f.surface, m),
            loops: _loops(f.loops, m.determinant < 0),
            reversed: f.reversed,
            origin: f.origin,
          ),
      ],
    );
  }

  /// Materialize affine spline control points for intersections between bodies
  /// with different placements. Planar faces and all shared edges stay exact.
  static BrepSolid materializeSplineBody(BrepSolid solid) {
    final source = unwrapAffine(solid);
    if (source == null ||
        !source.solid.faces.any((f) => f.surface is BrepSplineSurface) ||
        source.solid.faces.any(
          (f) => f.surface is! BrepSplineSurface && f.surface is! BrepPlane,
        ))
      return solid;
    final m = source.transform;
    return BrepSolid(
      vertices: solid.vertices,
      edges: [
        for (final e in source.solid.edges)
          BrepEdge(
            v0: e.v0,
            v1: e.v1,
            curve: _mapCurve(e.curve, m.point),
            origin: e.origin,
          ),
      ],
      faces: [
        for (var i = 0; i < source.solid.faces.length; i++)
          BrepFace(
            surface: _mapSurface(source.solid.faces[i].surface, m.point),
            loops: solid.faces[i].loops,
            reversed:
                solid.faces[i].reversed !=
                (m.determinant < 0 &&
                    source.solid.faces[i].surface is BrepSplineSurface),
            origin: solid.faces[i].origin,
          ),
      ],
    );
  }

  static List<BrepLoop> _loops(List<BrepLoop> loops, bool reverse) => [
    for (final loop in loops)
      BrepLoop(
        reverse
            ? [
                for (final c in loop.coedges.reversed)
                  BrepCoedge(c.edge, reversed: !c.reversed),
              ]
            : loop.coedges,
      ),
  ];

  /// Returns the original coordinate system when an entire solid shares one
  /// affine map. Operations that commute with affine maps can work there.
  static ({BrepSolid solid, BrepAffineTransform transform})? unwrapAffine(
    BrepSolid s,
  ) {
    if (s.faces.isEmpty || s.faces.first.surface is! BrepAffineSurface)
      return null;
    final m = (s.faces.first.surface as BrepAffineSurface).transform;
    if (s.faces.any(
          (f) =>
              f.surface is! BrepAffineSurface ||
              !(f.surface as BrepAffineSurface).transform.sameAs(m),
        ) ||
        s.edges.any(
          (e) =>
              e.curve is! BrepAffineCurve ||
              !(e.curve as BrepAffineCurve).transform.sameAs(m),
        ))
      return null;
    return (
      transform: m,
      solid: BrepSolid(
        vertices: [for (final p in s.vertices) m.inversePoint(p)],
        edges: [
          for (final e in s.edges)
            BrepEdge(
              v0: e.v0,
              v1: e.v1,
              curve: (e.curve as BrepAffineCurve).base,
              origin: e.origin,
            ),
        ],
        faces: [
          for (final f in s.faces)
            BrepFace(
              surface: (f.surface as BrepAffineSurface).base,
              loops: _loops(f.loops, m.determinant < 0),
              reversed: f.reversed,
              origin: f.origin,
            ),
        ],
      ),
    );
  }

  static bool _finite(Vector3 p) =>
      p.x.isFinite && p.y.isFinite && p.z.isFinite;

  static Vector3 _curveNormal(
    Vector3 origin,
    Vector3 x,
    Vector3 y,
    Vector3 Function(Vector3) f,
  ) {
    final o = f(origin);
    return (f(origin + x) - o).cross(f(origin + y) - o).normalized();
  }

  static BrepCurve _mapCurve(BrepCurve c, Vector3 Function(Vector3) f) {
    if (c is BrepHelix)
      return BrepAffineCurve(c, BrepAffineTransform.fromFunction(f));
    if (c is BrepSurfaceIntersectionCurve)
      return c.transformed(BrepAffineTransform.fromFunction(f));
    if (c is BrepRuledQuadricCurve)
      return BrepAffineCurve(c, BrepAffineTransform.fromFunction(f));
    if (c is BrepPlaneSectionCurve)
      return BrepAffineCurve(c, BrepAffineTransform.fromFunction(f));
    if (c is BrepSplineCurve)
      return c.transformed(BrepAffineTransform.fromFunction(f));
    if (c is BrepAffineCurve)
      return BrepAffineCurve(
        c.base,
        c.transform.then(BrepAffineTransform.fromFunction(f)),
      );
    if (c is BrepIntersectionCurve) {
      var carrier = _mapSurface(c.carrier, f) as BrepCylinder;
      var branch = c.branch;
      final mappedY =
          f(c.carrier.origin + c.carrier.yDir) - f(c.carrier.origin);
      if (mappedY.dot(carrier.yDir) < 0) {
        carrier = BrepCylinder(
          origin: carrier.origin,
          axis: -carrier.axis,
          radius: carrier.radius,
          xDir: carrier.xDir,
        );
        branch = -branch;
      }
      return BrepIntersectionCurve(
        carrier: carrier,
        target: _mapSurface(c.target, f),
        branch: branch,
        t0: c.t0,
        t1: c.t1,
      );
    }
    if (c is BrepLine) return BrepLine(f(c.a), f(c.b));
    if (c is BrepCircle) {
      return BrepCircle(
        center: f(c.center),
        normal: _curveNormal(c.center, c.xDir, c.yDir, f),
        radius: (f(c.center + c.xDir * c.radius) - f(c.center)).length,
        xDir: (f(c.center + c.xDir) - f(c.center)).normalized(),
        t0: c.t0,
        t1: c.t1,
      );
    }
    if (c is BrepEllipse) {
      return BrepEllipse(
        center: f(c.center),
        normal: _curveNormal(c.center, c.xDir, c.yDir, f),
        xDir: (f(c.center + c.xDir) - f(c.center)).normalized(),
        rMajor: (f(c.center + c.xDir * c.rMajor) - f(c.center)).length,
        rMinor: (f(c.center + c.yDir * c.rMinor) - f(c.center)).length,
        t0: c.t0,
        t1: c.t1,
      );
    }
    throw BrepUnsupported('Cannot transform ${c.runtimeType}');
  }

  static BrepSurface _mapSurface(BrepSurface s, Vector3 Function(Vector3) f) {
    if (s is BrepOffsetSurface)
      return BrepAffineSurface(s, BrepAffineTransform.fromFunction(f));
    if (s is BrepSplineSurface)
      return s.transformed(BrepAffineTransform.fromFunction(f));
    if (s is BrepAffineSurface)
      return BrepAffineSurface(
        s.base,
        s.transform.then(BrepAffineTransform.fromFunction(f)),
      );
    if (s is BrepRuledSurface) {
      return BrepRuledSurface(_mapCurve(s.lower, f), _mapCurve(s.upper, f));
    }
    if (s is BrepPlane) {
      final o = f(Vector3.zero());
      final handedness = (f(Vector3(1, 0, 0)) - o).dot(
        (f(Vector3(0, 1, 0)) - o).cross(f(Vector3(0, 0, 1)) - o),
      );
      final normal =
          _curveNormal(s.origin, s.xDir, s.yDir, f) *
          (handedness < 0 ? -1.0 : 1.0);
      return BrepPlane(
        origin: f(s.origin),
        normal: normal,
        xDir: (f(s.origin + s.xDir) - f(s.origin)).normalized(),
      );
    }
    if (s is BrepCylinder) {
      return BrepCylinder(
        origin: f(s.origin),
        axis: (f(s.origin + s.axis) - f(s.origin)).normalized(),
        radius: (f(s.origin + s.xDir * s.radius) - f(s.origin)).length,
        xDir: (f(s.origin + s.xDir) - f(s.origin)).normalized(),
      );
    }
    if (s is BrepSphere) {
      return BrepSphere(
        center: f(s.center),
        radius: (f(s.center + Vector3(s.radius, 0, 0)) - f(s.center)).length,
      );
    }
    if (s is BrepTorus) {
      return BrepTorus(
        center: f(s.center),
        axis: (f(s.center + s.axis) - f(s.center)).normalized(),
        majorRadius:
            (f(s.center + s.xDir * s.majorRadius) - f(s.center)).length,
        minorRadius:
            s.minorRadius * (f(s.center + s.xDir) - f(s.center)).length,
        xDir: (f(s.center + s.xDir) - f(s.center)).normalized(),
      );
    }
    if (s is BrepCone) {
      return BrepCone(
        origin: f(s.origin),
        axis: (f(s.origin + s.axis) - f(s.origin)).normalized(),
        radius: (f(s.origin + s.xDir * s.radius) - f(s.origin)).length,
        topRadius:
            (f(s.origin + s.axis * s.height + s.xDir * s.topRadius) -
                    f(s.origin + s.axis * s.height))
                .length,
        height: (f(s.origin + s.axis * s.height) - f(s.origin)).length,
        xDir: (f(s.origin + s.xDir) - f(s.origin)).normalized(),
      );
    }
    throw BrepUnsupported('Cannot transform ${s.runtimeType}');
  }
}
