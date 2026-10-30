part of 'xform.dart';

bool _equalScale(double a, double b) =>
    (a - b).abs() <= 1e-10 * (a.abs() + b.abs());
bool _orthogonal(Vector3 a, Vector3 b) =>
    a.dot(b).abs() <= 1e-10 * a.length * b.length;
bool _similarity(BrepAffineTransform m) =>
    _equalScale(m.x.length, m.y.length) &&
    _equalScale(m.x.length, m.z.length) &&
    _orthogonal(m.x, m.y) &&
    _orthogonal(m.x, m.z) &&
    _orthogonal(m.y, m.z);

BrepCurve _editCurve(BrepCurve curve) {
  if (curve is! BrepAffineCurve) return curve;
  var c = curve.base, m = curve.transform;
  while (c is BrepAffineCurve) {
    m = c.transform.then(m);
    c = c.base;
  }
  if (c is BrepLine) return BrepLine(m.point(c.a), m.point(c.b));
  if (c is BrepSplineCurve) return c.transformed(m);
  if (c is BrepCircle || c is BrepEllipse) {
    final center = c is BrepCircle ? c.center : (c as BrepEllipse).center;
    final x = c is BrepCircle ? c.xDir : (c as BrepEllipse).xDir;
    final y = c is BrepCircle ? c.yDir : (c as BrepEllipse).yDir;
    final xx = m.vector(x), yy = m.vector(y);
    if (_orthogonal(xx, yy)) {
      final rx =
          (c is BrepCircle ? c.radius : (c as BrepEllipse).rMajor) * xx.length;
      final ry =
          (c is BrepCircle ? c.radius : (c as BrepEllipse).rMinor) * yy.length;
      if (_equalScale(rx, ry)) {
        return BrepCircle(
          center: m.point(center),
          normal: xx.cross(yy),
          xDir: xx.normalized(),
          radius: rx,
          t0: c.t0,
          t1: c.t1,
        );
      }
      return BrepEllipse(
        center: m.point(center),
        normal: xx.cross(yy),
        xDir: xx.normalized(),
        rMajor: rx,
        rMinor: ry,
        t0: c.t0,
        t1: c.t1,
      );
    }
  }
  return BrepAffineCurve(c, m);
}

BrepCurve _editBoundaryCurve(BrepCurve curve) {
  // Existing spline trims are already supported by the editors. Retain their
  // representation (and shared display sampling) while exposing the primitive
  // lines and conics that those editors must recognize directly.
  var base = curve;
  while (base is BrepAffineCurve) {
    base = base.base;
  }
  return base is BrepSplineCurve ? curve : _editCurve(curve);
}

bool _editableRail(BrepCurve curve) {
  while (curve is BrepAffineCurve) {
    curve = curve.base;
  }
  return curve is BrepLine ||
      curve is BrepCircle ||
      curve is BrepEllipse ||
      curve is BrepSplineCurve;
}

BrepSolid _materializeForEditing(BrepSolid solid, bool preserveAffineSweeps) {
  if (!solid.faces.any((f) => f.surface is BrepAffineSurface) &&
      !solid.edges.any((e) => e.curve is BrepAffineCurve)) {
    return solid;
  }
  final faces = <BrepFace>[];
  final reverseWires = <int>{};
  for (final face in solid.faces) {
    var surface = face.surface;
    var reversed = face.reversed;
    if (surface is BrepAffineSurface) {
      var base = surface.base, m = surface.transform;
      while (base is BrepAffineSurface) {
        m = base.transform.then(m);
        base = base.base;
      }
      if (base is BrepPlane) {
        surface = BrepPlane(
          origin: m.point(base.origin),
          normal: m.normal(base.normal),
          xDir: m.vector(base.xDir),
        );
      } else if (base is BrepCylinder || base is BrepCone) {
        final axis = base is BrepCylinder ? base.axis : (base as BrepCone).axis;
        final x = base is BrepCylinder ? base.xDir : (base as BrepCone).xDir;
        final y = base is BrepCylinder ? base.yDir : (base as BrepCone).yDir;
        final ax = m.vector(axis), xx = m.vector(x), yy = m.vector(y);
        if (_equalScale(xx.length, yy.length) &&
            _orthogonal(xx, yy) &&
            _orthogonal(xx, ax) &&
            _orthogonal(yy, ax)) {
          surface = BrepXform._mapSurface(base, m.point);
        } else if (base is BrepCylinder) {
          final cylinder = base;
          // A stretched cylinder is an exact translated ellipse (possibly
          // oblique), not a circular cylinder. Bound its axial rail using the
          // linear/circular original trims; unknown curved trims stay wrapped.
          final edges = face.loops
              .expand((l) => l.coedges)
              .map((co) => solid.edges[co.edge])
              .toList();
          bool simple(BrepCurve c) {
            while (c is BrepAffineCurve) {
              c = c.base;
            }
            return c is BrepLine ||
                c is BrepCircle && c.normal.cross(axis).length < 1e-8;
          }

          if (edges.isNotEmpty && edges.every((e) => simple(e.curve))) {
            final zs = [
              for (final e in edges)
                for (final i in [e.v0, e.v1])
                  (m.inversePoint(solid.vertices[i]) - cylinder.origin).dot(
                    axis,
                  ),
            ]..sort();
            final inverse = BrepAffineTransform.fromFunction(m.inversePoint);
            final circles = edges
                .map((e) => _editCurve(BrepAffineCurve(e.curve, inverse)))
                .whereType<BrepCircle>()
                .toList();
            double start(BrepCircle c) => cylinder
                .parameters(
                  c.evaluate(
                    c.normal.dot(axis) * (c.t1 - c.t0) >= 0 ? c.t0 : c.t1,
                  ),
                )!
                .$1;
            final first = circles.firstOrNull;
            final u0 = first == null ? 0.0 : start(first);
            final span = first == null ? 0.0 : (first.t1 - first.t0).abs();
            final matching =
                first != null &&
                circles.every(
                  (c) =>
                      ((c.t1 - c.t0).abs() - span).abs() < 1e-8 &&
                      cylinder
                              .evaluate(start(c), 0)
                              .distanceTo(cylinder.evaluate(u0, 0)) <
                          1e-7,
                );
            if (matching && zs.last - zs.first > 1e-9) {
              BrepCurve rail(double z) => _editCurve(
                BrepAffineCurve(
                  BrepCircle(
                    center: cylinder.origin + axis * z,
                    normal: axis,
                    xDir: x,
                    radius: cylinder.radius,
                    t0: u0,
                    t1: u0 + span,
                  ),
                  m,
                ),
              );
              surface = BrepRuledSurface(rail(zs.first), rail(zs.last));
              reversed ^= m.determinant < 0;
              // Analytic cylinders can carry opposite parameter frames on
              // their cap circles. Their special cylindrical trim handling
              // accepts this, but a ruled wall needs a continuous UV wire.
              // Reverse both uses of a closed rim when the two rims would
              // otherwise wind the same way around the periodic chart.
              if (span >= 2 * math.pi - 1e-8 && face.loops.length == 1) {
                final uses = face.loops.single.coedges.where((co) {
                  final e = solid.edges[co.edge];
                  return e.v0 == e.v1 && e.curve is! BrepLine;
                }).toList();
                if (uses.length == 2) {
                  double direction(BrepCoedge co) {
                    final c = _editCurve(solid.edges[co.edge].curve);
                    final p = c.evaluate(c.t0);
                    final u = (surface as BrepRuledSurface).parameters(p).$1;
                    return surface.du(u, 0).dot(c.tangent(c.t0)) *
                        (c.t1 - c.t0).sign *
                        (co.reversed ? -1 : 1);
                  }

                  if (direction(uses[0]) * direction(uses[1]) > 0) {
                    reverseWires.add(uses[1].edge);
                  }
                }
              }
            }
          }
        }
      } else if (!preserveAffineSweeps &&
          base is BrepRuledSurface &&
          _editableRail(base.lower) &&
          _editableRail(base.upper)) {
        surface = BrepRuledSurface(
          _editCurve(BrepAffineCurve(base.lower, m)),
          _editCurve(BrepAffineCurve(base.upper, m)),
        );
        reversed ^= m.determinant < 0;
      } else if (!preserveAffineSweeps && base is BrepSplineSurface) {
        surface = base.transformed(m);
        reversed ^= m.determinant < 0;
      } else if (_similarity(m) && (base is BrepSphere || base is BrepTorus)) {
        surface = BrepXform._mapSurface(base, m.point);
      }
    }
    faces.add(
      BrepFace(
        surface: surface,
        loops: face.loops,
        reversed: reversed,
        origin: face.origin,
      ),
    );
  }
  if (reverseWires.isNotEmpty) {
    for (final face in faces) {
      face.loops = [
        for (final loop in face.loops)
          BrepLoop([
            for (final co in loop.coedges)
              BrepCoedge(
                co.edge,
                reversed: co.reversed ^ reverseWires.contains(co.edge),
              ),
          ]),
      ];
    }
  }
  return BrepSolid(
    vertices: solid.vertices,
    faces: faces,
    edges: [
      for (final e in solid.edges)
        BrepEdge(
          v0: e.v0,
          v1: e.v1,
          curve: _editBoundaryCurve(e.curve),
          origin: e.origin,
        ),
    ],
  );
}
