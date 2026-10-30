part of 'fillet.dart';

/// Remove internal planar patch boundaries before offsetting a face. Hollow
/// rims are often several strips; a blend must see their union, including its
/// hole, rather than treating each strip's diagonal as a material boundary.
/// Vertex/edge indices stay stable for the caller's selection.
BrepSolid _blendPlanarRegions(BrepSolid s) {
  final adjacency = <int, List<int>>{};
  for (var fi = 0; fi < s.faces.length; fi++) {
    for (final co in s.faces[fi].loops.expand((l) => l.coedges)) {
      (adjacency[co.edge] ??= []).add(fi);
    }
  }
  final visited = <int>{}, faces = <BrepFace>[];
  for (var seed = 0; seed < s.faces.length; seed++) {
    if (!visited.add(seed)) continue;
    final original = s.faces[seed], plane = original.surface;
    if (plane is! BrepPlane) {
      faces.add(original);
      continue;
    }
    final group = <int>{seed}, queue = <int>[seed];
    for (var q = 0; q < queue.length; q++) {
      for (final co in s.faces[queue[q]].loops.expand((l) => l.coedges)) {
        for (final next in adjacency[co.edge]!) {
          if (group.contains(next)) continue;
          final f = s.faces[next], p = f.surface;
          if (p is! BrepPlane ||
              original.outwardNormal(0, 0).dot(f.outwardNormal(0, 0)) <
                  1 - 1e-10 ||
              plane.height(p.origin).abs() > 1e-7) {
            continue;
          }
          group.add(next);
          visited.add(next);
          queue.add(next);
        }
      }
    }
    if (group.length == 1) {
      faces.add(original);
      continue;
    }
    final boundary = <BrepCoedge>[];
    for (final fi in group) {
      for (final co in s.faces[fi].loops.expand((l) => l.coedges)) {
        if (adjacency[co.edge]!.where(group.contains).length == 1) {
          boundary.add(co);
        }
      }
    }
    int start(BrepCoedge c) =>
        c.reversed ? s.edges[c.edge].v1 : s.edges[c.edge].v0;
    int end(BrepCoedge c) =>
        c.reversed ? s.edges[c.edge].v0 : s.edges[c.edge].v1;
    final loops = <BrepLoop>[];
    var valid = true;
    while (boundary.isNotEmpty && valid) {
      final chain = <BrepCoedge>[boundary.removeLast()];
      while (end(chain.last) != start(chain.first)) {
        final matches = boundary
            .where((c) => start(c) == end(chain.last))
            .toList();
        if (matches.length != 1) {
          valid = false;
          break;
        }
        chain.add(matches.single);
        boundary.remove(matches.single);
      }
      loops.add(BrepLoop(chain));
    }
    if (!valid) {
      faces.addAll(group.map((i) => s.faces[i]));
      continue;
    }
    double area(BrepLoop loop) {
      final pts = BrepTessellate.sampleLoop(s, loop, 32);
      var a = Vector3.zero();
      for (var i = 0; i < pts.length; i++) {
        a += (pts[i] - plane.origin).cross(
          pts[(i + 1) % pts.length] - plane.origin,
        );
      }
      return a.dot(plane.normal).abs();
    }

    loops.sort((a, b) => area(b).compareTo(area(a)));
    faces.add(
      BrepFace(
        surface: plane,
        loops: loops,
        reversed: original.reversed,
        origin: original.origin,
      ),
    );
  }
  return BrepSolid(vertices: s.vertices, edges: s.edges, faces: faces);
}

List<int> _creaseEdges(BrepSolid s, List<int> selected) => [
  for (final id in selected.toSet())
    if (id >= 0 && id < s.edges.length && !_tangentBlendEdge(s, id)) id,
];
bool _tangentBlendEdge(BrepSolid s, int id) {
  final adj = BrepFillet.facesForEdge(s, id);
  if (adj.isEmpty) return true; // Removed internal planar seam.
  if (adj.length != 2) return false;
  final c = s.edges[id].curve;
  for (final t in [.1, .5, .9]) {
    final p = c.evaluate(c.t0 + (c.t1 - c.t0) * t);
    final a = s.faces[adj[0]], b = s.faces[adj[1]];
    final ua = a.surface.parameters(p), ub = b.surface.parameters(p);
    if (ua == null ||
        ub == null ||
        a.outwardNormal(ua.$1, ua.$2).dot(b.outwardNormal(ub.$1, ub.$2)) <
            1 - 1e-8) {
      return false;
    }
  }
  return true;
}

BrepPlane? _blendEndPlane(
  BrepSolid s,
  int edge,
  int f1,
  int f2,
  int vertex,
  Vector3 axis,
) {
  BrepPlane? result;
  for (var fi = 0; fi < s.faces.length; fi++) {
    if (fi == f1 || fi == f2) continue;
    final f = s.faces[fi], p = f.surface;
    if (p is! BrepPlane || p.normal.dot(axis).abs() < 1e-8) continue;
    if (!f.loops
        .expand((l) => l.coedges)
        .any(
          (c) => s.edges[c.edge].v0 == vertex || s.edges[c.edge].v1 == vertex,
        )) {
      continue;
    }
    if (result != null && result.normal.dot(p.normal).abs() < 1 - 1e-8) {
      throw BrepUnsupported('Blend endpoint meets multiple end planes');
    }
    result = p;
  }
  return result;
}

Vector3 _clipBlendContact(Vector3 p, Vector3 axis, BrepPlane? end) =>
    end == null ? p : p - axis * (end.height(p) / end.normal.dot(axis));

BrepCurve _blendEndCurve(
  BrepCylinder cylinder,
  BrepPlane? end,
  Vector3 a,
  Vector3 b,
  double ua,
  double ub,
) {
  if (end == null || end.normal.dot(cylinder.axis).abs() > 1 - 1e-10) {
    return BrepCircle(
      center:
          cylinder.origin +
          cylinder.axis * ((a - cylinder.origin).dot(cylinder.axis)),
      normal: cylinder.axis,
      radius: cylinder.radius,
      xDir: cylinder.xDir,
      t0: ua,
      t1: ub,
    );
  }
  final ellipse = BrepIntersect.surfaces(
    end,
    cylinder,
  ).whereType<BrepEllipse>().single;
  final t0 = ellipse.closestParam(a);
  var t1 = ellipse.closestParam(b);
  while (t1 < t0) {
    t1 += 2 * math.pi;
  }
  final desired = _clipBlendContact(
    cylinder.evaluate((ua + ub) * .5, 0),
    cylinder.axis,
    end,
  );
  if (ellipse.evaluate((t0 + t1) * .5).distanceTo(desired) >
      ellipse.evaluate((t0 + t1 - 2 * math.pi) * .5).distanceTo(desired)) {
    t1 -= 2 * math.pi;
  }
  return BrepEllipse(
    center: ellipse.center,
    normal: ellipse.normal,
    xDir: ellipse.xDir,
    rMajor: ellipse.rMajor,
    rMinor: ellipse.rMinor,
    t0: t0,
    t1: t1,
  );
}

/// A convex body can be chamfered by its exact setback plane, including any
/// existing rounds that terminate at that plane. Require supporting tangents
/// on every face so a local chamfer never cuts a second arm of a concave body.
bool _isConvexBlendSolid(BrepSolid s) {
  final points = [
    for (final e in s.edges)
      for (final t in [0.0, .25, .5, .75, 1.0])
        e.curve.evaluate(e.curve.t0 + (e.curve.t1 - e.curve.t0) * t),
  ];
  for (final f in s.faces) {
    if (f.surface is! BrepPlane &&
        f.surface is! BrepCylinder &&
        f.surface is! BrepSphere) {
      return false;
    }
    for (final co in f.loops.expand((l) => l.coedges)) {
      final c = s.edges[co.edge].curve;
      for (final t in [.1, .5, .9]) {
        final p = c.evaluate(c.t0 + (c.t1 - c.t0) * t),
            uv = f.surface.parameters(p);
        if (uv == null) return false;
        final n = f.outwardNormal(uv.$1, uv.$2);
        if (points.any((v) => (v - p).dot(n) > 1e-6)) return false;
      }
    }
  }
  return true;
}

/// Rebuild a radial cylinder sector directly in native analytic geometry.
/// Its radial side planes trim the torus meridians exactly, without extending
/// an open arc into a full circle or running a surface-intersection search.
BrepSolid? _blendClippedCylinder(
  BrepSolid s,
  int edgeId,
  double size, {
  required bool round,
  ChamferOptions options = const ChamferOptions(),
}) {
  final rim = s.edges[edgeId].curve as BrepCircle;
  final cylinders = s.faces.where((f) => f.surface is BrepCylinder).toList();
  if (cylinders.length != 1 || cylinders.single.reversed) return null;
  final cyl = cylinders.single.surface as BrepCylinder;
  final capId = BrepFillet.facesForEdge(
    s,
    edgeId,
  ).firstWhere((id) => s.faces[id].surface is BrepPlane);
  final axis = s.faces[capId].outwardNormal(0, 0);
  if (axis.dot(cyl.axis).abs() < 1 - 1e-8) return null;
  var height = 0.0;
  var sideCount = 0;
  for (final f in s.faces) {
    if (identical(f, cylinders.single)) continue;
    final p = f.surface;
    if (p is! BrepPlane) return null;
    final alignment = p.normal.dot(axis).abs();
    if (alignment < 1e-8) {
      // Offset chord planes need a different endpoint intersection.
      if (p.height(rim.center).abs() > 1e-7) return null;
      sideCount++;
    } else if (alignment < 1 - 1e-8) {
      return null;
    }
  }
  if (sideCount == 0) return null;
  for (final v in s.vertices) {
    final p = v - rim.center, z = p.dot(axis);
    if ((p - axis * z).length > cyl.radius + 1e-7 || z > 1e-7) return null;
    height = math.max(height, -z);
  }
  final depth = round ? size : options.secondDistance(size, math.pi / 2);
  if (size > cyl.radius + 1e-7 || depth >= height - 1e-7) {
    throw BrepUnsupported('Blend consumes the cylinder sector');
  }
  final radius = math.max(0.0, cyl.radius - size);
  final span = (rim.t1 - rim.t0).abs();
  if (span > math.pi + 1e-7) return null;
  final forward = (rim.t1 - rim.t0) * rim.normal.dot(axis) > 0;
  final x = (rim.evaluate(forward ? rim.t0 : rim.t1) - rim.center).normalized();
  final y = axis.cross(x).normalized();
  final origin = rim.center - axis * height;
  Vector3 radial(double u) => x * math.cos(u) + y * math.sin(u);
  Vector3 point(double r, double z, double u) =>
      origin + axis * z + radial(u) * r;
  final profile = <(double, double)>[
    (0, 0),
    (cyl.radius, 0),
    (cyl.radius, height - depth),
    if (radius > 1e-8) (radius, height),
    (0, height),
  ];
  final b = BrepBuilder(), rings = <int?>[], starts = <int>[], ends = <int>[];
  for (final p in profile) {
    final a = b.vertex(point(p.$1, p.$2, 0));
    final z = p.$1 < 1e-8 ? a : b.vertex(point(p.$1, p.$2, span));
    starts.add(a);
    ends.add(z);
    rings.add(
      p.$1 < 1e-8
          ? null
          : b.edge(
              a,
              z,
              BrepCircle(
                center: origin + axis * p.$2,
                normal: axis,
                xDir: x,
                radius: p.$1,
                t0: 0,
                t1: span,
              ),
            ),
    );
  }
  final left = <int>[], right = <int>[];
  for (var i = 0; i < profile.length; i++) {
    final j = (i + 1) % profile.length;
    BrepCurve meridian(double u) {
      final a = point(profile[i].$1, profile[i].$2, u),
          z = point(profile[j].$1, profile[j].$2, u);
      if (i == 2 && round) {
        return BrepCircle(
          center: point(radius, height - depth, u),
          normal: radial(u).cross(axis),
          xDir: radial(u),
          radius: size,
          t0: 0,
          t1: math.pi / 2,
        );
      }
      return BrepLine(a, z);
    }

    final le = b.edge(starts[i], starts[j], meridian(0));
    final re = (rings[i] == null && rings[j] == null)
        ? le
        : b.edge(ends[i], ends[j], meridian(span));
    left.add(le);
    right.add(re);
    if (rings[i] == null && rings[j] == null) continue;
    BrepSurface surface;
    if (i == 0) {
      surface = BrepPlane(origin: origin, normal: -axis);
    } else if (i == 1) {
      surface = BrepCylinder(
        origin: origin,
        axis: axis,
        radius: cyl.radius,
        xDir: x,
      );
    } else if (i == 2) {
      surface = round
          ? (radius < 1e-8
                ? BrepSphere(
                    center: origin + axis * (height - depth),
                    radius: size,
                  )
                : BrepTorus(
                    center: origin + axis * (height - depth),
                    axis: axis,
                    majorRadius: radius,
                    minorRadius: size,
                    xDir: x,
                  ))
          : BrepCone(
              origin: origin + axis * (height - depth),
              axis: axis,
              radius: cyl.radius,
              topRadius: radius,
              height: depth,
              xDir: x,
            );
    } else {
      surface = BrepPlane(origin: rim.center, normal: axis);
    }
    b.face(surface, [
      BrepLoop([
        if (rings[i] != null) BrepCoedge(rings[i]!),
        BrepCoedge(re),
        if (rings[j] != null) BrepCoedge(rings[j]!, reversed: true),
        BrepCoedge(le, reversed: true),
      ]),
    ]);
  }
  b.face(BrepPlane(origin: origin, normal: -y), [
    BrepLoop([for (final e in left) BrepCoedge(e)]),
  ]);
  b.face(BrepPlane(origin: origin, normal: axis.cross(radial(span))), [
    BrepLoop([for (final e in right.reversed) BrepCoedge(e, reversed: true)]),
  ]);
  return b.build();
}

/// A setback can exactly consume a planar strip (for example a shell rim).
/// Split shared straight boundaries and cancel the zero-area strip, retaining
/// the adjoining curved face and inner wall rather than inventing a cap.
BrepSolid _resolveBlendLineContacts(BrepSolid s) {
  if (!s.faces.any((f) => f.surface is BrepPlane && f.loops.length > 1)) {
    return s;
  }
  final used = s.faces
      .expand((f) => f.loops)
      .expand((l) => l.coedges)
      .map((c) => c.edge)
      .toSet();
  final vertices = used
      .expand((id) => [s.edges[id].v0, s.edges[id].v1])
      .toSet();
  final b = BrepBuilder(), vm = <int, int>{};
  int vertex(int id) => vm.putIfAbsent(id, () {
    final p = s.vertices[id];
    for (var i = 0; i < b.vertices.length; i++) {
      if (b.vertices[i].distanceTo(p) < 1e-8) return i;
    }
    return b.vertex(p);
  });
  final chains = <int, List<BrepCoedge>>{}, lines = <(int, int), int>{};
  var shared = false;
  for (final id in used) {
    final e = s.edges[id], c = e.curve;
    if (c is! BrepLine) {
      chains[id] = [
        BrepCoedge(b.edge(vertex(e.v0), vertex(e.v1), c, origin: e.origin)),
      ];
      continue;
    }
    final a = s.vertices[e.v0], d = s.vertices[e.v1] - a;
    if (d.length2 < 1e-18) {
      chains[id] = [];
      continue;
    }
    final breaks = <(double, int)>[(0, e.v0), (1, e.v1)];
    for (final v in vertices) {
      if (v == e.v0 || v == e.v1) continue;
      final t = (s.vertices[v] - a).dot(d) / d.length2;
      if (t > 1e-8 &&
          t < 1 - 1e-8 &&
          (a + d * t).distanceTo(s.vertices[v]) < 1e-8) {
        breaks.add((t, v));
      }
    }
    breaks.sort((a, b) => a.$1.compareTo(b.$1));
    final chain = <BrepCoedge>[];
    for (var i = 1; i < breaks.length; i++) {
      final a = vertex(breaks[i - 1].$2), z = vertex(breaks[i].$2);
      if (a == z) continue;
      final key = a < z ? (a, z) : (z, a);
      final old = lines[key];
      if (old != null) shared = true;
      final edge =
          old ??
          b.edge(
            a,
            z,
            BrepLine(b.vertices[a], b.vertices[z]),
            origin: e.origin,
          );
      lines[key] = edge;
      chain.add(BrepCoedge(edge, reversed: b.edges[edge].v0 != a));
    }
    chains[id] = chain;
  }
  if (!shared) return s;
  var cancelledAny = false;
  for (final f in s.faces) {
    var loops = [
      for (final l in f.loops)
        BrepLoop([
          for (final co in l.coedges)
            for (final c
                in co.reversed ? chains[co.edge]!.reversed : chains[co.edge]!)
              BrepCoedge(c.edge, reversed: c.reversed != co.reversed),
        ]),
    ];
    if (f.surface is BrepPlane) {
      final all = loops.expand((l) => l.coedges).toList(),
          grouped = <int, List<BrepCoedge>>{};
      for (final c in all) {
        (grouped[c.edge] ??= []).add(c);
      }
      final cancelled = grouped.entries
          .where(
            (e) =>
                e.value.length == 2 &&
                e.value[0].reversed != e.value[1].reversed,
          )
          .map((e) => e.key)
          .toSet();
      if (cancelled.isNotEmpty) {
        cancelledAny = true;
        final remain = all.where((c) => !cancelled.contains(c.edge)).toList();
        loops = [];
        int start(BrepCoedge c) =>
            c.reversed ? b.edges[c.edge].v1 : b.edges[c.edge].v0;
        int end(BrepCoedge c) =>
            c.reversed ? b.edges[c.edge].v0 : b.edges[c.edge].v1;
        while (remain.isNotEmpty) {
          final wire = <BrepCoedge>[remain.removeLast()];
          while (end(wire.last) != start(wire.first)) {
            final next = remain
                .where((c) => start(c) == end(wire.last))
                .toList();
            if (next.length != 1) return s;
            wire.add(next.single);
            remain.remove(next.single);
          }
          if (wire.length >= 3) loops.add(BrepLoop(wire));
        }
        final temporary = BrepSolid(
          vertices: b.vertices,
          edges: b.edges,
          faces: [],
        );
        double area(BrepLoop l) {
          final pts = BrepTessellate.sampleLoop(temporary, l, 32);
          var a = Vector3.zero();
          for (var i = 0; i < pts.length; i++) {
            a += pts[i].cross(pts[(i + 1) % pts.length]);
          }
          return a.length2;
        }

        loops.sort((a, b) => area(b).compareTo(area(a)));
      }
    }
    if (loops.isNotEmpty) {
      b.face(f.surface, loops, reversed: f.reversed, origin: f.origin);
    }
  }
  if (!cancelledAny) return s;
  final out = b.build();
  return BrepHeal.topologyIsClosed(out) ? out : s;
}

Map<int, List<int>> _blendAdjacency(BrepSolid s) {
  final adjacency = <int, List<int>>{};
  for (var fi = 0; fi < s.faces.length; fi++) {
    for (final id
        in s.faces[fi].loops
            .expand((l) => l.coedges)
            .map((c) => c.edge)
            .toSet()) {
      (adjacency[id] ??= []).add(fi);
    }
  }
  return adjacency;
}
