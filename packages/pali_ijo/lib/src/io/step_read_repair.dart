part of 'step_read_analytic.dart';

extension _StepSourceRepair on StepReader {
  /// A corrupt cylinder can contain an isolated complete rim and a wire that
  /// retraces itself. There is no finite cylindrical wall defined by those
  /// bounds. Opt-in repair caps the isolated opening and sews the retraced
  /// wire. It changes the source shape and must always be reported to the user.
  BrepSolid _repairCollapsedCylinders(BrepSolid source, int bodyId) {
    if (!BrepHeal.topologyIsClosed(source)) return source;
    var solid = source;
    for (var fi = 0; fi < solid.faces.length; fi++) {
      final face = solid.faces[fi], surface = face.surface;
      if (surface is! BrepCylinder || face.loops.length != 2) continue;
      final rims = face.loops.where((l) => l.coedges.length == 1).toList();
      final wires = face.loops.where((l) => l.coedges.length == 2).toList();
      if (rims.length != 1 || wires.length != 1) continue;
      final rimEdge = rims.single.coedges.single.edge;
      final rim = solid.edges[rimEdge].curve;
      if (rim is! BrepCircle ||
          !rim.isClosed ||
          rim.normal.cross(surface.axis).length > 1e-8 ||
          (rim.radius - surface.radius).abs() > 1e-7 ||
          (surface.closestPoint(rim.evaluate(rim.t0)) - rim.evaluate(rim.t0))
                  .length >
              1e-7)
        continue;
      final aUse = wires.single.coedges[0], bUse = wires.single.coedges[1];
      if (aUse.edge == bUse.edge) continue; // A legitimate periodic seam.
      final a = solid.edges[aUse.edge], b = solid.edges[bUse.edge];
      final same = a.v0 == b.v0 && a.v1 == b.v1;
      final opposite = a.v0 == b.v1 && a.v1 == b.v0;
      if (!same && !opposite || a.v0 == a.v1) continue;
      if ((aUse.reversed != bUse.reversed) != same) continue;
      if (!List.generate(17, (i) => i / 16).every((t) {
        final p = a.curve.evaluate(a.curve.t0 + (a.curve.t1 - a.curve.t0) * t);
        final q = b.curve.evaluate(
          b.curve.t0 + (b.curve.t1 - b.curve.t0) * (same ? t : 1 - t),
        );
        return (p - q).length < 1e-7;
      }))
        continue;
      final rimNeighbours = <(int, int)>[];
      for (var other = 0; other < solid.faces.length; other++) {
        if (other == fi) continue;
        final f = solid.faces[other];
        for (var li = 0; li < f.loops.length; li++) {
          if (f.loops[li].coedges.any((c) => c.edge == rimEdge)) {
            rimNeighbours.add((other, li));
          }
        }
      }
      if (rimNeighbours.length != 1) continue;
      final (capIndex, capLoop) = rimNeighbours.single;
      final cap = solid.faces[capIndex];
      // Only fill a complete inner circle of its existing planar neighbour.
      // Never invent an outer cap, or erase an ordinary noncollapsed tube.
      if (cap.surface is! BrepPlane ||
          capLoop == 0 ||
          cap.loops[capLoop].coedges.length != 1 ||
          (cap.surface as BrepPlane).normal.cross(rim.normal).length > 1e-8)
        continue;
      final faces = <BrepFace>[];
      for (var other = 0; other < solid.faces.length; other++) {
        if (other == fi) continue;
        final f = solid.faces[other];
        faces.add(
          BrepFace(
            surface: f.surface,
            reversed: f.reversed,
            origin: f.origin,
            loops: [
              for (var li = 0; li < f.loops.length; li++)
                if (other != capIndex || li != capLoop)
                  BrepLoop([
                    for (final co in f.loops[li].coedges)
                      co.edge == bUse.edge
                          ? BrepCoedge(
                              aUse.edge,
                              reversed: co.reversed != !same,
                            )
                          : co,
                  ]),
            ],
          ),
        );
      }
      final candidate = BrepSolid(
        vertices: solid.vertices,
        edges: solid.edges,
        faces: faces,
      );
      if (!BrepHeal.topologyIsClosed(candidate)) continue;
      solid = candidate;
      fi--;
      _notes.add(
        'Repaired invalid STEP body #$bodyId: removed a collapsed cylindrical face, '
        'sewed its coincident edges, and capped its malformed ${(2 * rim.radius).toStringAsPrecision(8)} mm circular opening. '
        'This repair changes the source shape; the remaining native CAD surfaces are retained.',
      );
    }
    if (identical(solid, source)) return source;
    final vertexMap = <int, int>{}, edgeMap = <int, int>{};
    final vertices = <Vector3>[], edges = <BrepEdge>[];
    int vertex(int old) => vertexMap.putIfAbsent(old, () {
      vertices.add(solid.vertices[old]);
      return vertices.length - 1;
    });
    int edge(int old) => edgeMap.putIfAbsent(old, () {
      final e = solid.edges[old];
      edges.add(
        BrepEdge(
          v0: vertex(e.v0),
          v1: vertex(e.v1),
          curve: e.curve,
          origin: e.origin,
        ),
      );
      return edges.length - 1;
    });
    final faces = [
      for (final f in solid.faces)
        BrepFace(
          surface: f.surface,
          reversed: f.reversed,
          origin: f.origin,
          loops: [
            for (final l in f.loops)
              BrepLoop([
                for (final c in l.coedges)
                  BrepCoedge(edge(c.edge), reversed: c.reversed),
              ]),
          ],
        ),
    ];
    return BrepSolid(vertices: vertices, edges: edges, faces: faces);
  }
}
