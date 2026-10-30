/// Reads shared STEP topology, native analytic surfaces and rational B-splines.
/// Unsupported or incomplete bodies reject as a unit; faces and cavities are
/// never silently omitted. Display tessellation is a separate host operation.
library;

import 'dart:convert';
import 'dart:math' as math;
import 'dart:typed_data';

import 'package:vector_math/vector_math_64.dart';

import '../brep/geom.dart';
import '../ijo_log.dart';
import '../brep/solid.dart';
import '../brep/heal.dart';
import '../brep/spline_convert.dart';
import '../brep/xform.dart';
import 'step_syntax.dart';

part 'step_read_assembly.dart';
part 'step_read_repair.dart';

/// Raised when a STEP file is valid but describes geometry this kernel has no
/// representation for.
///
/// Held apart from [StepSyntaxException] because the two call for different
/// responses: a malformed file has invalid syntax, while unsupported surface families
/// may be valid geometry that requires another representation.
class StepBrepUnsupportedException implements Exception {
  StepBrepUnsupportedException(this.message, {this.entities = const {}});

  final String message;

  /// The STEP entity types that could not be mapped, and how many of each.
  /// Useful for a message that says what to change in the exporting system.
  final Map<String, int> entities;

  @override
  String toString() {
    if (entities.isEmpty) return message;
    final list =
        (entities.entries.toList()..sort((a, b) => b.value.compareTo(a.value)))
            .map((e) => e.value == 1 ? e.key : '${e.key} ×${e.value}')
            .join(', ');
    return '$message Found: $list.';
  }
}

/// What a STEP file yielded.
class StepBrepReadResult {
  StepBrepReadResult({
    required this.solids,
    required this.names,
    required this.unitName,
    required this.unitScale,
    required this.notes,
  });

  final List<BrepSolid> solids;

  /// Each solid's name, from its `MANIFOLD_SOLID_BREP`. Same length as
  /// [solids]; entries may be empty when the file named nothing.
  final List<String> names;

  /// The length unit the file declared, e.g. `mm` or `inch`.
  final String unitName;

  /// What one file unit is in millimetres. Already applied to [solids]; kept
  /// so a caller can report what it did.
  final double unitScale;

  /// Nonfatal import information. Missing boundaries and open shells fail.
  final List<String> notes;

  bool get isEmpty => solids.isEmpty;
}

/// Reads solids out of a STEP file.
class StepReader {
  StepReader._(
    this._file,
    this._scale,
    this._unitName,
    this._allowOpenShells,
    this._repairDegenerateFaces,
  );

  final bool _allowOpenShells;
  final bool _repairDegenerateFaces;

  final StepFile _file;
  final double _scale;
  final String _unitName;

  final _notes = <String>[];
  final _unsupported = <String, int>{};
  late final double _lengthTolerance = _readLengthTolerance();

  double _readLengthTolerance() {
    var tolerance = 1e-5;
    for (final id in _file.idsOf('UNCERTAINTY_MEASURE_WITH_UNIT')) {
      final measure = _file.typed(id, 'UNCERTAINTY_MEASURE_WITH_UNIT')!;
      final unit = measure.arg(1);
      final factor = unit is StepRef ? _lengthOf(_file, unit.id)?.$1 : null;
      final value = _number(measure.arg(0));
      if (factor != null && value != null && value.isFinite && value > 0) {
        tolerance = math.max(tolerance, value * factor);
      }
    }
    return tolerance;
  }

  /// Our vertex index for each STEP vertex instance, and for each welded
  /// position when a face gives points rather than vertices.
  final _vertexOfInstance = <int, int>{};
  final _vertexOfPosition = <String, int>{};

  /// Our edge index for each `EDGE_CURVE` instance. Shared, so the two faces
  /// either side of an edge refer to one edge as they do in the file.
  final _edgeOfInstance = <int, int>{};
  final _lineEdges = <(int, int), int>{};

  final _vertices = <Vector3>[];
  final _edges = <BrepEdge>[];

  static StepBrepReadResult read(
    Uint8List bytes, {
    bool allowOpenShells = false,
    bool repairDegenerateFaces = false,
    void Function(int completed, int total)? onProgress,
  }) {
    final file = StepFile.parse(_decode(bytes));
    final (scale, unitName) = _lengthUnit(file);
    return StepReader._(
      file,
      scale,
      unitName,
      allowOpenShells,
      repairDegenerateFaces,
    )._run(onProgress);
  }

  /// STEP is Latin-1 by the standard, and files in the wild are usually ASCII
  /// with the odd UTF-8 name. Trying UTF-8 first and falling back keeps both
  /// readable instead of failing on a byte over 0x7f.
  static String _decode(Uint8List bytes) {
    try {
      return utf8.decode(bytes);
    } on FormatException {
      return latin1.decode(bytes, allowInvalid: true);
    }
  }

  StepBrepReadResult _run(void Function(int, int)? onProgress) {
    final bodies = [
      ..._file.idsOf('MANIFOLD_SOLID_BREP'),
      ..._file.idsOf('FACETED_BREP'),
      ..._file.idsOf('BREP_WITH_VOIDS'),
      if (_allowOpenShells) ..._file.idsOf('SHELL_BASED_SURFACE_MODEL'),
    ].toSet().toList()..sort();

    if (bodies.isEmpty) {
      // Worth distinguishing: a surface model or a wireframe is a valid STEP
      // file that simply has no solid in it, and the fix is on the exporting
      // side.
      final kind = [
        if (_file.idsOf('SHELL_BASED_SURFACE_MODEL').isNotEmpty)
          'a surface model',
        if (_file.idsOf('GEOMETRIC_CURVE_SET').isNotEmpty) 'a wireframe',
        if (_file.idsOf('TESSELLATED_SHAPE_REPRESENTATION').isNotEmpty)
          'a tessellation',
      ];
      throw StepBrepUnsupportedException(
        kind.isEmpty
            ? 'This STEP file contains no solid body.'
            : 'This STEP file contains ${kind.join(' and ')} rather than a '
                  'solid body. Re-export it as a solid.',
      );
    }

    final solids = <BrepSolid>[];
    final names = <String>[];
    onProgress?.call(0, bodies.length);
    for (final id in bodies) {
      final record = _file.anyOf(id, [
        'MANIFOLD_SOLID_BREP',
        'FACETED_BREP',
        'BREP_WITH_VOIDS',
        if (_allowOpenShells) 'SHELL_BASED_SURFACE_MODEL',
      ])!;
      final solid = _solid(id, record);
      solids.add(solid);
      final name = record.arg(0);
      names.add(name is StepText ? name.value : '');
      onProgress?.call(solids.length, bodies.length);
    }

    if (solids.isEmpty) {
      throw StepBrepUnsupportedException(
        'None of the ${bodies.length} bodies in this STEP file could be read, '
        'because their surfaces are outside what pali ijo represents. It '
        'handles analytic and supported rational spline surfaces.',
        entities: Map.unmodifiable(_unsupported),
      );
    }
    if (_unsupported.isNotEmpty) {
      throw StepBrepUnsupportedException(
        'STEP contains unsupported geometry; no partial body was imported.',
        entities: Map.unmodifiable(_unsupported),
      );
    }

    final placed = _placeAssembly(bodies, solids, names);
    return StepBrepReadResult(
      solids: placed.$1,
      names: placed.$2,
      unitName: _unitName,
      unitScale: _scale,
      notes: List.unmodifiable(_notes),
    );
  }

  BrepSolid _solid(int id, StepRecord record) {
    _vertexOfInstance.clear();
    _vertexOfPosition.clear();
    _edgeOfInstance.clear();
    _lineEdges.clear();
    _vertices.clear();
    _edges.clear();

    // BREP_WITH_VOIDS carries an outer shell plus the voids inside it. The
    // voids are read too: a hollow part is one shell set, and dropping the
    // inner shells would fill the cavity in.
    final shells = <(int, bool)>[];
    final outer = record.arg(1);
    if (record.type == 'SHELL_BASED_SURFACE_MODEL' && outer is StepList) {
      for (final shell in outer.items) {
        if (shell is! StepRef)
          throw const FormatException('Invalid STEP surface shell');
        shells.add((shell.id, false));
      }
    } else if (outer is! StepRef) {
      throw StepBrepUnsupportedException(
        'Body #$id has no readable outer shell.',
      );
    } else {
      shells.add((outer.id, false));
    }
    if (record.type == 'BREP_WITH_VOIDS') {
      final voids = record.arg(2);
      if (voids is! StepList || voids.items.isEmpty) {
        throw StepBrepUnsupportedException(
          'Body #$id has no readable void shells.',
        );
      }
      for (final item in voids.items) {
        final shell = _file.follow(item, 'ORIENTED_CLOSED_SHELL');
        final ref = shell?.arg(2);
        final orientation = shell?.arg(3);
        if (ref is! StepRef || orientation is! StepEnum) {
          throw StepBrepUnsupportedException(
            'Body #$id has an unreadable void shell.',
          );
        }
        shells.add((ref.id, orientation.isFalse));
      }
    }

    final faces = <BrepFace>[];
    var attempted = 0;
    for (final (shellId, reverseShell) in shells) {
      final shell = _file.anyOf(shellId, [
        'CLOSED_SHELL',
        'OPEN_SHELL',
        'ORIENTED_CLOSED_SHELL',
      ]);
      if (shell == null) {
        throw StepBrepUnsupportedException(
          'Body #$id references a missing shell #$shellId.',
        );
      }
      final list = shell.arg(1);
      if (list is! StepList || list.items.isEmpty) {
        throw StepBrepUnsupportedException(
          'Body #$id has an unreadable shell #$shellId.',
        );
      }
      for (final item in list.items) {
        attempted++;
        if (item is! StepRef) continue;
        final face = _face(item.id);
        if (face != null)
          faces.add(
            reverseShell
                ? BrepFace(
                    surface: face.surface,
                    loops: [
                      for (final l in face.loops)
                        BrepLoop(_reverseLoop(l.coedges)),
                    ],
                    reversed: !face.reversed,
                  )
                : face,
          );
      }
    }

    if (faces.isEmpty || faces.length < attempted) {
      throw StepBrepUnsupportedException(
        'Body #$id has unreadable faces; refusing an incomplete solid.',
        entities: Map.unmodifiable(_unsupported),
      );
    }

    var solid = _compact(faces);
    if (_repairDegenerateFaces) solid = _repairCollapsedCylinders(solid, id);
    if (!BrepHeal.topologyIsClosed(solid)) {
      if (!_allowOpenShells) {
        throw StepBrepUnsupportedException(
          'Body #$id does not have a closed BREP boundary',
        );
      }
      _notes.add(
        'Body #$id contains open boundaries in the STEP source. All faces were retained; solid operations may be unavailable.',
      );
    }
    return solid;
  }

  /// Builds the solid, keeping only the edges and vertices that faces use.
  ///
  /// Collapsed periodic seams leave their edges behind, and an edge
  /// belonging to no face is not just clutter — the kernel's own operations
  /// walk the edge list to work out adjacency, so a stray edge reads as a
  /// crack in the topology.
  BrepSolid _compact(List<BrepFace> faces) {
    final edgeMap = <int, int>{};
    final vertexMap = <int, int>{};
    final edges = <BrepEdge>[];
    final vertices = <Vector3>[];

    int keepVertex(int old) => vertexMap.putIfAbsent(old, () {
      vertices.add(_vertices[old]);
      return vertices.length - 1;
    });

    for (final face in faces) {
      for (final loop in face.loops) {
        for (var i = 0; i < loop.coedges.length; i++) {
          final coedge = loop.coedges[i];
          final index = edgeMap.putIfAbsent(coedge.edge, () {
            final old = _edges[coedge.edge];
            edges.add(
              BrepEdge(
                v0: keepVertex(old.v0),
                v1: keepVertex(old.v1),
                curve: old.curve,
              ),
            );
            return edges.length - 1;
          });
          loop.coedges[i] = BrepCoedge(index, reversed: coedge.reversed);
        }
      }
    }

    return BrepSolid(vertices: vertices, edges: edges, faces: faces);
  }

  BrepFace? _face(int id) {
    final faceted = _file.typed(id, 'FACE');
    if (faceted != null) {
      final bounds = faceted.arg(1);
      final loops = <BrepLoop>[];
      if (bounds is! StepList) return null;
      for (final item in bounds.items) {
        final bound = item is StepRef
            ? _file.anyOf(item.id, ['FACE_OUTER_BOUND', 'FACE_BOUND'])
            : null;
        final poly = _file.follow(bound?.arg(1), 'POLY_LOOP');
        if (poly == null) return null;
        final sense = bound?.arg(2);
        final loop = _polyLoop(poly, sense is StepEnum && sense.isFalse);
        if (loop == null) return null;
        loops.add(loop);
      }
      if (loops.isEmpty) return null;
      final points = [
        for (final co in loops.first.coedges)
          _vertices[co.reversed ? _edges[co.edge].v1 : _edges[co.edge].v0],
      ];
      var normal = Vector3.zero();
      for (var i = 0; i < points.length; i++) {
        normal += points[i].cross(points[(i + 1) % points.length]);
      }
      if (normal.length < 1e-12) return null;
      final plane = BrepPlane(origin: points.first, normal: normal);
      if (points.any((p) => plane.height(p).abs() > 1e-6)) return null;
      return BrepFace(surface: plane, loops: loops);
    }
    final record = _file.anyOf(id, ['ADVANCED_FACE', 'FACE_SURFACE']);
    if (record == null) {
      _count(_file.typeName(id));
      return null;
    }

    final sameSense = record.arg(3);
    final reversed = sameSense is StepEnum && sameSense.isFalse;

    final bounds = record.arg(1);
    final loops = <BrepLoop>[];
    if (bounds is! StepList || bounds.items.isEmpty) return null;
    for (final item in bounds.items) {
      final loop = _bound(item);
      if (loop == null) return null;
      if (!_isSeam(loop)) loops.add(loop);
    }

    final surface = _surface(
      record.arg(2),
      boundary: [
        for (final loop in loops)
          for (final coedge in loop.coedges)
            for (var i = 0; i <= 4; i++)
              _edges[coedge.edge].curve.evaluate(
                _edges[coedge.edge].curve.t0 +
                    (_edges[coedge.edge].curve.t1 -
                            _edges[coedge.edge].curve.t0) *
                        i /
                        4,
              ),
      ],
    );
    if (surface == null) return null;

    if (surface is BrepSplineSurface && loops.isNotEmpty) {
      final u = surface.u0 + (surface.u1 - surface.u0) * .37,
          v = surface.v0 + (surface.v1 - surface.v0) * .53,
          point = surface.evaluate(u, v),
          normal = surface.normalAt(u, v);
      final plane = BrepPlane(origin: point, normal: normal);
      if (normal.length > .5 &&
          surface.points
              .expand((r) => r)
              .every((p) => plane.height(p).abs() < 1e-9))
        return BrepFace(
          surface: plane,
          loops: _planarLoops(loops, plane),
          reversed: reversed,
        );
    }
    if (surface is BrepSplineSurface && loops.isEmpty) {
      final native = BrepSplineConvert.ellipsoid(surface);
      if (native != null) {
        final u = surface.u0 + (surface.u1 - surface.u0) * .37,
            v = surface.v0 + (surface.v1 - surface.v0) * .53,
            p = surface.evaluate(u, v),
            uv = native.parameters(p)!;
        return BrepFace(
          surface: native,
          loops: [],
          reversed:
              reversed !=
              (surface.normalAt(u, v).dot(native.normalAt(uv.$1, uv.$2)) < 0),
        );
      }
    }
    return BrepFace(
      surface: surface,
      loops: surface is BrepPlane ? _planarLoops(loops, surface) : loops,
      reversed: reversed,
    );
  }

  // FACE_BOUND entries are an unordered set in STEP. In particular, many
  // exporters list an annulus's inner circle before its outer circle.
  List<BrepLoop> _planarLoops(List<BrepLoop> loops, BrepPlane plane) {
    if (loops.length < 2) return loops;
    double area(BrepLoop loop) {
      final points = <Vector3>[];
      for (final co in loop.coedges) {
        final c = _edges[co.edge].curve;
        final count = c is BrepLine ? 1 : 32;
        for (var i = 0; i < count; i++) {
          final fraction = co.reversed ? 1 - i / count : i / count;
          points.add(
            c.evaluate(c.t0 + (c.t1 - c.t0) * fraction) - plane.origin,
          );
        }
      }
      var total = 0.0;
      for (var i = 0; i < points.length; i++) {
        total += points[i]
            .cross(points[(i + 1) % points.length])
            .dot(plane.normal);
      }
      return total.abs();
    }

    final ordered = [for (final loop in loops) (loop, area(loop))]
      ..sort((a, b) => b.$2.compareTo(a.$2));
    return [for (final entry in ordered) entry.$1];
  }

  /// Whether [loop] is a seam rather than a trim.
  ///
  /// A whole sphere or torus has no boundary, but STEP requires every face to
  /// carry one, so both are written with a cut along the periodic parameter's
  /// zero line, traversed once each way. Such a loop encloses the entire
  /// surface and removes nothing — and taken literally it would, since the
  /// kernel treats a face's loops as trimming it, leaving a sphere measuring
  /// zero.
  ///
  /// The test is that *every* edge appears exactly twice with opposite senses.
  /// That distinguishes a seam from the ordinary case of a cylinder wall,
  /// whose loop also runs along a seam edge twice but closes over two rim
  /// circles that appear once each — and those do trim.
  static bool _isSeam(BrepLoop loop) {
    if (loop.coedges.isEmpty) return true;
    final senses = <int, int>{};
    for (final coedge in loop.coedges) {
      senses[coedge.edge] =
          (senses[coedge.edge] ?? 0) + (coedge.reversed ? -1 : 1);
    }
    if (senses.length * 2 != loop.coedges.length) return false;
    return senses.values.every((sum) => sum == 0);
  }

  BrepLoop? _bound(StepValue value) {
    final bound = value is StepRef
        ? _file.anyOf(value.id, ['FACE_OUTER_BOUND', 'FACE_BOUND'])
        : null;
    if (bound == null) return null;

    final orientation = bound.arg(2);
    final flip = orientation is StepEnum && orientation.isFalse;

    final loopRef = bound.arg(1);
    if (loopRef is! StepRef) return null;

    if (_file.typed(loopRef.id, 'VERTEX_LOOP') != null) {
      // A single degenerate point — a cone apex, a sphere pole. It bounds
      // nothing this kernel needs to store, since a face's loops here are
      // edge loops, and the surface's own degeneracy already handles it.
      return BrepLoop([]);
    }

    final poly = _file.typed(loopRef.id, 'POLY_LOOP');
    if (poly != null) return _polyLoop(poly, flip);

    final edgeLoop = _file.typed(loopRef.id, 'EDGE_LOOP');
    if (edgeLoop == null) {
      _count(_file.typeName(loopRef.id));
      return null;
    }

    final list = edgeLoop.arg(1);
    if (list is! StepList) return null;

    final coedges = <BrepCoedge>[];
    for (final item in list.items) {
      final oriented = _file.follow(item, 'ORIENTED_EDGE');
      if (oriented == null) return null;
      final edgeRef = oriented.arg(3);
      if (edgeRef is! StepRef) return null;
      final forward = oriented.arg(4);
      final reversed = forward is StepEnum && forward.isFalse;

      final index = _edge(edgeRef.id);
      if (index == null) return null;
      coedges.add(BrepCoedge(index, reversed: reversed));
    }
    if (coedges.isEmpty) return null;

    return BrepLoop(flip ? _reverseLoop(coedges) : coedges);
  }

  /// A faceted brep's loop, given as points rather than as shared topology.
  BrepLoop? _polyLoop(StepRecord poly, bool flip) {
    final list = poly.arg(1);
    if (list is! StepList) return null;
    final indices = <int>[];
    for (final item in list.items) {
      final p = _cartesian(item);
      if (p == null) return null;
      final index = _weld(p);
      if (indices.isEmpty || indices.last != index) indices.add(index);
    }
    if (indices.length > 1 && indices.first == indices.last) {
      indices.removeLast();
    }
    if (indices.length < 3) return null;

    final coedges = <BrepCoedge>[];
    for (var i = 0; i < indices.length; i++) {
      final a = indices[i];
      final b = indices[(i + 1) % indices.length];
      final index = _lineEdge(a, b);
      coedges.add(BrepCoedge(index, reversed: _edges[index].v0 != a));
    }
    return BrepLoop(flip ? _reverseLoop(coedges) : coedges);
  }

  static List<BrepCoedge> _reverseLoop(List<BrepCoedge> coedges) => [
    for (final c in coedges.reversed) BrepCoedge(c.edge, reversed: !c.reversed),
  ];

  int _lineEdge(int a, int b) {
    final key = a < b ? (a, b) : (b, a);
    final prior = _lineEdges[key];
    if (prior != null) return prior;
    _lineEdges[key] = _edges.length;
    _edges.add(
      BrepEdge(v0: a, v1: b, curve: BrepLine(_vertices[a], _vertices[b])),
    );
    return _edges.length - 1;
  }

  /// Maps an `EDGE_CURVE` to an edge index, reusing one already read.
  int? _edge(int id) {
    final existing = _edgeOfInstance[id];
    if (existing != null) return existing;

    final record = _file.typed(id, 'EDGE_CURVE');
    if (record == null) {
      _count(_file.typeName(id));
      return null;
    }

    final v0 = _vertex(record.arg(1));
    final v1 = _vertex(record.arg(2));
    if (v0 == null || v1 == null) return null;

    final sense = record.arg(4);
    final forward = !(sense is StepEnum && sense.isFalse);

    final curve = _curve(
      record.arg(3),
      _vertices[v0],
      _vertices[v1],
      forward: forward,
    );
    if (curve == null) return null;

    _edges.add(BrepEdge(v0: v0, v1: v1, curve: curve));
    final index = _edges.length - 1;
    _edgeOfInstance[id] = index;
    return index;
  }

  int? _vertex(StepValue? value) {
    if (value is! StepRef) return null;
    final existing = _vertexOfInstance[value.id];
    if (existing != null) return existing;

    final record = _file.typed(value.id, 'VERTEX_POINT');
    if (record == null) {
      _count(_file.typeName(value.id));
      return null;
    }
    final p = _cartesian(record.arg(1));
    if (p == null) return null;

    final index = _weld(p);
    _vertexOfInstance[value.id] = index;
    return index;
  }

  /// Adds a position, reusing a coincident one.
  ///
  /// STEP shares a vertex instance between the faces that meet there, so this
  /// mostly fires for faceted breps, which give positions per loop and would
  /// otherwise arrive as a bag of unconnected polygons. It also catches ends
  /// that are meant to be the same point and differ in the last bits.
  ///
  /// The neighbouring cells are checked as well as the point's own. A single
  /// rounded key is faster and wrong at the boundaries: two positions a
  /// nanometre apart either side of a cell edge would round differently and
  /// stay separate, leaving a hairline crack through the solid exactly where
  /// two faces were supposed to join.
  int _weld(Vector3 p) {
    final cx = (p.x * 1000).round();
    final cy = (p.y * 1000).round();
    final cz = (p.z * 1000).round();
    for (var dx = -1; dx <= 1; dx++) {
      for (var dy = -1; dy <= 1; dy++) {
        for (var dz = -1; dz <= 1; dz++) {
          final found = _vertexOfPosition['${cx + dx},${cy + dy},${cz + dz}'];
          if (found != null && _vertices[found].distanceToSquared(p) < 1e-12) {
            return found;
          }
        }
      }
    }
    _vertices.add(p);
    final index = _vertices.length - 1;
    _vertexOfPosition['$cx,$cy,$cz'] = index;
    return index;
  }

  Vector3? _cartesian(StepValue? value) {
    if (value is! StepRef) return null;
    final record = _file.typed(value.id, 'CARTESIAN_POINT');
    if (record == null) return null;
    final coords = record.arg(1);
    if (coords is! StepList || coords.items.length < 3) return null;
    final xyz = [
      for (final c in coords.items.take(3))
        c is StepNumber ? c.value : double.nan,
    ];
    if (xyz.any((v) => v.isNaN)) return null;
    return Vector3(xyz[0] * _scale, xyz[1] * _scale, xyz[2] * _scale);
  }

  Vector3? _dir(StepValue? value) {
    if (value is! StepRef) return null;
    final record = _file.typed(value.id, 'DIRECTION');
    if (record == null) return null;
    final coords = record.arg(1);
    if (coords is! StepList || coords.items.length < 3) return null;
    final xyz = [
      for (final c in coords.items.take(3))
        c is StepNumber ? c.value : double.nan,
    ];
    if (xyz.any((v) => v.isNaN)) return null;
    final v = Vector3(xyz[0], xyz[1], xyz[2]);
    return v.length < 1e-12 ? null : v.normalized();
  }

  /// An `AXIS2_PLACEMENT_3D` as origin, axis and reference direction.
  ///
  /// Both directions are optional in the schema. When the axis is absent it is
  /// +Z, and when the reference direction is absent any perpendicular will do
  /// — which is exactly what this kernel's own surfaces fall back to.
  ({Vector3 origin, Vector3 axis, Vector3 refDir})? _placement(
    StepValue? value,
  ) {
    if (value is! StepRef) return null;
    final record = _file.typed(value.id, 'AXIS2_PLACEMENT_3D');
    if (record == null) return null;
    final origin = _cartesian(record.arg(1));
    if (origin == null) return null;
    final axis = _dir(record.arg(2)) ?? Vector3(0, 0, 1);
    final refDir = _dir(record.arg(3)) ?? brepFrameU(axis);
    return (origin: origin, axis: axis, refDir: refDir);
  }

  double? _number(StepValue? value) {
    if (value is StepNumber) return value.value;
    if (value is StepTyped) return _number(value.single);
    return null;
  }

  List<double> _numbers(StepValue? value) {
    if (value is! StepList)
      throw const FormatException('Expected spline numeric list');
    return [
      for (final v in value.items)
        _number(v) ?? (throw const FormatException('Invalid spline number')),
    ];
  }

  List<double> _expandedKnots(StepValue? multiplicities, StepValue? values) {
    final m = _numbers(multiplicities), k = _numbers(values);
    if (m.length != k.length ||
        m.any((v) => v < 1 || v != v.roundToDouble() || v > 10000))
      throw const FormatException('Invalid spline knot multiplicities');
    return [
      for (var i = 0; i < k.length; i++)
        for (var j = 0; j < m[i]; j++) k[i],
    ];
  }

  List<Vector3> _controlPoints(StepValue? value) {
    if (value is! StepList)
      throw const FormatException('Expected spline control points');
    return [
      for (final v in value.items)
        _cartesian(v) ??
            (throw const FormatException('Invalid spline control point')),
    ];
  }

  BrepCurve? _untrimmedConic(int id) {
    final conic = _file.anyOf(id, ['CIRCLE', 'ELLIPSE']);
    if (conic == null) return null;
    final at = _placement(conic.arg(1)), major = _number(conic.arg(2));
    if (at == null || major == null) return null;
    if (conic.type == 'CIRCLE')
      return BrepCircle(
        center: at.origin,
        normal: at.axis,
        xDir: at.refDir,
        radius: major * _scale,
      );
    final minor = _number(conic.arg(3));
    if (minor == null) return null;
    return BrepEllipse(
      center: at.origin,
      normal: at.axis,
      xDir: at.refDir,
      rMajor: major * _scale,
      rMinor: minor * _scale,
    );
  }

  BrepSplineSurface? _splineSurface(int id) {
    final knots = _file.typed(id, 'B_SPLINE_SURFACE_WITH_KNOTS');
    if (knots == null) return null;
    final inherited = _file.typed(id, 'B_SPLINE_SURFACE');
    final base = inherited ?? knots,
        off = inherited == null ? 1 : 0,
        koff = inherited == null ? 8 : 0;
    final grid = base.arg(off + 2);
    if (grid is! StepList) throw const FormatException('Invalid spline grid');
    final weight = _file.typed(id, 'RATIONAL_B_SPLINE_SURFACE')?.arg(0);
    return BrepSplineSurface(
      degreeU: _number(base.arg(off))!.toInt(),
      degreeV: _number(base.arg(off + 1))!.toInt(),
      points: [for (final row in grid.items) _controlPoints(row)],
      knotsU: _expandedKnots(knots.arg(koff), knots.arg(koff + 2)),
      knotsV: _expandedKnots(knots.arg(koff + 1), knots.arg(koff + 3)),
      weights: weight is StepList
          ? [for (final row in weight.items) _numbers(row)]
          : null,
    );
  }

  BrepSplineCurve? _untrimmedSplineCurve(int id) {
    final knots = _file.typed(id, 'B_SPLINE_CURVE_WITH_KNOTS');
    if (knots == null) return null;
    final inherited = _file.typed(id, 'B_SPLINE_CURVE');
    final base = inherited ?? knots,
        off = inherited == null ? 1 : 0,
        koff = inherited == null ? 6 : 0;
    final weight = _file.typed(id, 'RATIONAL_B_SPLINE_CURVE')?.arg(0);
    return BrepSplineCurve(
      degree: _number(base.arg(off))!.toInt(),
      points: _controlPoints(base.arg(off + 1)),
      knots: _expandedKnots(knots.arg(koff), knots.arg(koff + 1)),
      weights: weight == null ? null : _numbers(weight),
    );
  }

  BrepSplineCurve? _splineCurve(int id, Vector3 p0, Vector3 p1, bool forward) {
    final c = _untrimmedSplineCurve(id);
    if (c == null) return null;
    var a = c.closestParam(p0), b = c.closestParam(p1);
    if (c.isClosed) {
      final period = c.knots[c.points.length] - c.knots[c.degree];
      if ((p0 - p1).length < 1e-8) {
        b = a + (forward ? period : -period);
      } else if (forward && b < a) {
        b += period;
      } else if (!forward && b > a) {
        b -= period;
      }
    }
    if ((c.evaluate(a) - p0).length > _lengthTolerance ||
        (c.evaluate(b) - p1).length > _lengthTolerance) {
      throw FormatException('Spline edge #$id endpoints are off its curve');
    }
    return c.trimmed(a, b);
  }

  BrepSurface? _surface(StepValue? value, {List<Vector3> boundary = const []}) {
    if (value is! StepRef) return null;
    final id = value.id;
    final spline = _splineSurface(id);
    if (spline != null) return spline;

    // Unwrap the wrappers a real file puts around a perfectly good surface.
    // The face's own bounds already say where the surface is used, so a
    // parametric trim on top of it carries no information this kernel needs.
    final trimmed = _file.typed(id, 'RECTANGULAR_TRIMMED_SURFACE');
    if (trimmed != null) return _surface(trimmed.arg(1), boundary: boundary);
    final bounded = _file.typed(id, 'CURVE_BOUNDED_SURFACE');
    if (bounded != null) return _surface(bounded.arg(1), boundary: boundary);

    final offset = _file.typed(id, 'OFFSET_SURFACE');
    if (offset != null) {
      final base = _surface(offset.arg(1), boundary: boundary);
      final distance = _number(offset.arg(2));
      if (base == null || distance == null) return null;
      final d = distance * _scale;
      if (base is BrepPlane) return base.offset(d);
      if (base is BrepCylinder) {
        if (base.radius + d <= 0)
          throw const FormatException('Collapsed offset cylinder');
        return BrepCylinder(
          origin: base.origin,
          axis: base.axis,
          xDir: base.xDir,
          radius: base.radius + d,
        );
      }
      if (base is BrepSphere) {
        if (base.radius + d <= 0)
          throw const FormatException('Collapsed offset sphere');
        return BrepSphere(center: base.center, radius: base.radius + d);
      }
      if (base is BrepTorus) {
        if (base.minorRadius + d <= 0)
          throw const FormatException('Collapsed offset torus');
        return BrepTorus(
          center: base.center,
          axis: base.axis,
          xDir: base.xDir,
          majorRadius: base.majorRadius,
          minorRadius: base.minorRadius + d,
        );
      }
      if (base is BrepCone) {
        final shift =
            d *
            math.sqrt(
              1 + math.pow((base.topRadius - base.radius) / base.height, 2),
            );
        return BrepCone(
          origin: base.origin,
          axis: base.axis,
          xDir: base.xDir,
          radius: base.radius + shift,
          topRadius: base.topRadius + shift,
          height: base.height,
        );
      }
      return BrepOffsetSurface(base, d);
    }
    final revolution = _file.typed(id, 'SURFACE_OF_REVOLUTION');
    if (revolution != null) {
      final basis = revolution.arg(1);
      final placement = _file.follow(revolution.arg(2), 'AXIS1_PLACEMENT');
      final origin = _cartesian(placement?.arg(1));
      final axis = _dir(placement?.arg(2));
      if (basis is! StepRef || origin == null || axis == null) return null;
      BrepSplineCurve? curve = _untrimmedSplineCurve(basis.id);
      final conic = _untrimmedConic(basis.id);
      if (conic != null) curve = BrepSplineConvert.curve(conic);
      if (curve == null) {
        _count(_file.typeName(basis.id));
        return null;
      }
      final ring = BrepSplineConvert.curve(
        BrepCircle(
          center: Vector3.zero(),
          normal: Vector3(0, 0, 1),
          xDir: Vector3(1, 0, 0),
          radius: 1,
        ),
      );
      return BrepSplineSurface(
        degreeU: ring.degree,
        degreeV: curve.degree,
        knotsU: ring.knots,
        knotsV: curve.knots,
        points: [
          for (final r in ring.points)
            [
              for (final p in curve.points)
                origin +
                    axis * (p - origin).dot(axis) +
                    ((p - origin) - axis * (p - origin).dot(axis)) * r.x +
                    axis.cross(p - origin) * r.y,
            ],
        ],
        weights: [
          for (final w in ring.weights) [for (final v in curve.weights) w * v],
        ],
      );
    }

    final extrusion = _file.typed(id, 'SURFACE_OF_LINEAR_EXTRUSION');
    if (extrusion != null) {
      final basis = extrusion.arg(1);
      final vector = _file.follow(extrusion.arg(2), 'VECTOR');
      final direction = _dir(vector?.arg(1));
      if (basis is! StepRef || direction == null || boundary.isEmpty)
        return null;
      BrepCurve? curve = _untrimmedSplineCurve(basis.id);
      final conic = _untrimmedConic(basis.id);
      if (conic != null)
        curve = conic is BrepEllipse ? BrepSplineConvert.curve(conic) : conic;
      if (curve == null) {
        _count(_file.typeName(basis.id));
        return null;
      }
      // Positive-weight spline control points bound the entire curve. A
      // circle's exact projected extrema avoid underestimating tilted bases.
      final cp = curve is BrepSplineCurve
          ? curve.points.map((p) => p.dot(direction)).toList()
          : [
              (curve as BrepCircle).center.dot(direction) -
                  curve.radius *
                      math.sqrt(
                        math.max(
                          0,
                          1 - math.pow(curve.normal.dot(direction), 2),
                        ),
                      ),
              curve.center.dot(direction) +
                  curve.radius *
                      math.sqrt(
                        math.max(
                          0,
                          1 - math.pow(curve.normal.dot(direction), 2),
                        ),
                      ),
            ];
      final bp = boundary.map((p) => p.dot(direction)).toList();
      final lo = bp.reduce(math.min) - cp.reduce(math.max) - 1;
      final hi = bp.reduce(math.max) - cp.reduce(math.min) + 1;
      BrepCurve shifted(double distance) {
        final m = BrepAffineTransform(
          Vector3(1, 0, 0),
          Vector3(0, 1, 0),
          Vector3(0, 0, 1),
          direction * distance,
        );
        return curve is BrepSplineCurve
            ? curve.transformed(m)
            : BrepAffineCurve(curve!, m);
      }

      return BrepRuledSurface(shifted(lo), shifted(hi));
    }

    final plane = _file.typed(id, 'PLANE');
    if (plane != null) {
      final at = _placement(plane.arg(1));
      if (at == null) return null;
      return BrepPlane(origin: at.origin, normal: at.axis, xDir: at.refDir);
    }

    final cylinder = _file.typed(id, 'CYLINDRICAL_SURFACE');
    if (cylinder != null) {
      final at = _placement(cylinder.arg(1));
      final radius = _number(cylinder.arg(2));
      if (at == null || radius == null) return null;
      return BrepCylinder(
        origin: at.origin,
        axis: at.axis,
        radius: radius * _scale,
        xDir: at.refDir,
      );
    }

    final sphere = _file.typed(id, 'SPHERICAL_SURFACE');
    if (sphere != null) {
      final at = _placement(sphere.arg(1));
      final radius = _number(sphere.arg(2));
      if (at == null || radius == null) return null;
      // This kernel's sphere is axis-free, which costs nothing: a sphere is
      // the same shape however it is turned.
      return BrepSphere(center: at.origin, radius: radius * _scale);
    }

    final torus = _file.typed(id, 'TOROIDAL_SURFACE');
    if (torus != null) {
      final at = _placement(torus.arg(1));
      final major = _number(torus.arg(2));
      final minor = _number(torus.arg(3));
      if (at == null || major == null || minor == null) return null;
      return BrepTorus(
        center: at.origin,
        axis: at.axis,
        majorRadius: major * _scale,
        minorRadius: minor * _scale,
        xDir: at.refDir,
      );
    }

    final cone = _file.typed(id, 'CONICAL_SURFACE');
    if (cone != null) {
      final at = _placement(cone.arg(1));
      final radius = _number(cone.arg(2));
      final semiAngle = _number(cone.arg(3));
      if (at == null || radius == null || semiAngle == null) return null;
      // STEP's cone is infinite: a radius at the placement plane and a
      // half-angle. This kernel's is a frustum, so any height describes the
      // same cone as long as the top radius follows the angle. One unit keeps
      // the numbers legible.
      return BrepCone(
        origin: at.origin,
        axis: at.axis,
        radius: radius * _scale,
        topRadius: radius * _scale + math.tan(semiAngle),
        height: 1,
        xDir: at.refDir,
      );
    }

    _count(_file.typeName(id));
    return null;
  }

  // Parent surfaces define the locus; the 3D representation selects its finite
  // branch. Verify the entire sampled branch before recovering native geometry.
  BrepIntersectionCurve? _recoverIntersection(
    BrepCylinder carrier,
    BrepSurface target,
    BrepCurve representation,
  ) {
    if (target is BrepCylinder && carrier.axis.cross(target.axis).length < 1e-8)
      return null;
    final points = [
      for (var i = 0; i <= 32; i++)
        representation.evaluate(
          representation.t0 + (representation.t1 - representation.t0) * i / 32,
        ),
    ];
    final parameters = <double>[];
    for (final p in points) {
      var u = carrier.parameters(p)!.$1;
      if (parameters.isNotEmpty) {
        while (u - parameters.last > math.pi) u -= 2 * math.pi;
        while (u - parameters.last < -math.pi) u += 2 * math.pi;
      }
      parameters.add(u);
    }
    if (parameters.last <= parameters.first + 1e-12) return null;
    for (var i = 1; i < parameters.length; i++) {
      if (parameters[i] < parameters[i - 1] - 1e-8) return null;
    }
    for (final branch in [-1, 1]) {
      final exact = BrepIntersectionCurve(
        carrier: carrier,
        target: target,
        branch: branch,
        t0: parameters.first,
        t1: parameters.last,
      );
      try {
        if (List.generate(points.length, (i) => i).every(
          (i) => (exact.evaluate(parameters[i]) - points[i]).length < 1e-6,
        ))
          return exact;
      } on StateError {
        continue;
      }
    }
    return null;
  }

  BrepSurfaceIntersectionCurve? _recoverSplineIntersection(
    BrepSplineSurface a,
    BrepSplineSurface b,
    BrepCurve representation,
  ) {
    // Refine branch anchors until between-anchor evaluations also agree with
    // the supplied 3D curve. A fixed anchor count can misparameterize short,
    // strongly curved or unevenly parameterized imported edges.
    for (var count = 64; count <= 1024; count *= 2) {
      final anchors = <Vector4>[];
      try {
        for (var i = 0; i <= count; i++) {
          final p = representation.evaluate(
                representation.t0 +
                    (representation.t1 - representation.t0) * i / count,
              ),
              u = a.parameters(p),
              v = b.parameters(p);
          final x = Vector4(
            (u.$1 - a.u0) / (a.u1 - a.u0),
            (u.$2 - a.v0) / (a.v1 - a.v0),
            (v.$1 - b.u0) / (b.u1 - b.u0),
            (v.$2 - b.v0) / (b.v1 - b.v0),
          );
          final corrected = BrepSurfaceIntersectionCurve.correct(
            a,
            b,
            x,
            x,
            BrepSurfaceIntersectionCurve.direction(a, b, x),
            1e-9,
          );
          if (corrected == null ||
              (BrepSurfaceIntersectionCurve.point(a, corrected.x, corrected.y) -
                          p)
                      .length >
                  1e-6) {
            IjoLog.line(
              'STEP intersection recovery sample $i/$count failed: p=$p uv=$x corrected=$corrected',
            );
            return null;
          }
          anchors.add(corrected);
        }
        final curve = BrepSurfaceIntersectionCurve(
          first: a,
          second: b,
          anchors: anchors,
        );
        var matches = true;
        for (var i = 0; i < count; i++) {
          var t = (i + .5) / count;
          final p = representation.evaluate(
            representation.t0 + (representation.t1 - representation.t0) * t,
          );
          final speed =
              (curve.evaluate((i + 1) / count) - curve.evaluate(i / count))
                  .length *
              count;
          // STEP does not require the recovered curve to have the supplied
          // 3D curve's parameterization. Verify geometric agreement within
          // the corresponding bounded anchor interval instead.
          for (var k = 0; k < 8; k++) {
            final q = curve.evaluate(t);
            if ((q - p).length <= 1e-6 || speed < 1e-14) break;
            t = (t - (q - p).dot(curve.tangent(t)) / speed).clamp(
              i / count,
              (i + 1) / count,
            );
          }
          if ((curve.evaluate(t) - p).length > 1e-6) {
            matches = false;
            break;
          }
        }
        if (matches) return curve;
      } on BrepUnsupported catch (e) {
        IjoLog.line('STEP intersection recovery rejected: $e');
        return null;
      }
    }
    IjoLog.line(
      'STEP intersection recovery could not match branch parameterization',
    );
    return null;
  }

  BrepCurve? _curve(
    StepValue? value,
    Vector3 p0,
    Vector3 p1, {
    required bool forward,
  }) {
    if (value is! StepRef) return null;
    final id = value.id;
    final spline = _splineCurve(id, p0, p1, forward);
    if (spline != null) return spline;

    final intersection = _file.typed(id, 'INTERSECTION_CURVE');
    if (intersection != null) {
      final representation = _curve(
        intersection.arg(1),
        p0,
        p1,
        forward: forward,
      );
      if (representation == null) return null;
      final associated = intersection.arg(2);
      if (associated is StepList && associated.items.length == 2) {
        BrepSurface? parent(StepValue value) {
          final pcurve = _file.follow(value, 'PCURVE');
          return _surface(pcurve?.arg(1) ?? value);
        }

        final a = parent(associated.items[0]), b = parent(associated.items[1]);
        if (a is BrepSplineSurface && b is BrepSplineSurface) {
          final exact = _recoverSplineIntersection(a, b, representation);
          if (exact != null) return exact;
        }
        if (a is BrepCylinder && (b is BrepCylinder || b is BrepSphere)) {
          final exact = _recoverIntersection(a, b!, representation);
          if (exact != null) return exact;
        }
        if (b is BrepCylinder && (a is BrepCylinder || a is BrepSphere)) {
          final exact = _recoverIntersection(b, a!, representation);
          if (exact != null) return exact;
        }
      }
      return representation;
    }

    // As with surfaces: unwrap. `SURFACE_CURVE` and `SEAM_CURVE` are how a
    // real file attaches parameter-space curves to a 3D one, and they are
    // everywhere — the 3D curve inside is the one that matters here.
    final surfaceCurve = _file.anyOf(id, ['SEAM_CURVE', 'SURFACE_CURVE']);
    if (surfaceCurve != null) {
      return _curve(surfaceCurve.arg(1), p0, p1, forward: forward);
    }
    final trimmed = _file.typed(id, 'TRIMMED_CURVE');
    if (trimmed != null) {
      final sense = trimmed.arg(4);
      return _curve(
        trimmed.arg(1),
        p0,
        p1,
        forward: forward != (sense is StepEnum && sense.isFalse),
      );
    }

    if (_file.typed(id, 'LINE') != null) {
      // Taken from the vertices rather than from the LINE's own point and
      // vector, since those describe an unbounded line and the edge is the
      // piece between its ends.
      if (p0.distanceToSquared(p1) < 1e-18) return null;
      return BrepLine(p0, p1);
    }

    final circle = _file.typed(id, 'CIRCLE');
    if (circle != null) {
      final at = _placement(circle.arg(1));
      final radius = _number(circle.arg(2));
      if (at == null || radius == null) return null;
      final (t0, t1) = _arcRange(
        at.origin,
        at.axis,
        at.refDir,
        p0,
        p1,
        ascending: forward,
        scaleU: radius * _scale,
        scaleV: radius * _scale,
      );
      return BrepCircle(
        center: at.origin,
        normal: at.axis,
        xDir: at.refDir,
        radius: radius * _scale,
        t0: t0,
        t1: t1,
      );
    }

    final hyperbola = _file.typed(id, 'HYPERBOLA');
    if (hyperbola != null) {
      final at = _placement(hyperbola.arg(1));
      final major = _number(hyperbola.arg(2));
      final minor = _number(hyperbola.arg(3));
      if (at == null ||
          major == null ||
          minor == null ||
          major <= 0 ||
          minor <= 0)
        return null;
      final x = at.refDir.normalized(),
          y = at.axis.cross(at.refDir).normalized();
      final a = major * _scale, b = minor * _scale;
      double parameter(Vector3 p) {
        final v = (p - at.origin).dot(y) / b;
        return v.sign * math.log(v.abs() + math.sqrt(v * v + 1));
      }

      double cosh(double t) => (math.exp(t) + math.exp(-t)) / 2;
      double sinh(double t) => (math.exp(t) - math.exp(-t)) / 2;
      Vector3 point(double t) =>
          at.origin + x * (a * cosh(t)) + y * (b * sinh(t));
      final start = parameter(p0), end = parameter(p1);
      if ((point(start) - p0).length > _lengthTolerance ||
          (point(end) - p1).length > _lengthTolerance) {
        throw const FormatException(
          'Hyperbola edge endpoints are off its curve',
        );
      }
      final count = math.max(1, (end - start).abs().ceil());
      final points = <Vector3>[point(start)], weights = <double>[1];
      for (var i = 0; i < count; i++) {
        final lo = start + (end - start) * i / count;
        final hi = start + (end - start) * (i + 1) / count;
        final w = cosh((hi - lo) / 2);
        points.add(at.origin + (point((lo + hi) / 2) - at.origin) / w);
        points.add(point(hi));
        weights.addAll([w, 1]);
      }
      return BrepSplineCurve(
        degree: 2,
        points: points,
        weights: weights,
        knots: [
          0,
          0,
          0,
          for (var i = 1; i < count; i++) ...[i.toDouble(), i.toDouble()],
          count.toDouble(),
          count.toDouble(),
          count.toDouble(),
        ],
      );
    }

    final ellipse = _file.typed(id, 'ELLIPSE');
    if (ellipse != null) {
      final at = _placement(ellipse.arg(1));
      final major = _number(ellipse.arg(2));
      final minor = _number(ellipse.arg(3));
      if (at == null || major == null || minor == null) return null;
      final (t0, t1) = _arcRange(
        at.origin,
        at.axis,
        at.refDir,
        p0,
        p1,
        ascending: forward,
        scaleU: major * _scale,
        scaleV: minor * _scale,
      );
      return BrepEllipse(
        center: at.origin,
        normal: at.axis,
        xDir: at.refDir,
        rMajor: major * _scale,
        rMinor: minor * _scale,
        t0: t0,
        t1: t1,
      );
    }

    _count(_file.typeName(id));
    return null;
  }

  /// The parameter range of the arc from [p0] to [p1].
  ///
  /// STEP states a whole circle or ellipse and leaves the edge's two vertices
  /// to say which piece of it is meant. This kernel does the opposite: the
  /// range lives on the curve. Reconciling them means measuring the angle of
  /// each endpoint in the curve's own frame.
  ///
  /// [ascending] is the edge's `same_sense`, and it decides which way round
  /// the arc goes. A STEP circle always runs anticlockwise about its axis, so
  /// an edge marked against the curve is a clockwise sweep — which this kernel
  /// writes as a **descending** range, `t1 < t0`, rather than by flipping the
  /// axis. Keeping the frame as the file states it is what makes the round
  /// trip exact; flipping it instead moves the same arc to a different
  /// parametrisation, and a face measured over that domain comes out wrong.
  (double, double) _arcRange(
    Vector3 center,
    Vector3 axis,
    Vector3 refDir,
    Vector3 p0,
    Vector3 p1, {
    required bool ascending,
    required double scaleU,
    required double scaleV,
  }) {
    final x = refDir - axis * refDir.dot(axis);
    final xDir = x.length < 1e-12 ? brepFrameU(axis) : x.normalized();
    final yDir = axis.cross(xDir).normalized();
    const twoPi = 2 * math.pi;

    double angle(Vector3 p) {
      final d = p - center;
      final a = math.atan2(d.dot(yDir) / scaleV, d.dot(xDir) / scaleU);
      return a < 0 ? a + twoPi : a;
    }

    final t0 = angle(p0);

    // Coincident ends mean the whole curve, which is how both a closed rim and
    // a synthesised seam arrive.
    if (p0.distanceToSquared(p1) < 1e-12) {
      return ascending ? (t0, t0 + twoPi) : (t0, t0 - twoPi);
    }

    var t1 = angle(p1);
    if (ascending && t1 <= t0) t1 += twoPi;
    if (!ascending && t1 >= t0) t1 -= twoPi;
    return (t0, t1);
  }

  void _count(String type) {
    _unsupported[type] = (_unsupported[type] ?? 0) + 1;
  }

  /// What one file unit is in millimetres, and what to call it.
  ///
  /// A STEP file states its units, and getting this wrong is the classic way
  /// to import a part at 25.4 times its size. The declaration hangs off the
  /// geometric context; when a file has none, millimetres is the convention.
  static (double, String) _lengthUnit(StepFile file) {
    for (final id in file.idsOf('GLOBAL_UNIT_ASSIGNED_CONTEXT')) {
      final record = file.typed(id, 'GLOBAL_UNIT_ASSIGNED_CONTEXT')!;
      final units = record.arg(0);
      if (units is! StepList) continue;
      for (final item in units.items) {
        if (item is! StepRef) continue;
        final found = _lengthOf(file, item.id);
        if (found != null) return found;
      }
    }
    return (1.0, 'mm');
  }

  static (double, String)? _lengthOf(StepFile file, int id) {
    if (file.typed(id, 'LENGTH_UNIT') == null) return null;

    final si = file.typed(id, 'SI_UNIT');
    if (si != null) {
      final name = si.arg(1);
      if (name is StepEnum && name.name == 'METRE') {
        final prefix = si.arg(0);
        final factor = prefix is StepEnum ? _prefix(prefix.name) : 1.0;
        // The kernel works in millimetres, so a metre is 1000 of them.
        return (
          factor * 1000,
          _prefixName(prefix is StepEnum ? prefix.name : null),
        );
      }
      return null;
    }

    // An imperial file: CONVERSION_BASED_UNIT gives a name and a factor
    // expressed in an SI unit.
    final conversion = file.typed(id, 'CONVERSION_BASED_UNIT');
    if (conversion != null) {
      final label = conversion.arg(0);
      final factorRef = conversion.arg(1);
      final measure =
          file.follow(factorRef, 'LENGTH_MEASURE_WITH_UNIT') ??
          file.follow(factorRef, 'MEASURE_WITH_UNIT');
      if (measure != null) {
        final raw = measure.arg(0);
        final value = raw is StepTyped
            ? (raw.single is StepNumber
                  ? (raw.single as StepNumber).value
                  : null)
            : (raw is StepNumber ? raw.value : null);
        final baseRef = measure.arg(1);
        final base = baseRef is StepRef ? _lengthOf(file, baseRef.id) : null;
        if (value != null && base != null) {
          return (
            value * base.$1,
            label is StepText && label.value.isNotEmpty ? label.value : 'unit',
          );
        }
      }
    }
    return null;
  }

  static double _prefix(String name) => switch (name) {
    'EXA' => 1e18,
    'PETA' => 1e15,
    'TERA' => 1e12,
    'GIGA' => 1e9,
    'MEGA' => 1e6,
    'KILO' => 1e3,
    'HECTO' => 1e2,
    'DECA' => 1e1,
    'DECI' => 1e-1,
    'CENTI' => 1e-2,
    'MILLI' => 1e-3,
    'MICRO' => 1e-6,
    'NANO' => 1e-9,
    'PICO' => 1e-12,
    'FEMTO' => 1e-15,
    'ATTO' => 1e-18,
    _ => 1.0,
  };

  static String _prefixName(String? prefix) => switch (prefix) {
    'MILLI' => 'mm',
    'CENTI' => 'cm',
    'DECI' => 'dm',
    'KILO' => 'km',
    'MICRO' => 'µm',
    null => 'm',
    _ => 'm',
  };
}
