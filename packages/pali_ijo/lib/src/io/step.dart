/// STEP import/export for BREP solids and explicitly faceted mesh bodies.
/// Analytic parts retain their surfaces. The mesh read API derives triangles
/// from the same BREP reader, while readSolid/readSolids retain the topology.
library;

import 'dart:typed_data';

import 'package:vector_math/vector_math_64.dart';

import '../brep/solid.dart';
import '../brep/primitives.dart';
import '../brep/tessellate.dart';
import '../mesh/mesh_data.dart';
import 'solid_json.dart';
import 'step_analytic.dart';
export 'step_analytic.dart'
    show StepBrepUnsupportedException, StepBrepReadResult, StepSyntaxException;
import 'step_read.dart';
import 'step_syntax.dart';

export 'step_read.dart' show StepUnsupportedException, StepReadResult;
export 'step_syntax.dart' show StepSyntaxException;

class StepPart {
  const StepPart(MeshData mesh, {this.name}) : _mesh = mesh, solid = null;
  const StepPart.solid(this.solid, {this.name}) : _mesh = null;

  final MeshData? _mesh;
  final BrepSolid? solid;
  MeshData get mesh => _mesh ?? BrepTessellate.run(solid!);
  final String? name;
}

class IjoStep {
  IjoStep._();

  static const schema = 'AUTOMOTIVE_DESIGN { 1 0 10303 214 1 1 1 1 }';

  static Uint8List writeMesh(
    MeshData mesh, {
    String name = 'Part',
    DateTime? timestamp,
  }) => write(mesh, name: name, timestamp: timestamp);

  static Uint8List write(
    MeshData mesh, {
    String name = 'Part',
    DateTime? timestamp,
  }) => writeAll([StepPart(mesh, name: name)], timestamp: timestamp);

  static Uint8List writeAll(
    List<StepPart> parts, {
    String fileName = 'part.step',
    DateTime? timestamp,
  }) {
    if (parts.isEmpty) {
      throw ArgumentError.value(parts, 'parts', 'nothing to write');
    }
    if (parts.any((part) => part.solid != null)) {
      return IjoStepBrep.writeAll(
        [
          for (final part in parts)
            StepBrepPart(
              part.solid ??
                  BrepPrimitives.tryPlanarSolid(part.mesh) ??
                  (throw BrepUnsupported(
                    'Mesh body has no closed planar boundary',
                  )),
              name: part.name,
            ),
        ],
        fileName: fileName,
        timestamp: timestamp,
      );
    }
    final w = StepWriter();

    final appContext = w.add(
      'APPLICATION_CONTEXT',
      StepWriter.text('core data for automotive mechanical design processes'),
    );
    w.add(
      'APPLICATION_PROTOCOL_DEFINITION',
      '${StepWriter.text('international standard')},'
          '${StepWriter.text('automotive_design')},2000,#$appContext',
    );
    final lengthUnit = w.addComplex([
      'NAMED_UNIT(*)',
      'SI_UNIT(.MILLI.,.METRE.)',
      'LENGTH_UNIT()',
    ]);
    final angleUnit = w.addComplex([
      'NAMED_UNIT(*)',
      'PLANE_ANGLE_UNIT()',
      r'SI_UNIT($,.RADIAN.)',
    ]);
    final solidAngleUnit = w.addComplex([
      'NAMED_UNIT(*)',
      r'SI_UNIT($,.STERADIAN.)',
      'SOLID_ANGLE_UNIT()',
    ]);
    final uncertainty = w.add(
      'UNCERTAINTY_MEASURE_WITH_UNIT',
      'LENGTH_MEASURE(1.E-07),#$lengthUnit,'
          '${StepWriter.text('distance_accuracy_value')},'
          '${StepWriter.text('confusion accuracy')}',
    );
    final context = w.addComplex([
      'GEOMETRIC_REPRESENTATION_CONTEXT(3)',
      'GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT((#$uncertainty))',
      'GLOBAL_UNIT_ASSIGNED_CONTEXT('
          '(#$lengthUnit,#$angleUnit,#$solidAngleUnit))',
      'REPRESENTATION_CONTEXT(${StepWriter.text('')},'
          '${StepWriter.text('3D')})',
    ]);

    final bodies = <int>[];
    for (var i = 0; i < parts.length; i++) {
      final part = parts[i];
      final name = part.name ?? (parts.length == 1 ? 'Part' : 'Part ${i + 1}');
      bodies.add(_writeSolid(w, part.mesh, name));
    }

    final origin = _placement(
      w,
      Vector3.zero(),
      Vector3(0, 0, 1),
      Vector3(1, 0, 0),
    );
    final shape = w.add(
      'FACETED_BREP_SHAPE_REPRESENTATION',
      '${StepWriter.text('')},'
          '${StepWriter.refs([origin, ...bodies])},#$context',
    );

    final productContext = w.add(
      'PRODUCT_CONTEXT',
      '${StepWriter.text('')},#$appContext,'
          '${StepWriter.text('mechanical')}',
    );
    final definitionContext = w.add(
      'PRODUCT_DEFINITION_CONTEXT',
      '${StepWriter.text('part definition')},#$appContext,'
          '${StepWriter.text('design')}',
    );
    final product = w.add(
      'PRODUCT',
      '${StepWriter.text(parts.first.name ?? 'Part')},'
          '${StepWriter.text(parts.first.name ?? 'Part')},'
          '${StepWriter.text('')},(#$productContext)',
    );
    final formation = w.add(
      'PRODUCT_DEFINITION_FORMATION',
      '${StepWriter.text('')},${StepWriter.text('')},#$product',
    );
    final definition = w.add(
      'PRODUCT_DEFINITION',
      '${StepWriter.text('design')},${StepWriter.text('')},'
          '#$formation,#$definitionContext',
    );
    final shapeDef = w.add(
      'PRODUCT_DEFINITION_SHAPE',
      '${StepWriter.text('')},${StepWriter.text('')},#$definition',
    );
    w.add('SHAPE_DEFINITION_REPRESENTATION', '#$shapeDef,#$shape');

    return Uint8List.fromList(
      w.bytes(
        description: 'pali ijo faceted solid',
        fileName: fileName,
        originator: 'pali ijo',
        schema: schema,
        timestamp: timestamp,
      ),
    );
  }

  static int _writeSolid(StepWriter w, MeshData mesh, String name) {
    final points = <int>[];
    for (final v in mesh.vertices) {
      points.add(
        w.dedupe(
          'CARTESIAN_POINT',
          '${StepWriter.text('')},(${StepWriter.num(v.x)},${StepWriter.num(v.y)},${StepWriter.num(v.z)})',
        ),
      );
    }
    final faces = <int>[];
    for (final f in mesh.faces) {
      if (f.length < 3) continue;
      final loop = w.add(
        'POLY_LOOP',
        '${StepWriter.text('')},${StepWriter.refs([for (final i in f) points[i]])}',
      );
      final bound = w.add(
        'FACE_OUTER_BOUND',
        '${StepWriter.text('')},#$loop,.T.',
      );
      faces.add(w.add('FACE', '${StepWriter.text('')},(#$bound)'));
    }
    final shell = w.add(
      'CLOSED_SHELL',
      '${StepWriter.text('')},${StepWriter.refs(faces)}',
    );
    return w.add('FACETED_BREP', '${StepWriter.text(name)},#$shell');
  }

  static int _placement(StepWriter w, Vector3 o, Vector3 z, Vector3 x) {
    final origin = w.dedupe(
      'CARTESIAN_POINT',
      '${StepWriter.text('')},(${StepWriter.num(o.x)},${StepWriter.num(o.y)},${StepWriter.num(o.z)})',
    );
    final axis = w.dedupe(
      'DIRECTION',
      '${StepWriter.text('')},(${StepWriter.num(z.x)},${StepWriter.num(z.y)},${StepWriter.num(z.z)})',
    );
    final ref = w.dedupe(
      'DIRECTION',
      '${StepWriter.text('')},(${StepWriter.num(x.x)},${StepWriter.num(x.y)},${StepWriter.num(x.z)})',
    );
    return w.add(
      'AXIS2_PLACEMENT_3D',
      '${StepWriter.text('')},#$origin,#$axis,#$ref',
    );
  }

  static Uint8List writeSolid(
    BrepSolid solid, {
    String name = 'Part',
    DateTime? timestamp,
  }) => IjoStepBrep.write(solid, name: name, timestamp: timestamp);

  static BrepSolid readSolid(Uint8List bytes) => IjoStepBrep.read(bytes);

  /// By default, incomplete shells are rejected. [allowOpenShells] retains
  /// readable open bodies for inspection and reports them in the result notes;
  /// unsupported or missing faces are still rejected.
  /// [repairDegenerateFaces] opts into capping isolated openings of collapsed
  /// cylindrical faces and sewing their retraced wires. This changes geometry;
  /// callers must present the returned repair notes to the user.
  /// [onProgress] counts source bodies after syntax parsing; assembly instances
  /// are placed afterwards. It is not a time estimate for the whole import.
  static StepBrepReadResult readSolids(
    Uint8List bytes, {
    bool allowOpenShells = false,
    bool repairDegenerateFaces = false,
    void Function(int completed, int total)? onProgress,
  }) => IjoStepBrep.readAll(
    bytes,
    allowOpenShells: allowOpenShells,
    repairDegenerateFaces: repairDegenerateFaces,
    onProgress: onProgress,
  );

  static MeshData read(Uint8List bytes) {
    try {
      final analytic = IjoStepBrep.readAll(bytes);
      if (analytic.solids.isNotEmpty) {
        return BrepTessellate.run(analytic.solids.first);
      }
    } catch (_) {}
    return readAll(bytes).solids.first;
  }

  static StepReadResult readAll(Uint8List bytes) {
    final result = readSolids(bytes);
    return StepReadResult(
      solids: [for (final solid in result.solids) BrepTessellate.run(solid)],
      names: result.names,
      unitName: result.unitName,
      unitScale: result.unitScale,
      notes: result.notes,
    );
  }

  static Object writableSolid(IjoDocument document) =>
      document.solidOrNull ?? document.mesh;
}
