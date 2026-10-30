import 'dart:convert';
import 'dart:typed_data';

import 'package:vector_math/vector_math_64.dart';
import 'package:xml/xml_events.dart';

import '../mesh/mesh_data.dart';
import 'zip.dart';

part 'threemf_read.dart';

/// One mesh object inside a 3MF file. The core spec allows many of these
/// in one archive, the same way STEP can hold several solids.
class ThreeMfPart {
  const ThreeMfPart(this.mesh, {this.name});

  final MeshData mesh;
  final String? name;
}

/// 3D Manufacturing Format — a ZIP of XML, and the right default for solids.
///
/// This is the format to use when a mesh will be read again. Against STL it
/// wins on every count that matters to a kernel: vertices are shared so
/// topology survives, coordinates are decimal text rather than float32, and
/// the file states its units instead of leaving millimetres to be guessed.
/// It is a *lossless mesh*
/// relative to the tessellation — not analytic CAD (STEP's job).
///
/// Reads core meshes, build instances and production-extension references,
/// preserving component transforms. CAD surfaces belong in STEP; 3MF stays mesh.
class IjoThreeMf {
  IjoThreeMf._();

  static const _modelPath = '3D/3dmodel.model';
  static const _modelNamespace =
      'http://schemas.microsoft.com/3dmanufacturing/core/2015/02';
  static const _relationshipType =
      'http://schemas.microsoft.com/3dmanufacturing/2013/01/3dmodel';

  static const _units = {
    'micron',
    'millimeter',
    'centimeter',
    'inch',
    'foot',
    'meter',
  };

  /// Writes a single-object 3MF.
  ///
  /// [unit] is one of the units the specification allows: `micron`,
  /// `millimeter`, `centimeter`, `inch`, `foot`, `meter`. The kernel works in
  /// millimetres, so that is the default and callers should have a reason to
  /// change it.
  static Uint8List write(
    MeshData mesh, {
    String name = 'pali_ijo',
    String unit = 'millimeter',
  }) => writeAll([ThreeMfPart(mesh, name: name)], unit: unit);

  /// Writes every part as its own `<object>` in one 3MF.
  static Uint8List writeAll(
    List<ThreeMfPart> parts, {
    String name = 'pali_ijo',
    String unit = 'millimeter',
  }) {
    if (parts.isEmpty) {
      throw ArgumentError.value(parts, 'parts', 'nothing to write');
    }
    if (!_units.contains(unit)) {
      throw ArgumentError.value(
        unit,
        'unit',
        '3MF allows only ${_units.join(', ')}',
      );
    }

    final archive = ZipArchive()
      ..addText('[Content_Types].xml', _contentTypes)
      ..addText('_rels/.rels', _rels)
      ..addText(_modelPath, _model(parts, title: name, unit: unit));
    return archive.toBytes();
  }

  /// Reads the first object's mesh.
  ///
  /// Coordinates come back in millimetres regardless of the file's declared
  /// unit, since the kernel has one length unit and converting on the way in
  /// is the only place the conversion can be done once. Multi-object files
  /// should use [readAll].
  static MeshData read(Uint8List bytes) {
    final parts = readAll(bytes);
    if (parts.isEmpty) {
      throw const FormatException('3MF model has no triangles');
    }
    return parts.first.mesh;
  }

  /// Reads the build's placed instances, resolving component objects and
  /// external model parts. Uninstantiated resources are not extra bodies.
  static List<ThreeMfPart> readAll(Uint8List bytes) =>
      _ThreeMfReader(ZipArchive.fromBytes(bytes)).read();

  static const _contentTypes =
      '<?xml version="1.0" encoding="UTF-8"?>\n'
      '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
      '<Default Extension="rels" '
      'ContentType="application/vnd.openxmlformats-package.relationships+xml" />'
      '<Default Extension="model" '
      'ContentType="application/vnd.ms-package.3dmanufacturing-3dmodel+xml" />'
      '</Types>\n';

  static const _rels =
      '<?xml version="1.0" encoding="UTF-8"?>\n'
      '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
      '<Relationship Target="/$_modelPath" Id="rel0" Type="$_relationshipType" />'
      '</Relationships>\n';

  static String _model(
    List<ThreeMfPart> parts, {
    required String title,
    required String unit,
  }) {
    String f(double x) {
      final s = x.toStringAsFixed(6);
      return s == '-0.000000' ? '0.000000' : s;
    }

    final out = StringBuffer()
      ..writeln('<?xml version="1.0" encoding="UTF-8"?>')
      ..writeln(
        '<model unit="$unit" xml:lang="en-US" xmlns="$_modelNamespace">',
      )
      ..writeln('  <metadata name="Application">pali ijo</metadata>')
      ..writeln('  <metadata name="Title">${_escape(title)}</metadata>')
      ..writeln('  <resources>');

    for (var i = 0; i < parts.length; i++) {
      final part = parts[i];
      final id = i + 1;
      final label = part.name;
      out.write('    <object id="$id" type="model"');
      if (label != null && label.isNotEmpty) {
        out.write(' name="${_escape(label)}"');
      }
      out
        ..writeln('>')
        ..writeln('      <mesh>')
        ..writeln('        <vertices>');
      for (final v in part.mesh.vertices) {
        out.writeln(
          '          <vertex x="${f(v.x)}" y="${f(v.y)}" z="${f(v.z)}" />',
        );
      }
      out
        ..writeln('        </vertices>')
        ..writeln('        <triangles>');
      for (final face in part.mesh.faces) {
        if (face.length < 3) continue;
        for (var t = 1; t + 1 < face.length; t++) {
          out.writeln(
            '          <triangle v1="${face[0]}" v2="${face[t]}" '
            'v3="${face[t + 1]}" />',
          );
        }
      }
      out
        ..writeln('        </triangles>')
        ..writeln('      </mesh>')
        ..writeln('    </object>');
    }

    out.writeln('  </resources>');
    out.writeln('  <build>');
    for (var i = 0; i < parts.length; i++) {
      out.writeln('    <item objectid="${i + 1}" />');
    }
    out
      ..writeln('  </build>')
      ..writeln('</model>');
    return out.toString();
  }

  static String _escape(String s) => const HtmlEscape().convert(s);
}
