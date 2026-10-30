part of 'threemf.dart';

const _productionNamespace =
    'http://schemas.microsoft.com/3dmanufacturing/production/2015/06';

class _ThreeMfObject {
  _ThreeMfObject(this.id, this.name);
  final int id;
  final String? name;
  final vertices = <Vector3>[];
  final faces = <List<int>>[];
  final components = <_ThreeMfReference>[];
}

class _ThreeMfReference {
  _ThreeMfReference(this.id, this.path, this.transform);
  final int id;
  final String? path;
  final Matrix4 transform;
}

class _ThreeMfModel {
  final objects = <int, _ThreeMfObject>{};
  final build = <_ThreeMfReference>[];
  double scale = 1;
  bool hasBuild = false;
}

class _ThreeMfReader {
  _ThreeMfReader(this.archive);
  final ZipArchive archive;
  final models = <String, _ThreeMfModel>{};

  List<ThreeMfPart> read() {
    String? rootPath;
    final relationships = archive.readText('_rels/.rels');
    if (relationships != null) {
      for (final e in parseEvents(relationships, validateNesting: true)) {
        if (e is! XmlStartElementEvent || e.localName != 'Relationship') {
          continue;
        }
        final attrs = {for (final a in e.attributes) a.localName: a.value};
        if (attrs['Type'] == IjoThreeMf._relationshipType) {
          if (attrs['TargetMode'] == 'External') {
            throw const FormatException('3MF model must be inside the package');
          }
          rootPath = _path('', attrs['Target'] ?? '');
          break;
        }
      }
    }
    rootPath ??= archive.names.contains(IjoThreeMf._modelPath)
        ? IjoThreeMf._modelPath
        : archive.names
              .where((n) => n.toLowerCase().endsWith('.model'))
              .firstOrNull;
    if (rootPath == null) {
      throw const FormatException('Not a 3MF: no model part');
    }
    final root = _model(rootPath);
    // Tolerate older mesh-only producers without a build section. An explicit
    // empty build is not permission to import unused library resources.
    final build = root.hasBuild
        ? root.build
        : [
            for (final o in root.objects.values)
              _ThreeMfReference(o.id, null, Matrix4.identity()),
          ];
    final result = <ThreeMfPart>[];
    for (final item in build) {
      final vertices = <Vector3>[];
      final faces = <List<int>>[];
      final path = item.path == null ? rootPath : _path(rootPath, item.path!);
      final label = _model(path).objects[item.id]?.name;
      _append(path, item.id, item.transform, {}, vertices, faces);
      if (faces.isEmpty) {
        throw const FormatException('3MF build item has no triangles');
      }
      result.add(
        ThreeMfPart(
          MeshData(vertices: vertices, faces: faces),
          name: label,
        ),
      );
    }
    if (result.isEmpty) {
      throw const FormatException('3MF model has no triangles');
    }
    return result;
  }

  void _append(
    String path,
    int id,
    Matrix4 transform,
    Set<(String, int)> active,
    List<Vector3> vertices,
    List<List<int>> faces,
  ) {
    final key = (path, id);
    if (!active.add(key)) {
      throw FormatException('Cyclic 3MF component $path#$id');
    }
    final model = _model(path), object = model.objects[id];
    if (object == null) throw FormatException('Missing 3MF object $path#$id');
    final base = vertices.length;
    for (final p in object.vertices) {
      final v = transform.transformed3(p);
      if (!v.x.isFinite || !v.y.isFinite || !v.z.isFinite) {
        throw const FormatException('Non-finite 3MF transformed coordinate');
      }
      vertices.add(v);
    }
    final mirrored = transform.determinant() < 0;
    for (final f in object.faces) {
      faces.add([
        base + f[0],
        base + f[mirrored ? 2 : 1],
        base + f[mirrored ? 1 : 2],
      ]);
    }
    for (final child in object.components) {
      final childPath = child.path == null ? path : _path(path, child.path!);
      _append(
        childPath,
        child.id,
        transform.multiplied(child.transform),
        active,
        vertices,
        faces,
      );
    }
    active.remove(key);
  }

  _ThreeMfModel _model(String path) => models.putIfAbsent(path, () {
    final xml = archive.readText(path);
    if (xml == null) throw FormatException('Missing 3MF model part "$path"');
    final model = _ThreeMfModel();
    _ThreeMfObject? object;
    var seenModel = false;
    // Pull parsing avoids retaining a DOM with millions of vertex/triangle
    // nodes. Only numeric geometry and component references survive parsing.
    for (final event in parseEvents(
      xml,
      validateNesting: true,
      validateDocument: true,
      withParent: true,
    )) {
      if (event is XmlEndElementEvent && event.localName == 'object') {
        object = null;
        continue;
      }
      if (event is! XmlStartElementEvent) continue;
      if (event.localName.toLowerCase() == 'beamlattice') {
        throw UnsupportedError('3MF beam lattice geometry is not supported');
      }
      if (event.namespaceUri != null &&
          event.namespaceUri != IjoThreeMf._modelNamespace) {
        continue;
      }
      final attrs = {
        for (final a in event.attributes)
          if (a.namespacePrefix == null) a.localName: a.value,
      };
      switch (event.localName) {
        case 'model':
          seenModel = true;
          model.scale = switch (attrs['unit'] ?? 'millimeter') {
            'micron' => .001,
            'millimeter' => 1,
            'centimeter' => 10,
            'inch' => 25.4,
            'foot' => 304.8,
            'meter' => 1000,
            _ => throw FormatException('Unknown 3MF unit "${attrs['unit']}"'),
          };
        case 'object':
          final id = _integer(attrs, 'id');
          if (model.objects.containsKey(id)) {
            throw FormatException('Duplicate 3MF object $id');
          }
          object = _ThreeMfObject(id, attrs['name']);
          model.objects[id] = object;
        case 'mesh':
          if (object == null) {
            object = _ThreeMfObject(0, null);
            model.objects[0] = object;
          }
        case 'vertex':
          if (object == null) {
            throw const FormatException('3MF vertex outside object');
          }
          object.vertices.add(
            Vector3(
              _number(attrs, 'x') * model.scale,
              _number(attrs, 'y') * model.scale,
              _number(attrs, 'z') * model.scale,
            ),
          );
        case 'triangle':
          if (object == null) {
            throw const FormatException('3MF triangle outside object');
          }
          final indices = [
            for (final k in ['v1', 'v2', 'v3']) _integer(attrs, k),
          ];
          if (indices.any((i) => i < 0 || i >= object!.vertices.length)) {
            throw const FormatException(
              '3MF triangle references an undefined vertex',
            );
          }
          object.faces.add(indices);
        case 'build':
          model.hasBuild = true;
        case 'component':
        case 'item':
          String? componentPath;
          for (final a in event.attributes) {
            if (a.localName == 'path' &&
                a.namespaceUri == _productionNamespace) {
              componentPath = a.value;
            }
          }
          final ref = _ThreeMfReference(
            _integer(attrs, 'objectid'),
            componentPath,
            _transform(attrs['transform'], model.scale),
          );
          if (event.localName == 'component') {
            if (object == null) {
              throw const FormatException('3MF component outside object');
            }
            object.components.add(ref);
          } else {
            model.build.add(ref);
          }
      }
      if (event.isSelfClosing && event.localName == 'object') object = null;
    }
    if (!seenModel) throw FormatException('Not a 3MF model: "$path"');
    for (final object in model.objects.values) {
      if (object.components.isNotEmpty && object.faces.isNotEmpty) {
        throw const FormatException(
          '3MF object has both a mesh and components',
        );
      }
    }
    return model;
  });

  static int _integer(Map<String, String> attrs, String key) {
    final value = int.tryParse(attrs[key] ?? '');
    if (value == null) throw FormatException('Invalid 3MF integer "$key"');
    return value;
  }

  static double _number(Map<String, String> attrs, String key) {
    final value = double.tryParse(attrs[key] ?? '');
    if (value == null || !value.isFinite) {
      throw FormatException('Invalid 3MF coordinate "$key"');
    }
    return value;
  }

  static Matrix4 _transform(String? value, double scale) {
    if (value == null) return Matrix4.identity();
    final values = value
        .trim()
        .split(RegExp(r'\s+'))
        .map(double.tryParse)
        .toList();
    if (values.length != 12 || values.any((n) => n == null || !n.isFinite)) {
      throw const FormatException('Invalid 3MF transform');
    }
    final m = values.cast<double>();
    // 3MF row-vector affine layout transposed into vector_math's column vectors.
    return Matrix4(
      m[0],
      m[1],
      m[2],
      0,
      m[3],
      m[4],
      m[5],
      0,
      m[6],
      m[7],
      m[8],
      0,
      m[9] * scale,
      m[10] * scale,
      m[11] * scale,
      1,
    );
  }

  static String _path(String current, String reference) {
    final uri = Uri.parse(reference);
    if (uri.hasScheme ||
        uri.hasAuthority ||
        uri.hasQuery ||
        uri.hasFragment ||
        reference.isEmpty) {
      throw FormatException('Invalid 3MF package reference "$reference"');
    }
    final base = Uri(path: '/$current');
    return Uri.decodeComponent(
      base.resolveUri(uri).path,
    ).replaceFirst(RegExp(r'^/'), '');
  }
}
