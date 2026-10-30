part of 'step_read_analytic.dart';

class _StepInstance {
  _StepInstance(this.child, this.transform, this.name);
  final int child;
  final BrepAffineTransform transform;
  final String name;
}

extension _StepAssembly on StepReader {
  (List<BrepSolid>, List<String>) _placeAssembly(
    List<int> ids,
    List<BrepSolid> definitions,
    List<String> names,
  ) {
    final relationships = _file
        .idsOf('REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION')
        .toList();
    if (relationships.isEmpty) return (definitions, names);
    final representatives = <int, int>{};
    int group(int id) {
      final parent = representatives[id];
      if (parent == null) return representatives[id] = id;
      return parent == id ? id : representatives[id] = group(parent);
    }

    // Geometry and product shape representations often describe the same
    // coordinate space in separate records (e.g. CATIA's assembly exports).
    for (final entry in _file.entities.entries) {
      if (_file.typed(
            entry.key,
            'REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION',
          ) !=
          null)
        continue;
      final relation = _file.anyOf(entry.key, [
        'SHAPE_REPRESENTATION_RELATIONSHIP',
        'REPRESENTATION_RELATIONSHIP',
      ]);
      final a = relation?.arg(2), b = relation?.arg(3);
      if (a is StepRef && b is StepRef)
        representatives[group(b.id)] = group(a.id);
    }
    final bodyIndex = {for (var i = 0; i < ids.length; i++) ids[i]: i};
    final bodiesAt = <int, Set<int>>{};
    for (final entry in _file.entities.entries) {
      for (final record in entry.value) {
        if (record.type != 'REPRESENTATION' &&
            !record.type.endsWith('SHAPE_REPRESENTATION'))
          continue;
        final items = record.arg(1);
        if (items is! StepList) continue;
        for (final item in items.items) {
          if (item is StepRef && bodyIndex.containsKey(item.id)) {
            (bodiesAt[group(entry.key)] ??= {}).add(item.id);
          }
        }
      }
    }
    final productShapes = <int, Set<int>>{};
    for (final id in _file.idsOf('SHAPE_DEFINITION_REPRESENTATION')) {
      final r = _file.typed(id, 'SHAPE_DEFINITION_REPRESENTATION')!;
      final definition = _file
          .follow(r.arg(0), 'PRODUCT_DEFINITION_SHAPE')
          ?.arg(2);
      final shape = r.arg(1);
      if (definition is StepRef && shape is StepRef) {
        (productShapes[definition.id] ??= {}).add(group(shape.id));
      }
    }
    final usages = <int, StepRecord>{};
    for (final id in _file.idsOf('CONTEXT_DEPENDENT_SHAPE_REPRESENTATION')) {
      final r = _file.typed(id, 'CONTEXT_DEPENDENT_SHAPE_REPRESENTATION')!;
      final relation = r.arg(0);
      final definition = _file.follow(r.arg(1), 'PRODUCT_DEFINITION_SHAPE');
      final usage = _file.follow(
        definition?.arg(2),
        'NEXT_ASSEMBLY_USAGE_OCCURRENCE',
      );
      if (relation is StepRef && usage != null) usages[relation.id] = usage;
    }
    final children = <int, List<_StepInstance>>{};
    final childGroups = <int>{};
    for (final id in relationships) {
      final derived = _file.typed(
        id,
        'REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION',
      )!;
      final base = _file.typed(id, 'REPRESENTATION_RELATIONSHIP') ?? derived;
      final r1 = base.arg(2), r2 = base.arg(3);
      final transformRef = derived.arg(derived.args.length == 1 ? 0 : 4);
      final transform = _file.follow(
        transformRef,
        'ITEM_DEFINED_TRANSFORMATION',
      );
      if (r1 is! StepRef || r2 is! StepRef || transform == null) {
        throw StepBrepUnsupportedException(
          'Unreadable STEP assembly placement #$id',
        );
      }
      final from = _assemblyFrame(transform.arg(2));
      final to = _assemblyFrame(transform.arg(3));
      var parent = group(r2.id), child = group(r1.id);
      var matrix = BrepAffineTransform.fromFunction(
        (p) => to.point(from.inversePoint(p)),
      );
      final usage = usages[id];
      final parentProduct = usage?.arg(3);
      if (parentProduct is StepRef &&
          productShapes[parentProduct.id]?.contains(group(r1.id)) == true &&
          productShapes[parentProduct.id]?.contains(group(r2.id)) != true) {
        parent = group(r1.id);
        child = group(r2.id);
        matrix = BrepAffineTransform.fromFunction(
          (p) => from.point(to.inversePoint(p)),
        );
      }
      final nameValue = usage?.arg(1);
      final name = nameValue is StepText ? nameValue.value.trim() : '';
      (children[parent] ??= []).add(_StepInstance(child, matrix, name));
      childGroups.add(child);
    }
    final roots = children.keys
        .where((id) => !childGroups.contains(id))
        .toList();
    if (roots.isEmpty) throw const FormatException('Cyclic STEP assembly');
    final solids = <BrepSolid>[], placedNames = <String>[];
    final emitted = <int>{}, active = <int>{}, visited = <int>{};
    final identity = BrepAffineTransform(
      Vector3(1, 0, 0),
      Vector3(0, 1, 0),
      Vector3(0, 0, 1),
      Vector3.zero(),
    );
    void visit(int rep, BrepAffineTransform placement, String label) {
      if (!active.add(rep)) throw const FormatException('Cyclic STEP assembly');
      visited.add(rep);
      for (final id in bodiesAt[rep] ?? const <int>{}) {
        final index = bodyIndex[id]!;
        solids.add(BrepXform.rigid(definitions[index], placement));
        placedNames.add(label.isEmpty ? names[index] : label);
        emitted.add(id);
      }
      for (final instance in children[rep] ?? const <_StepInstance>[]) {
        visit(
          instance.child,
          instance.transform.then(placement),
          instance.name.isEmpty ? label : instance.name,
        );
      }
      active.remove(rep);
    }

    for (final root in roots) visit(root, identity, '');
    if (children.keys.any((id) => !visited.contains(id))) {
      throw const FormatException('Cyclic STEP assembly');
    }
    // Preserve standalone bodies/surface models that have no assembly usage.
    for (var i = 0; i < ids.length; i++) {
      if (!emitted.contains(ids[i])) {
        solids.add(definitions[i]);
        placedNames.add(names[i]);
      }
    }
    _notes.add(
      'Placed ${solids.length} CAD bodies using STEP assembly instances.',
    );
    return (solids, placedNames);
  }

  BrepAffineTransform _assemblyFrame(StepValue? ref) {
    final frame = _placement(ref);
    if (frame == null)
      throw const FormatException('Missing STEP assembly coordinate frame');
    final z = frame.axis;
    final x = (frame.refDir - z * frame.refDir.dot(z)).normalized();
    return BrepAffineTransform(x, z.cross(x).normalized(), z, frame.origin);
  }
}
