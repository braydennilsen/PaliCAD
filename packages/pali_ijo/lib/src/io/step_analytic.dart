/// STEP BREP interchange with native analytic and rational spline geometry.
/// Conics and affine quadrics have exact rational representations. More general
/// bounded curves and ruled walls use adaptive spline fitting for interchange;
/// the original BREP is unchanged. Unsupported entities fail explicitly.
library;

import 'dart:math' as math;
import 'dart:typed_data';

import 'package:vector_math/vector_math_64.dart';

import '../brep/geom.dart';
import '../brep/solid.dart';
import '../brep/spline_convert.dart';
import '../brep/tessellate.dart';
import '../brep/heal.dart';
import '../brep/boolean.dart';
import '../measure/properties.dart';
import 'step_read_analytic.dart';
import 'step_syntax.dart';

export 'step_read_analytic.dart'
    show StepBrepUnsupportedException, StepBrepReadResult;
export 'step_syntax.dart' show StepSyntaxException;

/// One solid on its way into a STEP file.
class StepBrepPart {
  const StepBrepPart(this.solid, {this.name});

  final BrepSolid solid;

  /// The name the receiving system will show. Defaults to `Part`, numbered
  /// when a file holds several.
  final String? name;
}

/// Reads and writes STEP.
class IjoStepBrep {
  IjoStepBrep._();

  /// The schema written into `FILE_SCHEMA`.
  ///
  /// AP214 rather than AP203 because it is a superset for these entities and
  /// is what most systems expect to see; the entities used are common to both.
  static const schema = 'AUTOMOTIVE_DESIGN { 1 0 10303 214 1 1 1 1 }';

  /// Millimetres, matching the kernel's own unit.
  static const _lengthUnit = 'SI_UNIT(.MILLI.,.METRE.)';

  /// Writes one solid.
  static Uint8List write(
    BrepSolid solid, {
    String name = 'Part',
    DateTime? timestamp,
  }) => writeAll([StepBrepPart(solid, name: name)], timestamp: timestamp);

  /// Writes several solids into one file.
  ///
  /// They share one shape representation, each as its own
  /// `MANIFOLD_SOLID_BREP` carrying its own name. That is the form every
  /// reader handles; a true assembly structure, with occurrences and their own
  /// placements, is a different feature and this kernel has no scene graph to
  /// build one from — placement belongs to the host.
  static Uint8List writeAll(
    List<StepBrepPart> parts, {
    String fileName = 'part.step',
    DateTime? timestamp,
  }) {
    if (parts.isEmpty) {
      throw ArgumentError.value(parts, 'parts', 'nothing to write');
    }
    for (final part in parts) {
      if (part.solid.faces.isEmpty) {
        throw ArgumentError(
          'Solid "${part.name ?? 'Part'}" has no faces, so it is not a solid '
          'STEP can describe',
        );
      }
    }

    final w = StepWriter();
    final g = _Geometry(w);

    // The product and context scaffolding. Verbose and almost entirely fixed:
    // AP203 requires a product with a definition before a shape may be
    // attached to it, so a bare body is not expressible.
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
      _lengthUnit,
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
    // The tolerance a receiving system should assume. 1e-7 mm is well below
    // anything the kernel's own 1e-6 comparisons treat as distinct.
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
      final part = parts[i],
          name = part.name ?? (parts.length == 1 ? 'Part' : 'Part ${i + 1}');
      if (!BrepHeal.topologyIsClosed(part.solid))
        throw BrepUnsupported('STEP export requires a closed BREP');
      final solids = brepSplitSolids(part.solid);
      for (var j = 0; j < solids.length; j++)
        bodies.add(
          g.solid(solids[j], solids.length == 1 ? name : '$name ${j + 1}'),
        );
    }

    final origin = g.placement(
      Vector3.zero(),
      Vector3(0, 0, 1),
      Vector3(1, 0, 0),
    );
    final shape = w.add(
      'ADVANCED_BREP_SHAPE_REPRESENTATION',
      '${StepWriter.text('')},'
          '${StepWriter.refs([origin, ...bodies])},#$context',
    );

    // One product per file, named after the first part.
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
    final productName = parts.length == 1
        ? (parts.first.name ?? 'Part')
        : _stem(fileName);
    final product = w.add(
      'PRODUCT',
      '${StepWriter.text(productName)},${StepWriter.text(productName)},'
          '${StepWriter.text('')},(#$productContext)',
    );
    final formation = w.add(
      'PRODUCT_DEFINITION_FORMATION_WITH_SPECIFIED_SOURCE',
      '${StepWriter.text('')},${StepWriter.text('')},#$product,'
          '.NOT_KNOWN.',
    );
    final definition = w.add(
      'PRODUCT_DEFINITION',
      '${StepWriter.text('design')},${StepWriter.text('')},'
          '#$formation,#$definitionContext',
    );
    final definitionShape = w.add(
      'PRODUCT_DEFINITION_SHAPE',
      '${StepWriter.text('')},${StepWriter.text('')},#$definition',
    );
    w.add('SHAPE_DEFINITION_REPRESENTATION', '#$definitionShape,#$shape');

    return Uint8List.fromList(
      w.bytes(
        description:
            'pali ijo BREP; spline approximation target 1E-7 mm where needed',
        fileName: fileName,
        originator: 'pali ijo',
        schema: schema,
        timestamp: timestamp,
      ),
    );
  }

  /// Reads the first solid in [bytes].
  ///
  /// Throws [StepBrepUnsupportedException] when the file's geometry is outside
  /// what this kernel represents, and [StepSyntaxException] when it is not
  /// STEP at all.
  static BrepSolid read(Uint8List bytes) => readAll(bytes).solids.first;

  /// Reads every solid in [bytes], with the units and any refusals reported.
  static StepBrepReadResult readAll(
    Uint8List bytes, {
    bool allowOpenShells = false,
    bool repairDegenerateFaces = false,
    void Function(int completed, int total)? onProgress,
  }) => StepReader.read(
    bytes,
    allowOpenShells: allowOpenShells,
    repairDegenerateFaces: repairDegenerateFaces,
    onProgress: onProgress,
  );

  static String _stem(String path) {
    final slash = path.lastIndexOf(RegExp(r'[/\\]'));
    final base = slash < 0 ? path : path.substring(slash + 1);
    final dot = base.lastIndexOf('.');
    return dot <= 0 ? base : base.substring(0, dot);
  }
}

/// Emits the geometry and topology of solids into a [StepWriter].
///
/// Held as a class so points, directions and placements can be shared across
/// the whole file rather than per solid.
class _Geometry {
  _Geometry(this.w);

  final StepWriter w;

  /// `EDGE_CURVE` ids for the solid being written, by its edge index.
  var _edgeIds = <int, int>{};

  /// Points already written, bucketed by position to a micron.
  final _points = <String, List<(Vector3, int)>>{};

  /// Writes a point, reusing one that is already there.
  ///
  /// Coincidence is judged by distance, not by the text the numbers format to.
  /// Two ways of arriving at the same corner rarely agree to the last bit — a
  /// rim's `evaluate(0)` against the seam line's stored endpoint, say — and
  /// comparing text would write them as two points a femtometre apart. In a
  /// receiving system that is not a rounding difference, it is a crack: the
  /// two edges no longer share a vertex, so the shell is no longer closed and
  /// the import is rejected or silently repaired.
  int point(Vector3 p) {
    final cx = (p.x * 1000).round();
    final cy = (p.y * 1000).round();
    final cz = (p.z * 1000).round();
    for (var dx = -1; dx <= 1; dx++) {
      for (var dy = -1; dy <= 1; dy++) {
        for (var dz = -1; dz <= 1; dz++) {
          final bucket = _points['${cx + dx},${cy + dy},${cz + dz}'];
          if (bucket == null) continue;
          for (final (at, id) in bucket) {
            if (at.distanceToSquared(p) < 1e-18) return id;
          }
        }
      }
    }
    final id = w.add(
      'CARTESIAN_POINT',
      '${StepWriter.text('')},'
          '(${StepWriter.num(p.x)},${StepWriter.num(p.y)},'
          '${StepWriter.num(p.z)})',
    );
    (_points['$cx,$cy,$cz'] ??= []).add((p, id));
    return id;
  }

  int direction(Vector3 d) {
    final n = d.normalized();
    return w.dedupe(
      'DIRECTION',
      '${StepWriter.text('')},'
          '(${StepWriter.num(n.x)},${StepWriter.num(n.y)},'
          '${StepWriter.num(n.z)})',
    );
  }

  /// `AXIS2_PLACEMENT_3D`: an origin, a local Z, and a local X.
  ///
  /// Every analytic surface and curve in STEP is positioned by one of these,
  /// which is convenient, because it is also how this kernel stores them.
  int placement(Vector3 origin, Vector3 axis, Vector3 refDir) {
    // Orthogonalise, since a placement whose axes are not perpendicular is
    // interpreted by projecting anyway, and being explicit avoids differences
    // between readers.
    final z = axis.normalized();
    var x = refDir - z * refDir.dot(z);
    if (x.length < 1e-12) x = brepFrameU(z);
    return w.dedupe(
      'AXIS2_PLACEMENT_3D',
      '${StepWriter.text('')},#${point(origin)},#${direction(z)},'
          '#${direction(x.normalized())}',
    );
  }

  int vertex(Vector3 p) =>
      w.dedupe('VERTEX_POINT', '${StepWriter.text('')},#${point(p)}');

  /// Writes one solid and returns its `MANIFOLD_SOLID_BREP` id.
  int solid(BrepSolid s, String name) {
    _edgeIds = {};

    for (var i = 0; i < s.edges.length; i++) {
      final c = s.edges[i].curve, a = c.evaluate(c.t0);
      if ((c.evaluate(c.t1) - a).length < 1e-9 &&
          [0.25, .5, .75].every(
            (t) => (c.evaluate(c.t0 + (c.t1 - c.t0) * t) - a).length < 1e-9,
          ))
        continue;
      _edgeIds[i] = _edgeCurve(s, i);
    }

    final adjacency = <int, List<int>>{};
    for (var i = 0; i < s.faces.length; i++) {
      for (final l in s.faces[i].loops) {
        for (final co in l.coedges) {
          (adjacency[co.edge] ??= []).add(i);
        }
      }
    }
    final seen = <int>{}, groups = <List<int>>[];
    for (var i = 0; i < s.faces.length; i++) {
      if (!seen.add(i)) continue;
      final stack = [i], group = <int>[];
      while (stack.isNotEmpty) {
        final f = stack.removeLast();
        group.add(f);
        for (final l in s.faces[f].loops) {
          for (final co in l.coedges) {
            for (final n in adjacency[co.edge]!) {
              if (seen.add(n)) stack.add(n);
            }
          }
        }
      }
      groups.add(group);
    }
    int? outer;
    final voids = <int>[];
    for (final group in groups) {
      final sample = BrepSolid(
        vertices: s.vertices,
        edges: s.edges,
        faces: [for (final i in group) s.faces[i]],
      );
      final inward = groups.length > 1 && IjoMeasure.volumeSolid(sample) < 0;
      final faces = [
        for (final i in group)
          _face(
            s,
            inward
                ? BrepFace(
                    surface: s.faces[i].surface,
                    loops: [
                      for (final loop in s.faces[i].loops)
                        BrepLoop([
                          for (final c in loop.coedges.reversed)
                            BrepCoedge(c.edge, reversed: !c.reversed),
                        ]),
                    ],
                    reversed: !s.faces[i].reversed,
                  )
                : s.faces[i],
          ),
      ];
      final shell = w.add('CLOSED_SHELL', "'',${StepWriter.refs(faces)}");
      if (inward) {
        voids.add(w.add('ORIENTED_CLOSED_SHELL', "'',*,#$shell,.F."));
      } else {
        if (outer != null)
          throw BrepUnsupported('STEP part has multiple exterior shells');
        outer = shell;
      }
    }
    if (outer == null) throw BrepUnsupported('STEP part has no exterior shell');
    return voids.isEmpty
        ? w.add('MANIFOLD_SOLID_BREP', '${StepWriter.text(name)},#$outer')
        : w.add(
            'BREP_WITH_VOIDS',
            '${StepWriter.text(name)},#$outer,${StepWriter.refs(voids)}',
          );
  }

  int _edgeCurve(BrepSolid s, int index) {
    final edge = s.edges[index];
    final curve = _curve(edge.curve);

    // The vertices are taken from where the curve actually starts and ends,
    // not from the solid's vertex list.
    //
    // The two disagree more often than they should: a fillet can leave an
    // edge whose ends are recorded as one vertex while its arc sweeps 90°,
    // which is harmless inside the kernel — operations read the range off the
    // curve — and unreadable once written, because in STEP the vertices are
    // the *only* record of which part of a circle an edge covers. Writing the
    // geometry's own endpoints makes the file self-consistent. Positions are
    // deduplicated, so edges meeting at a corner still share one vertex and
    // the topology stays connected.
    final start = edge.curve.evaluate(edge.curve.t0);
    final end = edge.curve.evaluate(edge.curve.t1);

    // An edge that comes back to where it started gets one vertex, not two
    // that happen to coincide. A full turn does not evaluate to bit-identical
    // ends — `sin(2π)` is not zero in floating point — and left alone that
    // writes a circular rim as an edge with a gap of 1e-16 in it, which is
    // both untrue and the kind of thing a strict reader refuses.
    final a = vertex(start);
    final b = start.distanceToSquared(end) < 1e-18 ? a : vertex(end);

    // `same_sense` carries the sweep direction. A STEP circle always runs
    // anticlockwise about its placement axis, while this kernel says
    // clockwise with a descending range — so a descending range is exactly a
    // curve the edge runs against.
    final ascending =
        edge.curve is BrepAffineCurve ||
        edge.curve is BrepPlaneSectionCurve ||
        edge.curve is BrepSurfaceIntersectionCurve ||
        edge.curve is BrepRuledQuadricCurve ||
        edge.curve is BrepIntersectionCurve ||
        edge.curve.t1 >= edge.curve.t0;

    return w.add(
      'EDGE_CURVE',
      '${StepWriter.text('')},#$a,#$b,#$curve,'
          '${ascending ? '.T.' : '.F.'}',
    );
  }

  (List<int>, List<double>) _knots(List<double> expanded) {
    final mult = <int>[], values = <double>[];
    for (final x in expanded) {
      if (values.isNotEmpty && values.last == x) {
        mult[mult.length - 1]++;
      } else {
        values.add(x);
        mult.add(1);
      }
    }
    return (mult, values);
  }

  String _numbers(List<double> xs) => '(${xs.map(StepWriter.num).join(',')})';
  int _splineCurve(BrepSplineCurve c) {
    final k = _knots(c.knots);
    return w.addComplex([
      'BOUNDED_CURVE()',
      'B_SPLINE_CURVE(${c.degree},${StepWriter.refs(c.points.map(point).toList())},.UNSPECIFIED.,${(c.points.first - c.points.last).length < 1e-9 ? '.T.' : '.F.'},.U.)',
      'B_SPLINE_CURVE_WITH_KNOTS((${k.$1.join(',')}),${_numbers(k.$2)},.UNSPECIFIED.)',
      'CURVE()',
      'GEOMETRIC_REPRESENTATION_ITEM()',
      'RATIONAL_B_SPLINE_CURVE(${_numbers(c.weights)})',
      "REPRESENTATION_ITEM('')",
    ]);
  }

  int _splineSurface(BrepSplineSurface s) {
    final u = _knots(s.knotsU), v = _knots(s.knotsV);
    final points =
        '(${s.points.map((r) => StepWriter.refs(r.map(point).toList())).join(',')})';
    final weights = '(${s.weights.map(_numbers).join(',')})';
    return w.addComplex([
      'BOUNDED_SURFACE()',
      'B_SPLINE_SURFACE(${s.degreeU},${s.degreeV},$points,.UNSPECIFIED.,${List.generate(9, (i) => i / 8).every((t) => (s.evaluate(s.u0, s.v0 + (s.v1 - s.v0) * t) - s.evaluate(s.u1, s.v0 + (s.v1 - s.v0) * t)).length < 1e-8) ? '.T.' : '.F.'},${List.generate(9, (i) => i / 8).every((t) => (s.evaluate(s.u0 + (s.u1 - s.u0) * t, s.v0) - s.evaluate(s.u0 + (s.u1 - s.u0) * t, s.v1)).length < 1e-8) ? '.T.' : '.F.'},.U.)',
      'B_SPLINE_SURFACE_WITH_KNOTS((${u.$1.join(',')}),(${v.$1.join(',')}),${_numbers(u.$2)},${_numbers(v.$2)},.UNSPECIFIED.)',
      'GEOMETRIC_REPRESENTATION_ITEM()',
      'RATIONAL_B_SPLINE_SURFACE($weights)',
      "REPRESENTATION_ITEM('')",
      'SURFACE()',
    ]);
  }

  int _curve(BrepCurve curve) {
    switch (curve) {
      case BrepHelix():
        return _splineCurve(BrepSplineConvert.curve(curve));
      case BrepSplineCurve():
        return _splineCurve(curve);
      case BrepAffineCurve():
        return _splineCurve(BrepSplineConvert.curve(curve));
      case BrepIntersectionCurve():
        final representation = _splineCurve(BrepSplineConvert.curve(curve));
        return w.add(
          'INTERSECTION_CURVE',
          "'',#$representation,(#${_surface(curve.carrier)},#${_surface(curve.target)}),.CURVE_3D.",
        );
      case BrepSurfaceIntersectionCurve():
        final representation = _splineCurve(BrepSplineConvert.curve(curve));
        return w.add(
          'INTERSECTION_CURVE',
          "'',#$representation,(#${_surface(curve.first)},#${_surface(curve.second)}),.CURVE_3D.",
        );
      case BrepRuledQuadricCurve():
      case BrepPlaneSectionCurve():
        return _splineCurve(BrepSplineConvert.curve(curve));
      case BrepLine line:
        final d = line.b - line.a;
        final length = d.length;
        if (length < 1e-12) {
          throw ArgumentError('A zero-length line cannot be written to STEP');
        }
        final vector = w.dedupe(
          'VECTOR',
          '${StepWriter.text('')},#${direction(d)},'
              '${StepWriter.num(length)}',
        );
        return w.dedupe(
          'LINE',
          '${StepWriter.text('')},#${point(line.a)},#$vector',
        );
      case BrepCircle circle:
        return w.dedupe(
          'CIRCLE',
          '${StepWriter.text('')},'
              '#${placement(circle.center, circle.normal, circle.xDir)},'
              '${StepWriter.num(circle.radius)}',
        );
      case BrepEllipse ellipse:
        return w.dedupe(
          'ELLIPSE',
          '${StepWriter.text('')},'
              '#${placement(ellipse.center, ellipse.normal, ellipse.xDir)},'
              '${StepWriter.num(ellipse.rMajor)},'
              '${StepWriter.num(ellipse.rMinor)}',
        );
    }
  }

  int _surface(BrepSurface surface) {
    switch (surface) {
      case BrepOffsetSurface():
        final base = surface.base is BrepRuledSurface
            ? BrepSplineConvert.surface(surface.base)
            : surface.base;
        return w.add(
          'OFFSET_SURFACE',
          "'',#${_surface(base)},${StepWriter.num(surface.distance)},.U.",
        );
      case BrepSplineSurface():
        return _splineSurface(surface);
      case BrepAffineSurface():
        throw BrepUnsupported(
          'STEP export of affine surfaces is not implemented',
        );
      case BrepRuledSurface():
        throw BrepUnsupported(
          'STEP export of a ruled loft surface is not implemented',
        );
      case BrepPlane plane:
        return w.dedupe(
          'PLANE',
          '${StepWriter.text('')},'
              '#${placement(plane.origin, plane.normal, plane.xDir)}',
        );
      case BrepCylinder cylinder:
        return w.dedupe(
          'CYLINDRICAL_SURFACE',
          '${StepWriter.text('')},'
              '#${placement(cylinder.origin, cylinder.axis, cylinder.xDir)},'
              '${StepWriter.num(cylinder.radius)}',
        );
      case BrepSphere sphere:
        return w.dedupe(
          'SPHERICAL_SURFACE',
          '${StepWriter.text('')},'
              '#${placement(sphere.center, Vector3(0, 0, 1), Vector3(1, 0, 0))},'
              '${StepWriter.num(sphere.radius)}',
        );
      case BrepTorus torus:
        return w.dedupe(
          'TOROIDAL_SURFACE',
          '${StepWriter.text('')},'
              '#${placement(torus.center, torus.axis, torus.xDir)},'
              '${StepWriter.num(torus.majorRadius)},'
              '${StepWriter.num(torus.minorRadius)}',
        );
      case BrepCone cone:
        return _cone(cone);
    }
  }

  /// A cone, which needs more care than the others.
  ///
  /// STEP's `CONICAL_SURFACE` is an infinite cone stated as a radius at the
  /// placement plane plus a half-angle, and the schema requires that angle to
  /// be strictly between 0 and 90°. This kernel's cone is a frustum given by
  /// two radii and a height, which can express two things STEP's cannot: a
  /// cone that narrows along its axis (a negative angle) and one that does not
  /// taper at all (a zero angle). Both are handled here rather than left to a
  /// reader to reject.
  int _cone(BrepCone cone) {
    final dr = cone.topRadius - cone.radius;
    if (dr.abs() < 1e-9 || cone.height.abs() < 1e-12) {
      // No taper: it is a cylinder, whatever it is called.
      return w.dedupe(
        'CYLINDRICAL_SURFACE',
        '${StepWriter.text('')},'
            '#${placement(cone.origin, cone.axis, cone.xDir)},'
            '${StepWriter.num(cone.radius)}',
      );
    }

    // Point the axis the way the cone widens, so the half-angle is positive.
    final widensUp = dr > 0;
    final apexward = widensUp ? cone.axis : -cone.axis;
    final base = widensUp ? cone.origin : cone.origin + cone.axis * cone.height;
    final baseRadius = widensUp ? cone.radius : cone.topRadius;
    final semiAngle = math.atan(dr.abs() / cone.height.abs());

    if (semiAngle >= math.pi / 2 - 1e-9) {
      throw ArgumentError(
        'A cone with a half-angle of ${(semiAngle * 180 / math.pi).toStringAsFixed(1)}° is flat, and STEP has no way to say that',
      );
    }

    return w.dedupe(
      'CONICAL_SURFACE',
      '${StepWriter.text('')},'
          '#${placement(base, apexward, cone.xDir)},'
          '${StepWriter.num(baseRadius)},${StepWriter.num(semiAngle)}',
    );
  }

  int _face(BrepSolid s, BrepFace face) {
    final original = face.surface;
    final converted =
        original is BrepAffineSurface && original.base is BrepPlane
        ? BrepPlane(
            origin: original.transform.point(
              (original.base as BrepPlane).origin,
            ),
            normal: original.transform.normal(
              (original.base as BrepPlane).normal,
            ),
            xDir: original.transform.vector((original.base as BrepPlane).xDir),
          )
        : original is BrepAffineSurface || original is BrepRuledSurface
        ? BrepSplineConvert.surface(
            original,
            boundary: [
              for (final loop in face.loops)
                ...BrepTessellate.sampleLoop(s, loop, 64),
            ],
          )
        : original;
    var reversed = face.reversed;
    if (!identical(original, converted) && converted is BrepSplineSurface) {
      final u = converted.u0 + (converted.u1 - converted.u0) * .37,
          v = converted.v0 + (converted.v1 - converted.v0) * .53,
          p = converted.evaluate(u, v),
          uv = original.parameters(p)!;
      if (converted.normalAt(u, v).dot(original.normalAt(uv.$1, uv.$2)) < 0)
        reversed = !reversed;
    }
    final surface = _surface(converted);
    final loops = face.loops.isEmpty
        ? [_seamLoop(converted)]
        : _bounds(s, face);
    return w.add(
      'ADVANCED_FACE',
      "'',${StepWriter.refs(loops)},#$surface,${reversed ? '.F.' : '.T.'}",
    );
  }

  /// The face's bounds, outer first.
  ///
  /// STEP distinguishes the outer bound from the holes inside it, and the
  /// kernel does not — a face just has loops. For a planar face the outer one
  /// is the one enclosing the largest area, which is cheap to work out and
  /// right. For a curved face the question barely arises: a cylinder wall or a
  /// blend has one loop, and where a hole exists the first loop is the outer
  /// one by construction.
  List<int> _bounds(BrepSolid s, BrepFace face) {
    final order = List<int>.generate(face.loops.length, (i) => i);
    if (face.loops.length > 1 && face.surface is BrepPlane) {
      final plane = face.surface as BrepPlane;
      final areas = [for (final loop in face.loops) _loopArea(s, loop, plane)];
      order.sort((a, b) => areas[b].compareTo(areas[a]));
    }

    // Corrected seam edges, shared across this face's loops. A seam is used
    // twice, once each way, and both uses must name the same edge.
    final seams = <int, int>{};

    final out = <int>[];
    for (var i = 0; i < order.length; i++) {
      final loop = _edgeLoop(s, face, face.loops[order[i]], seams);
      out.add(
        w.add(
          i == 0 ? 'FACE_OUTER_BOUND' : 'FACE_BOUND',
          '${StepWriter.text('')},#$loop,.T.',
        ),
      );
    }
    return out;
  }

  double _loopArea(BrepSolid s, BrepLoop loop, BrepPlane plane) {
    final pts = <Vector2>[];
    for (final coedge in loop.coedges) {
      final edge = s.edges[coedge.edge];
      final v = coedge.reversed ? edge.v1 : edge.v0;
      final uv = plane.parameters(s.vertices[v]);
      if (uv != null) pts.add(Vector2(uv.$1, uv.$2));
    }
    if (pts.length < 3) return 0;
    var sum = 0.0;
    for (var i = 0; i < pts.length; i++) {
      final a = pts[i];
      final b = pts[(i + 1) % pts.length];
      sum += a.x * b.y - b.x * a.y;
    }
    return sum.abs() / 2;
  }

  int _edgeLoop(
    BrepSolid s,
    BrepFace face,
    BrepLoop loop,
    Map<int, int> seams,
  ) {
    final oriented = <int>[];
    for (final coedge in loop.coedges) {
      if (!_edgeIds.containsKey(coedge.edge)) continue;
      final edge = s.edges[coedge.edge];
      final onSurface = _snapToSurface(face.surface, edge.curve);
      final id = onSurface == null
          ? _edgeIds[coedge.edge]!
          : seams.putIfAbsent(coedge.edge, () => _seamEdge(onSurface));
      oriented.add(
        w.add(
          'ORIENTED_EDGE',
          '${StepWriter.text('')},*,*,#$id,'
              '${coedge.reversed ? '.F.' : '.T.'}',
        ),
      );
    }
    if (oriented.isEmpty) {
      final c = s.edges[loop.coedges.first.edge].curve;
      return w.add('VERTEX_LOOP', "'',#${vertex(c.evaluate(c.t0))}");
    }
    return w.add(
      'EDGE_LOOP',
      '${StepWriter.text('')},${StepWriter.refs(oriented)}',
    );
  }

  int _seamEdge(BrepLine line) => w.add(
    'EDGE_CURVE',
    '${StepWriter.text('')},#${vertex(line.a)},#${vertex(line.b)},'
        '#${_curve(line)},.T.',
  );

  /// A seam line moved onto the surface it bounds, or null when the edge is
  /// already where it should be.
  ///
  /// A cylinder wall needs a seam, and this kernel stores that seam as the
  /// cylinder's **axis**. Inside the kernel that is a neat trick rather than a
  /// mistake: the axis has no angle, so its image in parameter space is
  /// exactly the u = 0 line, which is what the face's domain needs. In STEP it
  /// is a lie. A boundary edge is required to lie on its face's surface, and
  /// written literally the file describes a cylinder bounded by a line running
  /// up its middle — which this reader would take back in and get right again,
  /// while any other CAD system would reject the shell or quietly repair it.
  ///
  /// So the endpoints are projected out onto the surface. For a point on the
  /// axis that lands on the seam, because both this kernel's cylinder and its
  /// cone resolve the ambiguous direction the same way, towards the surface's
  /// own reference direction.
  static BrepLine? _snapToSurface(BrepSurface surface, BrepCurve curve) {
    if (curve is! BrepLine) return null;
    // Only where a seam is possible. On a plane an off-surface edge is a real
    // defect and moving it would hide that.
    if (surface is! BrepCylinder && surface is! BrepCone) return null;

    final a = surface.closestPoint(curve.a);
    final b = surface.closestPoint(curve.b);
    const tol = 1e-9;
    if (a.distanceToSquared(curve.a) < tol &&
        b.distanceToSquared(curve.b) < tol) {
      return null;
    }
    return a.distanceToSquared(b) < 1e-18 ? null : BrepLine(a, b);
  }

  /// A boundary for a face that has none, by cutting the surface open.
  ///
  /// A whole sphere or a whole torus needs no edges to be a closed solid, and
  /// this kernel stores them with none. STEP has no such thing: every face in
  /// a shell must carry at least one bound. The way round it, which is what
  /// every CAD system does, is a **seam**: an edge along the periodic
  /// parameter's zero line, traversed once in each direction, so it bounds the
  /// parametric rectangle without removing any surface.
  int _splineSeam(BrepSplineSurface s) {
    // Opposite coincident boundaries share an edge. Collapsed pole boundaries
    // contribute no finite edge; meridians carry the sphere's seam.
    final curves = [
      BrepSplineCurve(
        degree: s.degreeU,
        points: [for (final r in s.points) r.first],
        knots: s.knotsU,
        weights: [for (final r in s.weights) r.first],
      ),
      BrepSplineCurve(
        degree: s.degreeV,
        points: s.points.last,
        knots: s.knotsV,
        weights: s.weights.last,
      ),
      BrepSplineCurve(
        degree: s.degreeU,
        points: [for (final r in s.points) r.last],
        knots: s.knotsU,
        weights: [for (final r in s.weights) r.last],
      ),
      BrepSplineCurve(
        degree: s.degreeV,
        points: s.points.first,
        knots: s.knotsV,
        weights: s.weights.first,
      ),
    ];
    final stored = <(BrepSplineCurve, int)>[], uses = <(int, bool)>[];
    for (var i = 0; i < curves.length; i++) {
      final c = curves[i], a = c.evaluate(c.t0), b = c.evaluate(c.t1);
      if (c.points.every((p) => (p - a).length < 1e-9)) continue;
      int? id;
      for (final prior in stored) {
        if ([0.0, .25, .5, .75, 1.0].every(
          (t) =>
              (c.evaluate(c.t0 + (c.t1 - c.t0) * t) -
                      prior.$1.evaluate(
                        prior.$1.t0 + (prior.$1.t1 - prior.$1.t0) * t,
                      ))
                  .length <
              1e-8,
        )) {
          id = prior.$2;
          break;
        }
      }
      id ??= w.add(
        'EDGE_CURVE',
        "'',#${vertex(a)},#${vertex(b)},#${_curve(c)},.T.",
      );
      stored.add((c, id));
      uses.add((id, i < 2));
    }
    return _boundOf(uses);
  }

  int _seamLoop(BrepSurface surface) {
    switch (surface) {
      case BrepOffsetSurface():
        throw BrepUnsupported(
          'An offset surface must have explicit boundaries',
        );
      case BrepSplineSurface():
        return _splineSeam(surface);
      case BrepAffineSurface():
        throw BrepUnsupported(
          'STEP export of affine surfaces is not implemented',
        );
      case BrepRuledSurface():
        throw BrepUnsupported(
          'A ruled loft patch must have explicit boundaries',
        );
      case BrepSphere sphere:
        // Poles are degenerate in longitude, so the meridian at u = 0 is the
        // whole boundary: up one side of it, back down the other.
        final north = sphere.center + Vector3(0, 0, sphere.radius);
        final south = sphere.center - Vector3(0, 0, sphere.radius);
        final meridian = BrepCircle(
          center: sphere.center,
          normal: Vector3(0, 1, 0),
          xDir: Vector3(0, 0, 1),
          radius: sphere.radius,
          t0: 0,
          t1: math.pi,
        );
        final edge = w.add(
          'EDGE_CURVE',
          '${StepWriter.text('')},#${vertex(north)},#${vertex(south)},'
              '#${_curve(meridian)},.T.',
        );
        return _boundOf([(edge, true), (edge, false)]);

      case BrepTorus torus:
        // Periodic both ways, so it takes two seams. The loop runs around the
        // tube, around the ring, and back along each — the boundary of the
        // parameter square with both pairs of sides identified.
        final seamPoint =
            torus.center + torus.xDir * (torus.majorRadius + torus.minorRadius);
        final v = vertex(seamPoint);
        final tube = BrepCircle(
          center: torus.center + torus.xDir * torus.majorRadius,
          normal: -torus.yDir,
          xDir: torus.xDir,
          radius: torus.minorRadius,
        );
        final ring = BrepCircle(
          center: torus.center,
          normal: torus.axis,
          xDir: torus.xDir,
          radius: torus.majorRadius + torus.minorRadius,
        );
        final tubeEdge = w.add(
          'EDGE_CURVE',
          '${StepWriter.text('')},#$v,#$v,#${_curve(tube)},.T.',
        );
        final ringEdge = w.add(
          'EDGE_CURVE',
          '${StepWriter.text('')},#$v,#$v,#${_curve(ring)},.T.',
        );
        return _boundOf([
          (ringEdge, true),
          (tubeEdge, true),
          (ringEdge, false),
          (tubeEdge, false),
        ]);

      case BrepPlane _:
      case BrepCylinder _:
      case BrepCone _:
        // An unbounded plane is not a face, and a cylinder or cone with no
        // loops has no ends, so neither can appear in a closed shell. Reaching
        // here means the solid was already broken.
        throw ArgumentError(
          'A ${surface.runtimeType} face with no boundary cannot be part of a '
          'closed solid, so there is nothing to write',
        );
    }
  }

  int _boundOf(List<(int, bool)> edges) {
    final oriented = <int>[
      for (final (edge, forward) in edges)
        w.add(
          'ORIENTED_EDGE',
          '${StepWriter.text('')},*,*,#$edge,${forward ? '.T.' : '.F.'}',
        ),
    ];
    final loop = w.add(
      'EDGE_LOOP',
      '${StepWriter.text('')},${StepWriter.refs(oriented)}',
    );
    return w.add('FACE_OUTER_BOUND', '${StepWriter.text('')},#$loop,.T.');
  }
}
