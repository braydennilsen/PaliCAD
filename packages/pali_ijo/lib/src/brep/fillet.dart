import 'dart:math' as math;

import 'package:vector_math/vector_math_64.dart';

import '../ijo_log.dart';
import '../mesh/mesh_seal.dart';
import 'geom.dart';
import 'chamfer_options.dart';
import 'profile_blend.dart';
import 'classify.dart';
import 'convex_planes.dart';
import 'intersect.dart';
import 'slice.dart';
import 'heal.dart';
import 'primitives.dart';
import 'provenance.dart';
import 'solid.dart';
import 'tessellate.dart';
import 'xform.dart';

part 'blend_planar_regions.dart';

/// Analytic fillet / chamfer. Plane–plane → cylinder / plane;
/// plane–cylinder circular rim → torus / cone.
class BrepFillet {
  BrepFillet._();

  /// Rounds [edgeIds] to [radius]. With [inverse], cuts a circular groove
  /// centered on the selected edge instead of a tangent rounding.
  ///
  /// Provenance is carried across, which is what makes "select the fillet"
  /// answerable: the blend surfaces are tori and cylinders that existed in
  /// neither the input nor anywhere else, so they come back tagged as created
  /// by this fillet while the faces it merely re-trimmed keep their tags.
  static BrepSolid fillet(
    BrepSolid solid,
    List<int> edgeIds, {
    required double radius,
    required bool inverse,
  }) => BrepProvenance.inheritFrom(
    BrepHeal.mergeCoplanar(
      BrepHeal.sewClosedRims(
        _resolveBlendLineContacts(
          _fillet(solid, edgeIds, radius: radius, inverse: inverse),
        ),
      ),
    ),
    solid,
    op: 'fillet',
  );

  static BrepSolid _fillet(
    BrepSolid solid,
    List<int> edgeIds, {
    required double radius,
    required bool inverse,
  }) {
    if (!radius.isFinite || radius < 0) {
      throw BrepUnsupported('Fillet radius must be finite and nonnegative');
    }
    final transformed = BrepXform.unwrapAffine(solid);
    if (transformed != null) {
      final scale = _similarityScale(transformed.transform);
      if (scale != null) {
        return BrepXform.affine(
          _fillet(
            transformed.solid,
            edgeIds,
            radius: radius / scale,
            inverse: inverse,
          ),
          transformed.transform,
        );
      }
    }
    solid = BrepXform.materializeForEditing(solid);
    final r = radius;
    if (r < 1e-9) return solid;
    // A Boolean may partition a flat cap into a ring and a coplanar fill.
    // Blend the continuous face before its moved perimeter reaches that seam.
    solid = _blendPlanarRegions(BrepHeal.mergeFilledPlanarLoops(solid));
    edgeIds = _creaseEdges(solid, edgeIds);
    if (edgeIds.isEmpty) return solid;
    final unique = inverse
        ? edgeIds
              .where((id) => id >= 0 && id < solid.edges.length)
              .toSet()
              .toList()
        : completeCornerEdges(solid, edgeIds);
    IjoLog.line('brep fillet edges=${unique.length} r=$r inverse=$inverse');
    if (inverse) {
      final cap = _prismCapFace(solid, unique);
      if (cap != null) {
        return _blendPrismCap(solid, cap, r, round: true, inverse: true);
      }
    }
    if (!inverse) {
      final profile = BrepProfileBlend.tryRun(solid, edgeIds, r, round: true);
      if (profile != null) return profile;
      if (unique.length == 1) {
        final id = unique.first;
        if (id >= 0 && id < solid.edges.length) {
          return _filletOne(solid, id, r);
        }
      }
      final capUpright = _prismCapPlusUprights(solid, unique);
      if (capUpright != null) {
        return _blendPrismCap(
          solid,
          capUpright,
          r,
          round: true,
          filletUprights: true,
        );
      }
      final prismCap = _prismCapFace(solid, unique);
      if (prismCap != null) {
        return _blendPrismCap(solid, prismCap, r, round: true);
      }
      final aabbLoop = _aabbFaceLoop(solid, unique);
      if (aabbLoop != null) {
        return _blendAabbFaceLoop(solid, aabbLoop, r, round: true);
      }
    }
    return _walkBlend(solid, unique, r, inverse: inverse, round: true);
  }

  static BrepSolid chamfer(
    BrepSolid solid,
    List<int> edgeIds,
    double distance, {
    ChamferOptions options = const ChamferOptions(),
  }) => BrepProvenance.inheritFrom(
    _resolveBlendLineContacts(_chamfer(solid, edgeIds, distance, options)),
    solid,
    op: 'chamfer',
  );

  static BrepSolid _chamfer(
    BrepSolid solid,
    List<int> edgeIds,
    double distance,
    ChamferOptions options,
  ) {
    if (!distance.isFinite || distance < 0) {
      throw BrepUnsupported('Chamfer distance must be finite and nonnegative');
    }
    final transformed = BrepXform.unwrapAffine(solid);
    if (transformed != null) {
      final scale = _similarityScale(transformed.transform);
      if (scale != null) {
        return BrepXform.affine(
          _chamfer(
            transformed.solid,
            edgeIds,
            distance / scale,
            options.copyWith(distance2: options.distance2 / scale),
          ),
          transformed.transform,
        );
      }
    }
    solid = BrepXform.materializeForEditing(solid);
    final d = distance;
    if (d < 1e-9) return solid;
    solid = _blendPlanarRegions(BrepHeal.mergeFilledPlanarLoops(solid));
    edgeIds = _creaseEdges(solid, edgeIds);
    if (edgeIds.isEmpty) return solid;
    final profile = BrepProfileBlend.tryRun(
      solid,
      edgeIds,
      d,
      round: false,
      options: options,
    );
    if (profile != null) return profile;
    final unique = completeCornerEdges(solid, edgeIds);
    IjoLog.line('brep chamfer edges=${unique.length} d=$d');
    if (unique.length == 1) {
      final id = unique.first;
      if (id >= 0 && id < solid.edges.length) {
        return _chamferOne(solid, id, d, options: options);
      }
    }
    if (options.mode != ChamferMode.defaultMode) {
      return _walkBlend(
        solid,
        unique,
        d,
        inverse: false,
        round: false,
        options: options,
      );
    }
    final capUpright = _prismCapPlusUprights(solid, unique);
    if (capUpright != null) {
      return _blendPrismCap(
        solid,
        capUpright,
        d,
        round: false,
        filletUprights: true,
      );
    }
    final prismCap = _prismCapFace(solid, unique);
    if (prismCap != null) {
      return _blendPrismCap(solid, prismCap, d, round: false);
    }
    final aabbLoop = _aabbFaceLoop(solid, unique);
    if (aabbLoop != null) {
      return _blendAabbFaceLoop(solid, aabbLoop, d, round: false);
    }
    return _walkBlend(solid, unique, d, inverse: false, round: false);
  }

  /// One topology walk per connected chain. Plans every blend from the
  /// input solid, then rewrites; sequential [_filletOne] is only the
  /// fallback when a chain is a single isolated edge.
  static BrepSolid _walkBlend(
    BrepSolid solid,
    List<int> edgeIds,
    double r, {
    required bool inverse,
    required bool round,
    ChamferOptions options = const ChamferOptions(),
  }) {
    final valid = [
      for (final id in edgeIds)
        if (id >= 0 && id < solid.edges.length) id,
    ];
    if (valid.isEmpty) return solid;
    final extents = _aabbExtents(solid);
    if (round &&
        !inverse &&
        extents != null &&
        valid.toSet().length == solid.edges.length) {
      return _roundedAabb(extents.$1, extents.$2, r);
    }
    final corner = _aabbCornerTriple(solid, valid);
    if (corner != null && !inverse && round) {
      return _blendAabbCorner(solid, corner, r);
    }
    var result = solid;
    for (final originalChain in _edgeChains(solid, valid)) {
      // Each preceding chain rebuilds topology, so its output edge indices
      // cannot be addressed using indices from the original solid.
      final chain = [
        for (final id in originalChain)
          identical(result, solid)
              ? id
              : _matchCurve(result, solid.edges[id].curve, r),
      ];
      if (chain.any((id) => id < 0)) {
        throw BrepUnsupported(
          'Could not locate a remaining blend edge after rebuilding topology',
        );
      }
      if (chain.length == 1) {
        result = _blendOne(
          result,
          chain.first,
          r,
          inverse: inverse,
          round: round,
          options: options,
        );
        continue;
      }
      if (!inverse && options.mode == ChamferMode.defaultMode) {
        final capUpright = _prismCapPlusUprights(result, chain);
        if (capUpright != null) {
          result = _blendPrismCap(
            result,
            capUpright,
            r,
            round: round,
            filletUprights: true,
          );
          continue;
        }
        final prismCap = _prismCapFace(result, chain);
        if (prismCap != null) {
          result = _blendPrismCap(result, prismCap, r, round: round);
          continue;
        }
      }
      final seeds = [
        for (final id in chain)
          if (id >= 0 && id < result.edges.length) result.edges[id].curve,
      ];
      for (final curve in seeds) {
        final id = _matchCurve(result, curve, r);
        if (id < 0)
          throw BrepUnsupported(
            'A selected edge could not be located after the preceding blend',
          );
        result = _blendOne(
          result,
          id,
          r,
          inverse: inverse,
          round: round,
          options: options,
        );
      }
    }
    return result;
  }

  static BrepSolid _blendOne(
    BrepSolid s,
    int edgeId,
    double r, {
    required bool inverse,
    required bool round,
    ChamferOptions options = const ChamferOptions(),
  }) {
    if (!round) return _chamferOne(s, edgeId, r, options: options);
    return inverse
        ? _filletOne(s, edgeId, r, inverse: true)
        : _filletOne(s, edgeId, r);
  }

  static double? _similarityScale(BrepAffineTransform m) {
    final scale = m.x.length;
    if ((m.y.length - scale).abs() > 1e-10 * scale ||
        (m.z.length - scale).abs() > 1e-10 * scale ||
        m.x.dot(m.y).abs() > 1e-10 * scale * scale ||
        m.x.dot(m.z).abs() > 1e-10 * scale * scale ||
        m.y.dot(m.z).abs() > 1e-10 * scale * scale)
      return null;
    return scale;
  }

  static List<List<int>> _edgeChains(BrepSolid s, List<int> edgeIds) {
    final wanted = edgeIds.toSet();
    final atVertex = <int, List<int>>{};
    for (final id in wanted) {
      final e = s.edges[id];
      (atVertex[e.v0] ??= []).add(id);
      if (e.v1 != e.v0) (atVertex[e.v1] ??= []).add(id);
    }
    final adj = <int, List<int>>{
      for (final id in wanted)
        id: {
          ...atVertex[s.edges[id].v0]!,
          ...atVertex[s.edges[id].v1]!,
        }.where((other) => other != id).toList(),
    };
    final seen = <int>{};
    final chains = <List<int>>[];
    for (final start in wanted) {
      if (!seen.add(start)) continue;
      final chain = <int>[start];
      final stack = [start];
      while (stack.isNotEmpty) {
        final cur = stack.removeLast();
        for (final n in adj[cur] ?? const <int>[]) {
          if (seen.add(n)) {
            chain.add(n);
            stack.add(n);
          }
        }
      }
      chains.add(chain);
    }
    return chains;
  }

  /// Construct all twelve box rounds together so adjacent blends share their
  /// spherical corners instead of repeatedly trimming an already changed box.
  static BrepSolid _roundedAabb(Vector3 low, Vector3 high, double radius) {
    if (radius * 2 >=
        math.min(high.x - low.x, math.min(high.y - low.y, high.z - low.z))) {
      throw BrepUnsupported(
        'Fillet radius must be less than half the box width',
      );
    }
    final lo = low + Vector3.all(radius), hi = high - Vector3.all(radius);
    final axes = [Vector3(1, 0, 0), Vector3(0, 1, 0), Vector3(0, 0, 1)];
    final b = BrepBuilder();
    int vertex(Vector3 p) {
      for (var i = 0; i < b.vertices.length; i++) {
        if ((b.vertices[i] - p).length < 1e-8) return i;
      }
      return b.vertex(p);
    }

    BrepCoedge use(BrepCurve curve) {
      final p = curve.evaluate(curve.t0), q = curve.evaluate(curve.t1);
      final i = vertex(p),
          j = vertex(q),
          mid = curve.evaluate((curve.t0 + curve.t1) / 2);
      for (var k = 0; k < b.edges.length; k++) {
        final e = b.edges[k];
        if ((e.v0 == i && e.v1 == j || e.v0 == j && e.v1 == i) &&
            (e.curve.evaluate((e.curve.t0 + e.curve.t1) / 2) - mid).length <
                1e-8) {
          return BrepCoedge(k, reversed: e.v0 != i);
        }
      }
      return BrepCoedge(b.edge(i, j, curve));
    }

    void face(BrepSurface surface, List<BrepCurve> curves) {
      var loop = curves.map(use).toList();
      final pts = [for (final c in curves) c.evaluate(c.t0)];
      var area = Vector3.zero(), center = Vector3.zero();
      for (var i = 0; i < pts.length; i++) {
        area += pts[i].cross(pts[(i + 1) % pts.length]);
        center += pts[i];
      }
      final uv = surface.parameters(center / pts.length.toDouble())!;
      if (area.dot(surface.normalAt(uv.$1, uv.$2)) < 0) {
        loop = [
          for (final c in loop.reversed)
            BrepCoedge(c.edge, reversed: !c.reversed),
        ];
      }
      b.face(surface, [BrepLoop(loop)]);
    }

    BrepCurve arc(Vector3 center, Vector3 from, Vector3 to) => BrepCircle(
      center: center,
      normal: from.cross(to).normalized(),
      radius: radius,
      xDir: from,
      t0: 0,
      t1: math.pi / 2,
    );
    for (var axis = 0; axis < 3; axis++) {
      final j = (axis + 1) % 3, k = (axis + 2) % 3;
      for (final sign in [-1, 1]) {
        final normal = axes[axis] * sign.toDouble();
        final pts = <Vector3>[];
        for (final pair in [(0, 0), (1, 0), (1, 1), (0, 1)]) {
          final p = Vector3.zero();
          p[axis] = sign < 0 ? low[axis] : high[axis];
          p[j] = pair.$1 == 0 ? lo[j] : hi[j];
          p[k] = pair.$2 == 0 ? lo[k] : hi[k];
          pts.add(p);
        }
        face(BrepPlane(origin: pts.first, normal: normal), [
          for (var i = 0; i < 4; i++) BrepLine(pts[i], pts[(i + 1) % 4]),
        ]);
      }
      for (final sj in [-1, 1])
        for (final sk in [-1, 1]) {
          final n1 = axes[j] * sj.toDouble(), n2 = axes[k] * sk.toDouble();
          final c0 = lo.clone();
          c0[j] = sj < 0 ? lo[j] : hi[j];
          c0[k] = sk < 0 ? lo[k] : hi[k];
          final c1 = c0.clone()..[axis] = hi[axis];
          face(
            BrepCylinder(
              origin: c0,
              axis: axes[axis],
              radius: radius,
              xDir: n1,
            ),
            [
              BrepLine(c0 + n1 * radius, c1 + n1 * radius),
              arc(c1, n1, n2),
              BrepLine(c1 + n2 * radius, c0 + n2 * radius),
              arc(c0, n2, n1),
            ],
          );
        }
    }
    for (final sx in [-1, 1])
      for (final sy in [-1, 1])
        for (final sz in [-1, 1]) {
          final center = Vector3(
            sx < 0 ? lo.x : hi.x,
            sy < 0 ? lo.y : hi.y,
            sz < 0 ? lo.z : hi.z,
          );
          final ns = [
            axes[0] * sx.toDouble(),
            axes[1] * sy.toDouble(),
            axes[2] * sz.toDouble(),
          ];
          face(BrepSphere(center: center, radius: radius), [
            for (var i = 0; i < 3; i++) arc(center, ns[i], ns[(i + 1) % 3]),
          ]);
        }
    return b.build();
  }

  static List<int>? _aabbCornerTriple(BrepSolid s, List<int> edges) {
    if (edges.length != 3 || _aabbExtents(s) == null) return null;
    final counts = <int, int>{};
    for (final id in edges) {
      final e = s.edges[id];
      counts[e.v0] = (counts[e.v0] ?? 0) + 1;
      counts[e.v1] = (counts[e.v1] ?? 0) + 1;
    }
    final corner = counts.entries.where((e) => e.value == 3);
    if (corner.length != 1) return null;
    return edges;
  }

  /// Three concurrent AABB edges → two-edge spherical corner, then the
  /// remaining upright when that rewrite stays closed.
  static BrepSolid _blendAabbCorner(BrepSolid s, List<int> edges, double r) {
    BrepSolid? rounded;
    for (var i = 0; i < edges.length; i++) {
      final pair = [edges[i], edges[(i + 1) % edges.length]];
      final loop = _aabbFaceLoop(s, pair);
      if (loop == null) continue;
      final nSel =
          (loop.f ? 1 : 0) +
          (loop.r ? 1 : 0) +
          (loop.b ? 1 : 0) +
          (loop.l ? 1 : 0);
      if (nSel != 2) continue;
      final base = _blendAabbFaceLoop(s, loop, r, round: true);
      rounded ??= base;
      final leftover = edges[(i + 2) % edges.length];
      final seed = s.edges[leftover].curve;
      final id = _matchCurve(base, seed, r);
      if (id < 0 || !_canBlend(base, id)) continue;
      try {
        final next = _filletOne(base, id, r);
        if (MeshSeal.isClosedSolid(BrepTessellate.run(next))) {
          return next;
        }
      } on BrepUnsupported {
        continue;
      }
    }
    if (rounded != null) return rounded;
    var result = s;
    for (final id in edges) {
      if (id < 0 || id >= result.edges.length) continue;
      if (!_canBlend(result, id)) continue;
      result = _filletOne(result, id, r);
    }
    return result;
  }

  static (Vector3 min, Vector3 max)? _aabbExtents(BrepSolid s) {
    if (s.vertices.length != 8 || s.faces.length != 6) return null;
    for (final f in s.faces) {
      if (f.surface is! BrepPlane) return null;
    }
    var min = s.vertices.first.clone();
    var max = s.vertices.first.clone();
    for (final v in s.vertices) {
      min = Vector3(
        math.min(min.x, v.x),
        math.min(min.y, v.y),
        math.min(min.z, v.z),
      );
      max = Vector3(
        math.max(max.x, v.x),
        math.max(max.y, v.y),
        math.max(max.z, v.z),
      );
    }
    for (final v in s.vertices) {
      final onX = (v.x - min.x).abs() < 1e-7 || (v.x - max.x).abs() < 1e-7;
      final onY = (v.y - min.y).abs() < 1e-7 || (v.y - max.y).abs() < 1e-7;
      final onZ = (v.z - min.z).abs() < 1e-7 || (v.z - max.z).abs() < 1e-7;
      if (!(onX && onY && onZ)) return null;
    }
    return (min, max);
  }

  /// Edges of one AABB face, or null. [sides] is front/right/back/left of
  /// that face after mapping it to +Z.
  static ({int axis, bool high, bool f, bool r, bool b, bool l})? _aabbFaceLoop(
    BrepSolid s,
    List<int> edges,
  ) {
    if (edges.length < 2 || edges.length > 4) return null;
    final aabb = _aabbExtents(s);
    if (aabb == null) return null;
    var bestFace = -1;
    for (var fi = 0; fi < s.faces.length; fi++) {
      final ids = <int>{};
      for (final loop in s.faces[fi].loops) {
        for (final c in loop.coedges) {
          ids.add(c.edge);
        }
      }
      if (edges.every(ids.contains)) {
        bestFace = fi;
        break;
      }
    }
    if (bestFace < 0) return null;
    final pl = s.faces[bestFace].surface;
    if (pl is! BrepPlane) return null;
    final n = pl.normal;
    late final int axis;
    late final bool high;
    if (n.z.abs() > 0.9) {
      axis = 2;
      high = n.z > 0;
    } else if (n.y.abs() > 0.9) {
      axis = 1;
      high = n.y > 0;
    } else if (n.x.abs() > 0.9) {
      axis = 0;
      high = n.x > 0;
    } else {
      return null;
    }
    final min = aabb.$1;
    final max = aabb.$2;
    final sx = switch (axis) {
      0 => max.y - min.y,
      1 => max.x - min.x,
      _ => max.x - min.x,
    };
    final sy = switch (axis) {
      0 => max.z - min.z,
      1 => max.z - min.z,
      _ => max.y - min.y,
    };
    final sz = switch (axis) {
      0 => max.x - min.x,
      1 => max.y - min.y,
      _ => max.z - min.z,
    };
    Vector3 toLocal(Vector3 p) {
      switch (axis) {
        case 0:
          return Vector3(
            p.y - min.y,
            p.z - min.z,
            high ? p.x - min.x : max.x - p.x,
          );
        case 1:
          return Vector3(
            p.x - min.x,
            p.z - min.z,
            high ? p.y - min.y : max.y - p.y,
          );
        default:
          return Vector3(
            p.x - min.x,
            p.y - min.y,
            high ? p.z - min.z : max.z - p.z,
          );
      }
    }

    var f = false, rgt = false, b = false, l = false;
    for (final id in edges) {
      if (id < 0 || id >= s.edges.length) continue;
      final e = s.edges[id];
      final mid = (s.vertices[e.v0] + s.vertices[e.v1]) * 0.5;
      final q = toLocal(mid);
      if ((q.z - sz).abs() > 1e-4) continue;
      if (q.y.abs() < 1e-4) {
        f = true;
      } else if ((q.y - sy).abs() < 1e-4) {
        b = true;
      } else if (q.x.abs() < 1e-4) {
        l = true;
      } else if ((q.x - sx).abs() < 1e-4) {
        rgt = true;
      }
    }
    if (!(f || rgt || b || l)) return null;
    return (axis: axis, high: high, f: f, r: rgt, b: b, l: l);
  }

  static BrepSolid _blendAabbFaceLoop(
    BrepSolid s,
    ({int axis, bool high, bool f, bool r, bool b, bool l}) loop,
    double r, {
    required bool round,
  }) {
    final aabb = _aabbExtents(s)!;
    final min = aabb.$1;
    final max = aabb.$2;
    final dx = max.x - min.x;
    final dy = max.y - min.y;
    final dz = max.z - min.z;
    late final double sx;
    late final double sy;
    late final double sz;
    switch (loop.axis) {
      case 0:
        sx = dy;
        sy = dz;
        sz = dx;
      case 1:
        sx = dx;
        sy = dz;
        sz = dy;
      default:
        sx = dx;
        sy = dy;
        sz = dz;
    }
    final nSel =
        (loop.f ? 1 : 0) +
        (loop.r ? 1 : 0) +
        (loop.b ? 1 : 0) +
        (loop.l ? 1 : 0);
    final opposite =
        (loop.f && loop.b && !loop.l && !loop.r) ||
        (loop.l && loop.r && !loop.f && !loop.b);
    final cap = opposite
        ? math.min((loop.f || loop.b) ? sy : sx, sz * 2) / 2 - 1e-6
        : nSel >= 4
        ? math.min(sx, math.min(sy, sz * 2)) / 2 - 1e-6
        : math.min(sx, math.min(sy, sz)) - 1e-6;
    if (r >= cap || r >= sz - 1e-6) {
      throw BrepUnsupported(
        'Blend radius is too large for this face '
        '(max ${cap.toStringAsFixed(2)} mm)',
      );
    }
    BrepSolid local;
    if (nSel >= 4) {
      local = _canonicalFaceBlend(sx, sy, sz, r, round: round);
    } else if (opposite) {
      if (loop.f || loop.b) {
        local = _canonicalOppositeBlend(sx, sy, sz, r, round: round);
      } else {
        local = BrepXform.map(
          _canonicalOppositeBlend(sy, sx, sz, r, round: round),
          (p) => Vector3(p.y, sx - p.x, p.z),
          reverseWinding: true,
        );
      }
    } else if (nSel == 2) {
      local = _canonicalAdjacentFrontLeft(sx, sy, sz, r, round: round);
      if (loop.f && loop.l) {
        // identity
      } else if (loop.f && loop.r) {
        local = BrepXform.map(
          local,
          (p) => Vector3(sx - p.x, p.y, p.z),
          reverseWinding: true,
        );
      } else if (loop.b && loop.l) {
        local = BrepXform.map(
          local,
          (p) => Vector3(p.x, sy - p.y, p.z),
          reverseWinding: true,
        );
      } else if (loop.b && loop.r) {
        local = BrepXform.map(local, (p) => Vector3(sx - p.x, sy - p.y, p.z));
      } else {
        local = _canonicalFaceBlend(sx, sy, sz, r, round: round);
      }
    } else {
      // Three sides of the face: reuse the full loop (extra edge is
      // cheaper than a glitched sequential vertex).
      local = _canonicalFaceBlend(sx, sy, sz, r, round: round);
    }
    Vector3 toWorld(Vector3 p) {
      final x = p.x, y = p.y, z = loop.high ? p.z : sz - p.z;
      switch (loop.axis) {
        case 0:
          return Vector3(
            loop.high ? min.x + z : max.x - z,
            min.y + x,
            min.z + y,
          );
        case 1:
          return Vector3(
            min.x + x,
            loop.high ? min.y + z : max.y - z,
            min.z + y,
          );
        default:
          return Vector3(
            min.x + x,
            min.y + y,
            loop.high ? min.z + z : max.z - z,
          );
      }
    }

    return BrepXform.map(local, toWorld, reverseWinding: !loop.high);
  }

  /// Full line-loop of a planar prism cap (n ≥ 3), or null.
  static int? _prismCapFace(BrepSolid s, List<int> edges) {
    if (edges.length < 3) return null;
    var best = -1;
    for (var fi = 0; fi < s.faces.length; fi++) {
      final face = s.faces[fi];
      if (face.surface is! BrepPlane || face.loops.isEmpty) continue;
      final ids = <int>{for (final c in face.loops.first.coedges) c.edge};
      if (ids.length == edges.toSet().length && edges.every(ids.contains)) {
        best = fi;
        break;
      }
    }
    if (best < 0) return null;
    if (_prismLegs(s, best) == null) return null;
    return best;
  }

  /// Cap loop plus the wall leg at every cap vertex — the face-fillet set.
  static int? _prismCapPlusUprights(BrepSolid s, List<int> edges) {
    final selected = edges.toSet();
    if (selected.length < 6) return null;
    for (var fi = 0; fi < s.faces.length; fi++) {
      final face = s.faces[fi];
      if (face.surface is! BrepPlane || face.loops.isEmpty) continue;
      final capIds = <int>{for (final c in face.loops.first.coedges) c.edge};
      if (capIds.length < 3 || !capIds.every(selected.contains)) continue;
      final legs = _prismLegs(s, fi);
      if (legs == null) continue;
      final extra = selected.difference(capIds);
      if (extra.length == legs.length && extra.every(legs.contains)) {
        return fi;
      }
    }
    return null;
  }

  /// One non-cap edge per cap vertex, each nearly parallel to the face normal.
  static Set<int>? _prismLegs(BrepSolid s, int faceIndex) {
    final face = s.faces[faceIndex];
    if (face.surface is! BrepPlane || face.loops.isEmpty) return null;
    for (final c in face.loops.first.coedges) {
      if (s.edges[c.edge].curve is! BrepLine) return null;
    }
    final capIds = _loopVertexIds(s, face.loops.first);
    if (capIds.length < 3) return null;
    final nHat = face.outwardNormal(0, 0);
    final capSet = {for (final c in face.loops.first.coedges) c.edge};
    final legs = <int>{};
    for (final vi in capIds) {
      var leg = -1;
      for (var e = 0; e < s.edges.length; e++) {
        if (capSet.contains(e)) continue;
        final ed = s.edges[e];
        if (ed.v0 != vi && ed.v1 != vi) continue;
        if (leg >= 0) return null;
        leg = e;
      }
      if (leg < 0) return null;
      final ed = s.edges[leg];
      final other = ed.v0 == vi ? ed.v1 : ed.v0;
      final dir = (s.vertices[vi] - s.vertices[other]).normalized();
      if (dir.dot(nHat) < 0.85) return null;
      legs.add(leg);
    }
    return legs;
  }

  static List<int> _loopVertexIds(BrepSolid s, BrepLoop loop) {
    return [
      for (final c in loop.coedges)
        c.reversed ? s.edges[c.edge].v1 : s.edges[c.edge].v0,
    ];
  }

  static Vector2? _intersect2(Vector2 a, Vector2 b, Vector2 c, Vector2 d) {
    final r = b - a;
    final q = d - c;
    final den = r.x * q.y - r.y * q.x;
    if (den.abs() < 1e-14) return null;
    final t = ((c.x - a.x) * q.y - (c.y - a.y) * q.x) / den;
    return Vector2(a.x + r.x * t, a.y + r.y * t);
  }

  /// Rebuild a prism with one cap fully blended (any n ≥ 3).
  ///
  /// When [filletUprights] is set, each wall leg becomes a cylinder (or
  /// chamfer plane) too: n cap blends + n uprights + n spheres at the cap
  /// vertices. The opposite face stays planar; its edges stay lines, with
  /// quarter-circles only where the uprights meet it.
  static BrepSolid _blendPrismCap(
    BrepSolid s,
    int faceIndex,
    double r, {
    required bool round,
    bool filletUprights = false,
    bool inverse = false,
  }) {
    final face = s.faces[faceIndex];
    final plane = face.surface as BrepPlane;
    var nHat = face.outwardNormal(0, 0);
    final capIds = _loopVertexIds(s, face.loops.first);
    final n = capIds.length;
    final C = [for (final i in capIds) s.vertices[i].clone()];
    final capSet = {for (final c in face.loops.first.coedges) c.edge};
    final B = <Vector3>[];
    for (final vi in capIds) {
      for (var e = 0; e < s.edges.length; e++) {
        if (capSet.contains(e)) continue;
        final ed = s.edges[e];
        if (ed.v0 != vi && ed.v1 != vi) continue;
        B.add(s.vertices[ed.v0 == vi ? ed.v1 : ed.v0].clone());
        break;
      }
    }
    if (B.length != n) {
      throw BrepUnsupported('Cap fillet needs a prism wall at each vertex');
    }
    // This rebuild covers the local prism, not necessarily the whole body.
    // A cube may be attached through a hole in its opposite cap to a sphere
    // or another feature. Keep that cap's inner loops and the attached faces.
    final rebuiltFaces = <int>{faceIndex};
    for (final co in face.loops.first.coedges) {
      rebuiltFaces.addAll(facesForEdge(s, co.edge));
    }
    int? opposite;
    for (var fi = 0; fi < s.faces.length; fi++) {
      if (rebuiltFaces.contains(fi)) continue;
      final f = s.faces[fi];
      if (f.surface is! BrepPlane || f.loops.isEmpty) continue;
      final ids = _loopVertexIds(s, f.loops.first);
      if (ids.length == n &&
          B.every((p) => ids.any((id) => (s.vertices[id] - p).length < 1e-7))) {
        opposite = fi;
        break;
      }
    }
    if (opposite == null ||
        rebuiltFaces.length != n + 1 ||
        rebuiltFaces.any(
          (fi) =>
              s.faces[fi].loops.length != 1 ||
              (fi != faceIndex && s.faces[fi].loops.first.coedges.length != 4),
        )) {
      if (!round && BrepConvexPlanes.supports(s)) {
        // Existing bevels change wall topology. Chamfer the requested cap
        // boundary without re-chamfering completed upright corner edges.
        return _chamferConvexEdges(
          s,
          capSet.toList(),
          r,
          const ChamferOptions(),
        );
      }
      throw BrepUnsupported(
        'Face-loop blend needs intact prism walls and an opposite cap',
      );
    }
    rebuiltFaces.add(opposite);
    var area = Vector3.zero();
    for (var i = 0; i < n; i++) {
      area += C[i].cross(C[(i + 1) % n]);
    }
    if (area.dot(nHat) < 0) {
      final cr = C.reversed.toList();
      final br = B.reversed.toList();
      C
        ..clear()
        ..addAll(cr);
      B
        ..clear()
        ..addAll(br);
    }
    final origin = C.first;
    var u = plane.xDir.clone();
    var v = nHat.cross(u);
    if (v.length2 < 1e-16) {
      u = nHat.cross(Vector3(1, 0, 0));
      if (u.length2 < 1e-12) u = nHat.cross(Vector3(0, 1, 0));
      u.normalize();
      v = nHat.cross(u).normalized();
    } else {
      v.normalize();
    }
    final uv = [
      for (final p in C) Vector2((p - origin).dot(u), (p - origin).dot(v)),
    ];
    var uvArea = 0.0;
    for (var i = 0; i < n; i++) {
      final a = uv[i];
      final b = uv[(i + 1) % n];
      uvArea += a.x * b.y - b.x * a.y;
    }
    if (uvArea < 0) {
      u = -u;
      v = -v;
      for (var i = 0; i < n; i++) {
        uv[i] = Vector2((C[i] - origin).dot(u), (C[i] - origin).dot(v));
      }
    }
    final insets = <double>[];
    for (var i = 0; i < n; i++) {
      final e = C[(i + 1) % n] - C[i];
      final M = e.cross(nHat);
      if (M.length < 1e-12) {
        throw BrepUnsupported('Cap edge is degenerate');
      }
      M.normalize();
      final wallN = M;
      final beta = math.acos(nHat.dot(wallN).clamp(-1.0, 1.0));
      final alpha = math.pi - beta;
      final half = alpha * 0.5;
      if (half.abs() < 1e-6 || math.tan(half).abs() < 1e-9) {
        throw BrepUnsupported('Cap fillet angle is degenerate');
      }
      insets.add(r / math.tan(half));
    }
    final off = <(Vector2, Vector2)>[];
    for (var i = 0; i < n; i++) {
      final a = uv[i];
      final b = uv[(i + 1) % n];
      final e = b - a;
      final len = e.length;
      if (len < 1e-12) {
        throw BrepUnsupported('Cap fillet profile is degenerate');
      }
      final inward = Vector2(-e.y, e.x) / len * insets[i];
      off.add((a + inward, b + inward));
    }
    final insetUv = <Vector2>[];
    for (var i = 0; i < n; i++) {
      final prev = off[(i - 1 + n) % n];
      final cur = off[i];
      final hit = _intersect2(prev.$1, prev.$2, cur.$1, cur.$2);
      if (hit == null) {
        throw BrepUnsupported('Cap fillet inset failed');
      }
      insetUv.add(hit);
    }
    var insetArea = 0.0;
    for (var i = 0; i < n; i++) {
      final a = insetUv[i];
      final b = insetUv[(i + 1) % n];
      insetArea += a.x * b.y - b.x * a.y;
    }
    if (insetArea <= 1e-8) {
      throw BrepUnsupported('Blend radius is too large for this face');
    }
    var minH = double.infinity;
    for (var i = 0; i < n; i++) {
      final h = (C[i] - B[i]).dot(nHat);
      if (h < minH) minH = h;
    }
    if (r >= minH - 1e-6) {
      throw BrepUnsupported(
        'Blend radius is too large for this face '
        '(max ${(minH - 1e-6).toStringAsFixed(2)} mm)',
      );
    }
    final T = [for (final q in insetUv) origin + u * q.x + v * q.y];
    final Q = [for (final p in T) p - nHat * r];
    final S = [for (final p in C) p - nHat * r];
    if (inverse) {
      // The miter seams below require a convex, straight prism. Check every
      // inset against every wall; signed area alone misses collapsed offsets.
      for (var i = 0; i < n; i++) {
        final edge = C[(i + 1) % n] - C[i];
        final inward = nHat.cross(edge).normalized();
        final leg = C[i] - B[i];
        if (leg.cross(nHat).length > 1e-7 ||
            (T[(i + 1) % n] - T[i]).dot(edge) <= 1e-8 ||
            T.any((p) => (p - C[i]).dot(inward) < r - 1e-7)) {
          throw BrepUnsupported(
            'Inverse face fillet needs a convex straight prism and a radius that fits its face',
          );
        }
      }
      return _inversePrismCap(
        s,
        faceIndex,
        opposite,
        rebuiltFaces,
        C,
        B,
        T,
        S,
        nHat,
        r,
      );
    }
    final M = <Vector3>[];
    final w0 = <Vector3>[];
    final w1 = <Vector3>[];
    for (var i = 0; i < n; i++) {
      final e = C[(i + 1) % n] - C[i];
      final m = e.cross(nHat).normalized();
      M.add(m);
      w0.add(Q[i] + m * r);
      w1.add(Q[(i + 1) % n] + m * r);
    }
    final heights = [for (var i = 0; i < n; i++) (C[i] - B[i]).dot(nHat)];
    final U = [for (var i = 0; i < n; i++) T[i] - nHat * heights[i]];
    final bt0 = [for (var i = 0; i < n; i++) U[i] + M[i] * r];
    final bt1 = [for (var i = 0; i < n; i++) U[i] + M[(i - 1 + n) % n] * r];
    if (filletUprights) {
      for (var i = 0; i < n; i++) {
        final j = (i + 1) % n;
        final along = C[j] - C[i];
        if (along.length < 1e-12) {
          throw BrepUnsupported('Cap fillet axis is degenerate');
        }
        final remain = (bt1[j] - bt0[i]).dot(along.normalized());
        if (remain < 1e-6) {
          throw BrepUnsupported('Blend radius is too large for this face');
        }
      }
    }

    final bld = BrepBuilder();
    int vtx(Vector3 p) => bld.vertex(p);
    final bI = filletUprights ? const <int>[] : [for (final p in B) vtx(p)];
    final sI = filletUprights ? const <int>[] : [for (final p in S) vtx(p)];
    final tI = [for (final p in T) vtx(p)];
    final w0I = [for (final p in w0) vtx(p)];
    final w1I = [for (final p in w1) vtx(p)];
    final bt0I = filletUprights ? [for (final p in bt0) vtx(p)] : const <int>[];
    final bt1I = filletUprights ? [for (final p in bt1) vtx(p)] : const <int>[];

    final lineOf = <String, int>{};
    int line(int i, int j) {
      final lo = i < j ? i : j;
      final hi = i < j ? j : i;
      final k = '$lo-$hi';
      final existing = lineOf[k];
      if (existing != null) return existing;
      final id = bld.edge(i, j, BrepLine(bld.vertices[i], bld.vertices[j]));
      lineOf[k] = id;
      return id;
    }

    int arc(
      int i,
      int j,
      Vector3 center,
      Vector3 normal,
      Vector3 xDir,
      double t0,
      double t1,
    ) {
      if (!round) return line(i, j);
      return bld.edge(
        i,
        j,
        BrepCircle(
          center: center,
          normal: normal,
          radius: r,
          xDir: xDir,
          t0: t0,
          t1: t1,
        ),
      );
    }

    BrepLoop loop(List<(int, bool)> es) =>
        BrepLoop([for (final e in es) BrepCoedge(e.$1, reversed: e.$2)]);

    final eBot = <int>[];
    final eVert = <int>[];
    final eTop = <int>[];
    final eWall = <int>[];
    final eS0 = <int>[];
    final eS1 = <int>[];
    final eGen0 = <int>[];
    final eGen1 = <int>[];
    final aCap = <int>[];
    final aWall = <int>[];
    final aBot = <int>[];
    for (var i = 0; i < n; i++) {
      final j = (i + 1) % n;
      final iPrev = (i - 1 + n) % n;
      eTop.add(line(tI[i], tI[j]));
      eWall.add(line(w0I[i], w1I[i]));
      final ax = (Q[j] - Q[i]);
      if (ax.length < 1e-12) {
        throw BrepUnsupported('Cap fillet axis is degenerate');
      }
      final axis = ax.normalized();
      aCap.add(arc(tI[i], w0I[i], Q[i], axis, nHat, 0, math.pi / 2));
      aWall.add(arc(tI[j], w1I[i], Q[j], axis, nHat, 0, math.pi / 2));
      if (filletUprights) {
        eBot.add(line(bt0I[i], bt1I[j]));
        eGen0.add(line(bt0I[i], w0I[i]));
        eGen1.add(line(bt1I[i], w1I[iPrev]));
        aBot.add(
          arc(
            bt0I[i],
            bt1I[i],
            U[i],
            -nHat,
            M[i],
            0,
            _arcToward(M[i], M[iPrev], -nHat, B[i] - U[i]),
          ),
        );
      } else {
        eBot.add(line(bI[i], bI[j]));
        eVert.add(line(bI[i], sI[i]));
        eS0.add(line(sI[i], w0I[i]));
        eS1.add(line(w1I[i], sI[j]));
      }
    }

    // Bottom: outward -nHat, opposite of cap CCW.
    if (filletUprights) {
      bld.face(BrepPlane(origin: U.first, normal: -nHat), [
        loop([
          for (var i = n - 1; i >= 0; i--) ...[
            (eBot[i], true),
            (aBot[i], false),
          ],
        ]),
      ]);
    } else {
      bld.face(BrepPlane(origin: B.first, normal: -nHat), [
        loop([for (var i = n - 1; i >= 0; i--) (eBot[i], true)]),
      ]);
    }
    bld.face(BrepPlane(origin: T.first, normal: nHat), [
      loop([for (final e in eTop) (e, false)]),
    ]);

    for (var i = 0; i < n; i++) {
      final j = (i + 1) % n;
      if (filletUprights) {
        bld.face(BrepPlane(origin: C[i], normal: M[i]), [
          loop([
            (eBot[i], false),
            (eGen1[j], false),
            (eWall[i], true),
            (eGen0[i], true),
          ]),
        ]);
      } else {
        bld.face(BrepPlane(origin: C[i], normal: M[i]), [
          loop([
            (eBot[i], false),
            (eVert[j], false),
            (eS1[i], true),
            (eWall[i], true),
            (eS0[i], true),
            (eVert[i], true),
          ]),
        ]);
      }
      BrepSurface blendSurf;
      if (round) {
        blendSurf = BrepCylinder(
          origin: Q[i],
          axis: (Q[j] - Q[i]).normalized(),
          radius: r,
          xDir: nHat,
        );
      } else {
        final mid = (T[i] + w0[i]) * 0.5;
        var pn = (T[j] - T[i]).cross(w0[i] - T[i]);
        if (pn.length < 1e-12) pn = M[i] + nHat;
        pn.normalize();
        if (pn.dot(M[i] + nHat) < 0) pn = -pn;
        blendSurf = BrepPlane(origin: mid, normal: pn);
      }
      bld.face(blendSurf, [
        loop([
          (eTop[i], false),
          (aWall[i], false),
          (eWall[i], true),
          (aCap[i], true),
        ]),
      ]);
    }

    if (round) {
      for (var i = 0; i < n; i++) {
        final iPrev = (i - 1 + n) % n;
        final hArc = arc(
          w0I[i],
          w1I[iPrev],
          Q[i],
          nHat,
          M[i],
          0,
          _arcToward(M[i], M[iPrev], nHat, S[i] - Q[i]),
        );
        if (filletUprights) {
          bld.face(BrepSphere(center: Q[i], radius: r), [
            loop([(aWall[iPrev], false), (hArc, true), (aCap[i], true)]),
          ]);
          BrepSurface uprightSurf = BrepCylinder(
            origin: U[i],
            axis: nHat,
            radius: r,
            xDir: M[i],
          );
          bld.face(uprightSurf, [
            loop([
              (eGen0[i], false),
              (hArc, false),
              (eGen1[i], true),
              (aBot[i], true),
            ]),
          ]);
        } else {
          bld.face(BrepSphere(center: Q[i], radius: r), [
            loop([(aCap[i], false), (hArc, false), (aWall[iPrev], true)]),
          ]);
          bld.face(BrepPlane(origin: S[i], normal: nHat), [
            loop([(eS0[i], false), (hArc, false), (eS1[iPrev], false)]),
          ]);
        }
      }
    } else {
      for (var i = 0; i < n; i++) {
        final iPrev = (i - 1 + n) % n;
        final e0 = line(tI[i], w0I[i]);
        final e1 = line(w0I[i], w1I[iPrev]);
        // This edge already exists as aWall[iPrev], directed from the cap.
        final e2 = aWall[iPrev];
        var pn = (w0[i] - T[i]).cross(w1[iPrev] - T[i]);
        if (pn.length < 1e-12) continue;
        pn.normalize();
        if (pn.dot(nHat) < 0) pn = -pn;
        if (filletUprights) {
          bld.face(BrepPlane(origin: T[i], normal: pn), [
            loop([(e2, false), (e1, true), (e0, true)]),
          ]);
          var un = M[i] + M[iPrev];
          if (un.length < 1e-12) un = B[i] - U[i];
          un.normalize();
          bld.face(BrepPlane(origin: bt0[i], normal: un), [
            loop([
              (eGen0[i], false),
              (e1, false),
              (eGen1[i], true),
              (aBot[i], true),
            ]),
          ]);
        } else {
          bld.face(BrepPlane(origin: T[i], normal: pn), [
            loop([(e0, false), (e1, false), (e2, true)]),
          ]);
        }
      }
    }
    // Copy only untouched topology. All copied inner-cap edges are shared
    // with the attached faces through the same map, preserving their curves.
    final vertexMap = <int, int>{};
    final edgeMap = <int, int>{};
    int copyVertex(int id) =>
        vertexMap.putIfAbsent(id, () => bld.vertex(s.vertices[id]));
    int copyEdge(int id) => edgeMap.putIfAbsent(id, () {
      final e = s.edges[id];
      return bld.edge(
        copyVertex(e.v0),
        copyVertex(e.v1),
        e.curve,
        origin: e.origin,
      );
    });
    BrepLoop copyLoop(BrepLoop loop) => BrepLoop([
      for (final co in loop.coedges)
        BrepCoedge(copyEdge(co.edge), reversed: co.reversed),
    ]);
    bld.faces.first.loops.addAll(s.faces[opposite].loops.skip(1).map(copyLoop));
    for (var fi = 0; fi < s.faces.length; fi++) {
      if (rebuiltFaces.contains(fi)) continue;
      final f = s.faces[fi];
      bld.face(
        f.surface,
        f.loops.map(copyLoop).toList(),
        reversed: f.reversed,
        origin: f.origin,
      );
    }
    return bld.build();
  }

  /// Subtract equal-radius cylinders along a convex cap perimeter together.
  /// Adjacent cutters meet on their angle-bisector plane. That intersection
  /// is an ellipse from the inset cap corner to the lowered wall corner.
  static BrepSolid _inversePrismCap(
    BrepSolid s,
    int cap,
    int opposite,
    Set<int> rebuilt,
    List<Vector3> c,
    List<Vector3> b,
    List<Vector3> t,
    List<Vector3> w,
    Vector3 normal,
    double r,
  ) {
    final builder = BrepBuilder();
    final n = c.length;
    final bi = [for (final p in b) builder.vertex(p)];
    final ti = [for (final p in t) builder.vertex(p)];
    final wi = [for (final p in w) builder.vertex(p)];
    int line(int a, int z) =>
        builder.edge(a, z, BrepLine(builder.vertices[a], builder.vertices[z]));
    final bottom = [for (var i = 0; i < n; i++) line(bi[i], bi[(i + 1) % n])];
    final top = [for (var i = 0; i < n; i++) line(ti[i], ti[(i + 1) % n])];
    final wall = [for (var i = 0; i < n; i++) line(wi[i], wi[(i + 1) % n])];
    final upright = [for (var i = 0; i < n; i++) line(bi[i], wi[i])];
    final seams = <int>[];
    for (var i = 0; i < n; i++) {
      final bisector = t[i] - c[i];
      final x = bisector.normalized();
      seams.add(
        builder.edge(
          ti[i],
          wi[i],
          BrepEllipse(
            center: c[i],
            normal: x.cross(-normal),
            xDir: x,
            rMajor: bisector.length,
            rMinor: r,
            t0: 0,
            t1: math.pi / 2,
          ),
        ),
      );
    }
    BrepLoop loop(List<(int, bool)> edges) =>
        BrepLoop([for (final e in edges) BrepCoedge(e.$1, reversed: e.$2)]);
    builder.face(BrepPlane(origin: b.first, normal: -normal), [
      loop([for (final e in bottom.reversed) (e, true)]),
    ], origin: s.faces[opposite].origin);
    builder.face(BrepPlane(origin: t.first, normal: normal), [
      loop([for (final e in top) (e, false)]),
    ], origin: s.faces[cap].origin);
    for (var i = 0; i < n; i++) {
      final j = (i + 1) % n;
      final axis = (c[j] - c[i]).normalized();
      builder.face(BrepPlane(origin: c[i], normal: axis.cross(normal)), [
        loop([
          (bottom[i], false),
          (upright[j], false),
          (wall[i], true),
          (upright[i], true),
        ]),
      ]);
      builder.face(
        BrepCylinder(origin: c[i], axis: axis, radius: r, xDir: -normal),
        [
          loop([
            (top[i], true),
            (seams[i], false),
            (wall[i], false),
            (seams[j], true),
          ]),
        ],
        reversed: true,
      );
    }
    // Preserve inner loops on the opposite cap and any attached geometry.
    final vertices = <int, int>{};
    final edges = <int, int>{};
    int copyVertex(int id) =>
        vertices.putIfAbsent(id, () => builder.vertex(s.vertices[id]));
    int copyEdge(int id) => edges.putIfAbsent(id, () {
      final e = s.edges[id];
      return builder.edge(
        copyVertex(e.v0),
        copyVertex(e.v1),
        e.curve,
        origin: e.origin,
      );
    });
    BrepLoop copyLoop(BrepLoop l) => BrepLoop([
      for (final e in l.coedges)
        BrepCoedge(copyEdge(e.edge), reversed: e.reversed),
    ]);
    builder.faces.first.loops.addAll(
      s.faces[opposite].loops.skip(1).map(copyLoop),
    );
    for (var i = 0; i < s.faces.length; i++) {
      if (rebuilt.contains(i)) continue;
      final f = s.faces[i];
      builder.face(
        f.surface,
        f.loops.map(copyLoop).toList(),
        reversed: f.reversed,
        origin: f.origin,
      );
    }
    return builder.build();
  }

  /// Signed circle span from [from] to [to] (radians about [axis], xDir=[from]).
  /// Positive is CCW; negative is CW. Picks the sector that contains [via]
  /// so a 90° box corner is −π/2, not a 270° trip through the solid.
  static double _arcToward(
    Vector3 from,
    Vector3 to,
    Vector3 axis,
    Vector3 via,
  ) {
    final f = (from - axis * from.dot(axis));
    final t = (to - axis * to.dot(axis));
    final w = (via - axis * via.dot(axis));
    if (f.length < 1e-12 || t.length < 1e-12) return math.pi / 2;
    f.normalize();
    t.normalize();
    var ang = math.atan2(f.cross(t).dot(axis), f.dot(t));
    if (ang < 0) ang += math.pi * 2;
    if (w.length < 1e-12) {
      if (ang > math.pi) ang -= math.pi * 2;
      return ang.abs() < 1e-6 ? math.pi / 2 : ang;
    }
    w.normalize();
    var viaAng = math.atan2(f.cross(w).dot(axis), f.dot(w));
    if (viaAng < 0) viaAng += math.pi * 2;
    final onCcw = viaAng <= ang + 1e-6;
    final signed = onCcw ? ang : ang - math.pi * 2;
    return signed.abs() < 1e-6 ? math.pi / 2 : signed;
  }

  /// Two opposite +Z edges (front y=0 and back y=sy).
  static BrepSolid _canonicalOppositeBlend(
    double sx,
    double sy,
    double sz,
    double r, {
    required bool round,
  }) {
    var s = BrepPrimitives.box(sx, sy, sz);
    // Top-front = edge 4, top-back = edge 6 (see BrepPrimitives.box).
    if (round) {
      s = _filletOne(s, 4, r);
      final id = _matchCurve(
        s,
        BrepLine(Vector3(sx, sy, sz), Vector3(0, sy, sz)),
        r,
      );
      if (id >= 0) s = _filletOne(s, id, r);
    } else {
      s = _chamferOne(s, 4, r);
      final id = _matchCurve(
        s,
        BrepLine(Vector3(sx, sy, sz), Vector3(0, sy, sz)),
        r,
      );
      if (id >= 0) s = _chamferOne(s, id, r);
    }
    return s;
  }

  /// Front (y=0) and left (x=0) of the +Z face, with a spherical/chamfer corner.
  static BrepSolid _canonicalAdjacentFrontLeft(
    double sx,
    double sy,
    double sz,
    double r, {
    required bool round,
  }) {
    final b = BrepBuilder();
    int v(double x, double y, double z) => b.vertex(Vector3(x, y, z));
    final b0 = v(0, 0, 0);
    final b1 = v(sx, 0, 0);
    final b2 = v(sx, sy, 0);
    final b3 = v(0, sy, 0);
    final sFL = v(0, 0, sz - r);
    final fFR = v(sx, 0, sz - r);
    final tBR = v(sx, sy, sz);
    final lBL = v(0, sy, sz - r);
    final tFL = v(r, r, sz);
    final tFR = v(sx, r, sz);
    final tBL = v(r, sy, sz);
    final fFL = v(r, 0, sz - r);
    final lFL = v(0, r, sz - r);
    final vBR = tBR;

    final lineOf = <String, int>{};
    int line(int i, int j) {
      final lo = i < j ? i : j;
      final hi = i < j ? j : i;
      final k = '$lo-$hi';
      final existing = lineOf[k];
      if (existing != null) return existing;
      final id = b.edge(i, j, BrepLine(b.vertices[i], b.vertices[j]));
      lineOf[k] = id;
      return id;
    }

    int arc(int i, int j, Vector3 center, Vector3 normal, Vector3 xDir) {
      if (!round) return line(i, j);
      return b.edge(
        i,
        j,
        BrepCircle(
          center: center,
          normal: normal,
          radius: r,
          xDir: xDir,
          t0: 0,
          t1: math.pi / 2,
        ),
      );
    }

    final eBot01 = line(b0, b1);
    final eBot12 = line(b1, b2);
    final eBot23 = line(b2, b3);
    final eBot30 = line(b3, b0);
    final eV0 = line(b0, sFL);
    final eV1 = line(b1, fFR);
    final eV2 = line(b2, vBR);
    final eV3 = line(b3, lBL);
    final eTopFLFR = line(tFL, tFR);
    final eTopFRBR = line(tFR, tBR);
    final eTopBRBL = line(tBR, tBL);
    final eTopBLFL = line(tBL, tFL);
    final eFrontOff = line(fFL, fFR);
    final eLeftOff = line(lFL, lBL);
    final eShelfF = line(sFL, fFL);
    final eShelfL = line(sFL, lFL);

    final cFL = Vector3(r, r, sz - r);
    final cFR = Vector3(sx, r, sz - r);
    final cBL = Vector3(r, sy, sz - r);
    final aFrontL = arc(tFL, fFL, cFL, Vector3(1, 0, 0), Vector3(0, 0, 1));
    final aFrontR = arc(tFR, fFR, cFR, Vector3(1, 0, 0), Vector3(0, 0, 1));
    final aLeftF = arc(tFL, lFL, cFL, Vector3(0, -1, 0), Vector3(0, 0, 1));
    final aLeftB = arc(tBL, lBL, cBL, Vector3(0, -1, 0), Vector3(0, 0, 1));

    BrepLoop loop(List<(int, bool)> es) =>
        BrepLoop([for (final e in es) BrepCoedge(e.$1, reversed: e.$2)]);

    b.face(BrepPlane(origin: Vector3(0, 0, 0), normal: Vector3(0, 0, -1)), [
      loop([(eBot01, true), (eBot30, true), (eBot23, true), (eBot12, true)]),
    ]);
    b.face(BrepPlane(origin: Vector3(0, 0, sz), normal: Vector3(0, 0, 1)), [
      loop([
        (eTopFLFR, false),
        (eTopFRBR, false),
        (eTopBRBL, false),
        (eTopBLFL, false),
      ]),
    ]);
    b.face(BrepPlane(origin: Vector3(0, 0, 0), normal: Vector3(0, -1, 0)), [
      loop([
        (eBot01, false),
        (eV1, false),
        (eFrontOff, true),
        (eShelfF, true),
        (eV0, true),
      ]),
    ]);
    b.face(BrepPlane(origin: Vector3(sx, 0, 0), normal: Vector3(1, 0, 0)), [
      loop([
        (eBot12, false),
        (eV2, false),
        (eTopFRBR, true),
        (aFrontR, false),
        (eV1, true),
      ]),
    ]);
    b.face(BrepPlane(origin: Vector3(0, sy, 0), normal: Vector3(0, 1, 0)), [
      loop([
        (eBot23, false),
        (eV3, false),
        (aLeftB, true),
        (eTopBRBL, true),
        (eV2, true),
      ]),
    ]);
    b.face(BrepPlane(origin: Vector3(0, 0, 0), normal: Vector3(-1, 0, 0)), [
      loop([
        (eBot30, false),
        (eV0, false),
        (eShelfL, false),
        (eLeftOff, false),
        (eV3, true),
      ]),
    ]);

    BrepSurface frontSurf() {
      if (round) {
        return BrepCylinder(
          origin: Vector3(0, r, sz - r),
          axis: Vector3(1, 0, 0),
          radius: r,
        );
      }
      return BrepPlane(
        origin: Vector3(0, 0, sz - r),
        normal: Vector3(0, -1, 1).normalized(),
      );
    }

    BrepSurface leftSurf() {
      if (round) {
        return BrepCylinder(
          origin: Vector3(r, 0, sz - r),
          axis: Vector3(0, 1, 0),
          radius: r,
        );
      }
      return BrepPlane(
        origin: Vector3(0, 0, sz - r),
        normal: Vector3(-1, 0, 1).normalized(),
      );
    }

    b.face(frontSurf(), [
      loop([
        (eTopFLFR, false),
        (aFrontR, false),
        (eFrontOff, true),
        (aFrontL, true),
      ]),
    ]);
    b.face(leftSurf(), [
      loop([
        (eTopBLFL, false),
        (aLeftF, false),
        (eLeftOff, false),
        (aLeftB, true),
      ]),
    ]);

    if (round) {
      final hFL = b.edge(
        lFL,
        fFL,
        BrepCircle(
          center: cFL,
          normal: Vector3(0, 0, 1),
          radius: r,
          xDir: Vector3(1, 0, 0),
          t0: math.pi,
          t1: 1.5 * math.pi,
        ),
      );
      b.face(BrepSphere(center: cFL, radius: r), [
        loop([(aFrontL, false), (hFL, true), (aLeftF, true)]),
      ]);
      b.face(
        BrepPlane(origin: Vector3(0, 0, sz - r), normal: Vector3(0, 0, 1)),
        [
          loop([(eShelfF, false), (hFL, true), (eShelfL, true)]),
        ],
      );
    } else {
      final eFL = line(fFL, lFL);
      final eTFL = line(tFL, fFL);
      final eTL = line(tFL, lFL);
      var n = (b.vertices[fFL] - b.vertices[tFL]).cross(
        b.vertices[lFL] - b.vertices[tFL],
      );
      if (n.length > 1e-12) {
        n.normalize();
        b.face(BrepPlane(origin: b.vertices[tFL], normal: n), [
          loop([(eTFL, false), (eFL, false), (eTL, true)]),
        ]);
      }
      b.face(
        BrepPlane(origin: Vector3(0, 0, sz - r), normal: Vector3(0, 0, 1)),
        [
          loop([(eShelfF, false), (eFL, false), (eShelfL, true)]),
        ],
      );
    }

    return b.build();
  }

  /// Unit-origin box [0,sx]×[0,sy]×[0,sz] with the +Z face loop blended.
  static BrepSolid _canonicalFaceBlend(
    double sx,
    double sy,
    double sz,
    double r, {
    required bool round,
  }) {
    final b = BrepBuilder();
    int v(double x, double y, double z) => b.vertex(Vector3(x, y, z));
    final b0 = v(0, 0, 0);
    final b1 = v(sx, 0, 0);
    final b2 = v(sx, sy, 0);
    final b3 = v(0, sy, 0);
    final s0 = v(0, 0, sz - r);
    final s1 = v(sx, 0, sz - r);
    final s2 = v(sx, sy, sz - r);
    final s3 = v(0, sy, sz - r);
    final t0 = v(r, r, sz);
    final t1 = v(sx - r, r, sz);
    final t2 = v(sx - r, sy - r, sz);
    final t3 = v(r, sy - r, sz);
    final f0 = v(r, 0, sz - r);
    final f1 = v(sx - r, 0, sz - r);
    final k0 = v(sx - r, sy, sz - r);
    final k1 = v(r, sy, sz - r);
    final l0 = v(0, r, sz - r);
    final l1 = v(0, sy - r, sz - r);
    final p0 = v(sx, r, sz - r);
    final p1 = v(sx, sy - r, sz - r);

    (int, bool) sew(int i, int j, BrepCurve curve) {
      // Unique by endpoints; rebuild is from scratch so no duplicates yet.
      return (b.edge(i, j, curve), false);
    }

    final lineOf = <String, int>{};
    int line(int i, int j) {
      final lo = i < j ? i : j;
      final hi = i < j ? j : i;
      final k = '$lo-$hi';
      final existing = lineOf[k];
      if (existing != null) return existing;
      final id = b.edge(i, j, BrepLine(b.vertices[i], b.vertices[j]));
      lineOf[k] = id;
      return id;
    }

    int arc(
      int i,
      int j,
      Vector3 center,
      Vector3 normal,
      Vector3 xDir,
      double t0,
      double t1,
    ) {
      if (!round) {
        return line(i, j);
      }
      return sew(
        i,
        j,
        BrepCircle(
          center: center,
          normal: normal,
          radius: r,
          xDir: xDir,
          t0: t0,
          t1: t1,
        ),
      ).$1;
    }

    final eBot01 = line(b0, b1);
    final eBot12 = line(b1, b2);
    final eBot23 = line(b2, b3);
    final eBot30 = line(b3, b0);
    final eV0 = line(b0, s0);
    final eV1 = line(b1, s1);
    final eV2 = line(b2, s2);
    final eV3 = line(b3, s3);
    final eTop01 = line(t0, t1);
    final eTop12 = line(t1, t2);
    final eTop23 = line(t2, t3);
    final eTop30 = line(t3, t0);

    final eF01 = line(f0, f1);
    final eK01 = line(k0, k1);
    final eL01 = line(l0, l1);
    final eP01 = line(p0, p1);
    final eS0F0 = line(s0, f0);
    final eF1S1 = line(f1, s1);
    final eS1P0 = line(s1, p0);
    final eP1S2 = line(p1, s2);
    final eS2K0 = line(s2, k0);
    final eK1S3 = line(k1, s3);
    final eS3L1 = line(s3, l1);
    final eL0S0 = line(l0, s0);

    // Quarter-circles / chamfer lines at cylinder/plane ends.
    final cFL = Vector3(r, r, sz - r);
    final cFR = Vector3(sx - r, r, sz - r);
    final cBR = Vector3(sx - r, sy - r, sz - r);
    final cBL = Vector3(r, sy - r, sz - r);
    final aF0 = arc(
      t0,
      f0,
      cFL,
      Vector3(1, 0, 0),
      Vector3(0, 0, 1),
      0,
      math.pi / 2,
    );
    final aF1 = arc(
      t1,
      f1,
      cFR,
      Vector3(1, 0, 0),
      Vector3(0, 0, 1),
      0,
      math.pi / 2,
    );
    final aK0 = arc(
      t2,
      k0,
      cBR,
      Vector3(-1, 0, 0),
      Vector3(0, 0, 1),
      0,
      math.pi / 2,
    );
    final aK1 = arc(
      t3,
      k1,
      cBL,
      Vector3(-1, 0, 0),
      Vector3(0, 0, 1),
      0,
      math.pi / 2,
    );
    final aL0 = arc(
      t0,
      l0,
      cFL,
      Vector3(0, -1, 0),
      Vector3(0, 0, 1),
      0,
      math.pi / 2,
    );
    final aL1 = arc(
      t3,
      l1,
      cBL,
      Vector3(0, -1, 0),
      Vector3(0, 0, 1),
      0,
      math.pi / 2,
    );
    final aP0 = arc(
      t1,
      p0,
      cFR,
      Vector3(0, 1, 0),
      Vector3(0, 0, 1),
      0,
      math.pi / 2,
    );
    final aP1 = arc(
      t2,
      p1,
      cBR,
      Vector3(0, 1, 0),
      Vector3(0, 0, 1),
      0,
      math.pi / 2,
    );

    BrepLoop loop(List<(int, bool)> es) =>
        BrepLoop([for (final e in es) BrepCoedge(e.$1, reversed: e.$2)]);

    // Bottom.
    b.face(BrepPlane(origin: Vector3(0, 0, 0), normal: Vector3(0, 0, -1)), [
      loop([(eBot01, true), (eBot30, true), (eBot23, true), (eBot12, true)]),
    ]);
    // Top inset.
    b.face(BrepPlane(origin: Vector3(0, 0, sz), normal: Vector3(0, 0, 1)), [
      loop([
        (eTop01, false),
        (eTop12, false),
        (eTop23, false),
        (eTop30, false),
      ]),
    ]);
    // Front y=0.
    b.face(BrepPlane(origin: Vector3(0, 0, 0), normal: Vector3(0, -1, 0)), [
      loop([
        (eBot01, false),
        (eV1, false),
        (eF1S1, true),
        (eF01, true),
        (eS0F0, true),
        (eV0, true),
      ]),
    ]);
    // Back y=sy.
    b.face(BrepPlane(origin: Vector3(0, sy, 0), normal: Vector3(0, 1, 0)), [
      loop([
        (eBot23, false),
        (eV3, false),
        (eK1S3, true),
        (eK01, true),
        (eS2K0, true),
        (eV2, true),
      ]),
    ]);
    // Left x=0.
    b.face(BrepPlane(origin: Vector3(0, 0, 0), normal: Vector3(-1, 0, 0)), [
      loop([
        (eBot30, false),
        (eV0, false),
        (eL0S0, true),
        (eL01, false),
        (eS3L1, true),
        (eV3, true),
      ]),
    ]);
    // Right x=sx.
    b.face(BrepPlane(origin: Vector3(sx, 0, 0), normal: Vector3(1, 0, 0)), [
      loop([
        (eBot12, false),
        (eV2, false),
        (eP1S2, true),
        (eP01, true),
        (eS1P0, true),
        (eV1, true),
      ]),
    ]);

    BrepSurface blendSurf(Vector3 origin, Vector3 axis) {
      if (round) {
        return BrepCylinder(origin: origin, axis: axis, radius: r);
      }
      if (axis.x.abs() > 0.5) {
        // Front/back chamfer in YZ.
        final towardFront = origin.y < sy / 2;
        final n = Vector3(0, towardFront ? -1.0 : 1.0, 1).normalized();
        return BrepPlane(origin: origin, normal: n);
      }
      final towardLeft = origin.x < sx / 2;
      final n = Vector3(towardLeft ? -1.0 : 1.0, 0, 1).normalized();
      return BrepPlane(origin: origin, normal: n);
    }

    // Front blend: t0-t1-f1-f0
    b.face(blendSurf(Vector3(0, r, sz - r), Vector3(1, 0, 0)), [
      loop([(eTop01, false), (aF1, false), (eF01, true), (aF0, true)]),
    ]);
    // Back blend: t2-t3-k1-k0 — top edge t2-t3 is eTop23
    b.face(blendSurf(Vector3(0, sy - r, sz - r), Vector3(1, 0, 0)), [
      loop([(eTop23, false), (aK1, false), (eK01, true), (aK0, true)]),
    ]);
    // Left blend: t3-t0-l0-l1
    b.face(blendSurf(Vector3(r, 0, sz - r), Vector3(0, 1, 0)), [
      loop([(eTop30, false), (aL0, false), (eL01, false), (aL1, true)]),
    ]);
    // Right blend: t1-t2-p1-p0
    b.face(blendSurf(Vector3(sx - r, 0, sz - r), Vector3(0, 1, 0)), [
      loop([(eTop12, false), (aP1, false), (eP01, true), (aP0, true)]),
    ]);

    if (round) {
      int sphArc(
        int i,
        int j,
        Vector3 c,
        Vector3 n,
        Vector3 xDir,
        double t0,
        double t1,
      ) => sew(
        i,
        j,
        BrepCircle(center: c, normal: n, radius: r, xDir: xDir, t0: t0, t1: t1),
      ).$1;

      // 90° arcs in z=sz-r toward each vertical corner.
      final hFL = sphArc(
        l0,
        f0,
        cFL,
        Vector3(0, 0, 1),
        Vector3(1, 0, 0),
        math.pi,
        1.5 * math.pi,
      );
      final hFR = sphArc(
        f1,
        p0,
        cFR,
        Vector3(0, 0, 1),
        Vector3(1, 0, 0),
        1.5 * math.pi,
        math.pi * 2,
      );
      final hBR = sphArc(
        p1,
        k0,
        cBR,
        Vector3(0, 0, 1),
        Vector3(1, 0, 0),
        0,
        math.pi / 2,
      );
      final hBL = sphArc(
        k1,
        l1,
        cBL,
        Vector3(0, 0, 1),
        Vector3(1, 0, 0),
        math.pi / 2,
        math.pi,
      );

      void sphereCorner(
        Vector3 c,
        int eA,
        bool rA,
        int eB,
        bool rB,
        int eC,
        bool rC,
      ) {
        b.face(BrepSphere(center: c, radius: r), [
          loop([(eA, rA), (eB, rB), (eC, rC)]),
        ]);
      }

      // t0 -> f0 -> l0 -> t0
      sphereCorner(cFL, aF0, false, hFL, true, aL0, true);
      // t1 -> p0 -> f1 -> t1
      sphereCorner(cFR, aP0, false, hFR, true, aF1, true);
      // t2 -> k0 -> p1 -> t2
      sphereCorner(cBR, aK0, false, hBR, true, aP1, true);
      // t3 -> l1 -> k1 -> t3
      sphereCorner(cBL, aL1, false, hBL, true, aK1, true);

      void shelf(
        int eA,
        bool rA,
        int eB,
        bool rB,
        int eC,
        bool rC,
        Vector3 origin,
      ) {
        b.face(BrepPlane(origin: origin, normal: Vector3(0, 0, 1)), [
          loop([(eA, rA), (eB, rB), (eC, rC)]),
        ]);
      }

      shelf(eS0F0, false, hFL, true, eL0S0, false, b.vertices[s0]);
      shelf(eS1P0, false, hFR, true, eF1S1, false, b.vertices[s1]);
      shelf(eS2K0, false, hBR, true, eP1S2, false, b.vertices[s2]);
      shelf(eS3L1, false, hBL, true, eK1S3, false, b.vertices[s3]);
    } else {
      // Chamfer corners: triangle t-endA-endB.
      void tri(int i, int j, int k) {
        final e0 = line(i, j);
        final e1 = line(j, k);
        final e2 = line(k, i);
        final p = b.vertices[i];
        final q = b.vertices[j];
        final w = b.vertices[k];
        var n = (q - p).cross(w - p);
        if (n.length < 1e-12) return;
        n.normalize();
        b.face(BrepPlane(origin: p, normal: n), [
          loop([(e0, false), (e1, false), (e2, true)]),
        ]);
      }

      tri(t0, f0, l0);
      tri(t1, p0, f1);
      tri(t2, k0, p1);
      tri(t3, l1, k1);

      void shelf(int c, int a, int d) {
        final e0 = line(c, a);
        final e1 = line(a, d);
        final e2 = line(d, c);
        b.face(BrepPlane(origin: b.vertices[c], normal: Vector3(0, 0, 1)), [
          loop([(e0, false), (e1, false), (e2, true)]),
        ]);
      }

      shelf(s0, f0, l0);
      shelf(s1, p0, f1);
      shelf(s2, k0, p1);
      shelf(s3, l1, k1);
    }

    return b.build();
  }

  static int _matchCurve(BrepSolid s, BrepCurve orig, double inset) {
    final mid = orig.evaluate(0.5 * (orig.t0 + orig.t1));
    var best = -1;
    var bestD = double.infinity;
    for (var i = 0; i < s.edges.length; i++) {
      final c = s.edges[i].curve;
      if (c.runtimeType != orig.runtimeType) continue;
      final m = c.evaluate(0.5 * (c.t0 + c.t1));
      final d = (m - mid).length;
      if (d < bestD) {
        bestD = d;
        best = i;
      }
    }
    final maxD = math.max(1e-3, inset.abs() * 0.65);
    if (best < 0 || bestD > maxD) return -1;
    return best;
  }

  static bool _canBlend(BrepSolid s, int edge) {
    final adj = facesForEdge(s, edge);
    if (adj.length < 2) return false;
    final s1 = s.faces[adj[0]].surface;
    final s2 = s.faces[adj[1]].surface;
    if (s1 is BrepPlane && s2 is BrepPlane) return true;
    return s.edges[edge].curve is BrepCircle &&
        ((s1 is BrepPlane && s2 is BrepCylinder) ||
            (s1 is BrepCylinder && s2 is BrepPlane));
  }

  /// Adds the missing feature edge at any vertex where two of three
  /// incident creases are already selected. A cube face loop therefore
  /// picks up the four uprights (eight edges) and leaves the opposite
  /// face sharp.
  static List<int> completeCornerEdges(BrepSolid solid, List<int> edgeIds) {
    final selected = <int>{
      for (final id in edgeIds)
        if (id >= 0 && id < solid.edges.length) id,
    };
    if (selected.isEmpty) return const [];

    final atVertex = <int, List<int>>{};
    final adjacency = _blendAdjacency(solid);
    for (var id = 0; id < solid.edges.length; id++) {
      if (adjacency[id]?.length != 2) continue;
      final e = solid.edges[id];
      if (e.v0 == e.v1) continue;
      final faces = adjacency[id]!;
      final mid = e.curve.evaluate((e.curve.t0 + e.curve.t1) * .5);
      final normals = [
        for (final fi in faces)
          (() {
            final f = solid.faces[fi], uv = f.surface.parameters(mid);
            return uv == null ? Vector3.zero() : f.outwardNormal(uv.$1, uv.$2);
          })(),
      ];
      if (normals[0].dot(normals[1]) > 1 - 1e-8) continue;
      (atVertex[e.v0] ??= []).add(id);
      (atVertex[e.v1] ??= []).add(id);
    }

    var changed = true;
    while (changed) {
      changed = false;
      for (final incident in atVertex.values) {
        if (incident.length != 3) continue;
        var nSel = 0;
        var missing = 0;
        for (var i = 0; i < 3; i++) {
          if (selected.contains(incident[i])) {
            nSel++;
          } else {
            missing = i;
          }
        }
        if (nSel != 2) continue;
        selected.add(incident[missing]);
        changed = true;
      }
    }
    return selected.toList()..sort();
  }

  static List<int> resolveEdges(
    BrepSolid solid,
    List<(Vector3, Vector3)> chords,
  ) {
    final ids = <int>{};
    for (final (a, b) in chords) {
      var best = -1;
      var bestD = double.infinity;
      for (var i = 0; i < solid.edges.length; i++) {
        if (facesForEdge(solid, i).length != 2) continue;
        final c = solid.edges[i].curve;
        final d =
            (c.closestPoint(a) - a).length + (c.closestPoint(b) - b).length;
        if (d < bestD) {
          bestD = d;
          best = i;
        }
      }
      if (best >= 0 && bestD < 2.0) ids.add(best);
    }
    return ids.toList();
  }

  /// Convex plane–plane inset equals radius; use clearance to other loops.
  static double maxConvexRadius(BrepSolid s, List<int> edgeIds) {
    var limit = double.infinity;
    for (final id in edgeIds) {
      if (id < 0 || id >= s.edges.length) continue;
      final edgeLimit = _edgeConvexLimit(s, id);
      if (edgeLimit < limit) limit = edgeLimit;
    }
    if (!limit.isFinite || limit <= 0) return 0;
    return limit;
  }

  /// BREP-based search bound. Detailed offset checks still run on the
  /// requested blend; display triangle subdivisions never limit this value.
  static double? maxBlendDistance(BrepSolid s, List<int> ids) {
    if (ids.isEmpty) return null;
    final transformed = BrepXform.unwrapAffine(s);
    if (transformed != null) {
      final scale = _similarityScale(transformed.transform);
      if (scale != null) {
        final limit = maxBlendDistance(transformed.solid, ids);
        return limit == null ? null : limit * scale;
      }
    }
    s = _blendPlanarRegions(
      BrepHeal.mergeFilledPlanarLoops(BrepXform.materializeForEditing(s)),
    );
    ids = _creaseEdges(s, ids);
    if (ids.isEmpty) return null;
    final planar = maxPlanarRadius(s, ids);
    if (planar != null) return planar;
    var limit = double.infinity;
    for (final id in ids) {
      if (id < 0 || id >= s.edges.length) return null;
      limit = math.min(limit, _edgeConvexLimit(s, id));
      final curve = s.edges[id].curve;
      if (curve is BrepCircle) {
        for (final fi in facesForEdge(s, id)) {
          if (s.faces[fi].surface is BrepCylinder) {
            final frame = chamferFrame(s, id);
            if (frame != null && frame.$2.dot(frame.$1 - curve.center) < 0) {
              limit = math.min(limit, curve.radius);
            }
          }
        }
      }
    }
    if (!limit.isFinite) {
      final points = [
        for (final e in s.edges)
          for (var i = 0; i <= 16; i++)
            e.curve.evaluate(e.curve.t0 + (e.curve.t1 - e.curve.t0) * i / 16),
      ];
      if (points.isEmpty) return null;
      final lo = points.first.clone(), hi = lo.clone();
      for (final p in points) {
        for (var k = 0; k < 3; k++) {
          lo[k] = math.min(lo[k], p[k]);
          hi[k] = math.max(hi[k], p[k]);
        }
      }
      limit = (hi - lo).length * .5;
    }
    return limit.isFinite && limit > 0 ? limit : null;
  }

  /// Search bound for planar chamfer setbacks. Short edges left by an earlier
  /// bevel do not constrain a later bevel: measure travel across each face.
  /// The operation still checks each requested contact and the resulting solid.
  static double? maxChamferDistance(BrepSolid s, List<int> ids) {
    final planar = BrepXform.materializeForEditing(s);
    if (ids.isEmpty || !BrepConvexPlanes.supports(planar)) {
      return maxBlendDistance(s, ids);
    }
    var limit = double.infinity;
    for (final id in ids.toSet()) {
      final frame = chamferFrame(planar, id);
      if (frame == null) return maxBlendDistance(s, ids);
      for (final direction in [frame.$2, frame.$3]) {
        for (final face in planar.faces) {
          final n = face.outwardNormal(0, 0);
          final rate = n.dot(direction);
          if (rate <= 1e-10) continue;
          final remaining = n.dot(
            (face.surface as BrepPlane).origin - frame.$1,
          );
          limit = math.min(limit, remaining / rate);
        }
      }
    }
    return limit.isFinite && limit > 0 ? limit : maxBlendDistance(s, ids);
  }

  /// Local radius bound for straight edges between planar faces. Returns null
  /// when this calculation does not describe the selected geometry. Face trims,
  /// rather than display triangles, determine how far tangency can travel.
  static double? maxPlanarRadius(BrepSolid s, List<int> edgeIds) {
    if (edgeIds.isEmpty) return null;
    final selected = edgeIds.toSet();
    var limit = double.infinity;
    for (final id in selected) {
      if (id < 0 || id >= s.edges.length || s.edges[id].curve is! BrepLine) {
        return null;
      }
      final adj = facesForEdge(s, id);
      if (adj.length != 2 ||
          adj.any((fi) => s.faces[fi].surface is! BrepPlane)) {
        return null;
      }
      Vector3 normal(int fi) {
        final f = s.faces[fi];
        final n = (f.surface as BrepPlane).normal.normalized();
        return f.reversed ? -n : n;
      }

      final phi = math.acos(
        normal(adj[0]).dot(normal(adj[1])).clamp(-1.0, 1.0),
      );
      if (phi < 1e-3 || phi > math.pi - 1e-3) return null;
      final factor = math.tan((math.pi - phi) / 2);
      for (final fi in adj) {
        final setback = _faceInsetFromEdge(s, fi, id, selected: selected);
        limit = math.min(limit, setback * factor);
      }
    }
    return limit.isFinite && limit > 0 ? limit : null;
  }

  static double _edgeConvexLimit(BrepSolid s, int edgeId) {
    final adj = facesForEdge(s, edgeId);
    if (adj.length < 2) return 0;
    var limit = double.infinity;
    for (final fi in adj) {
      final d = _faceInsetFromEdge(s, fi, edgeId);
      if (d < limit) limit = d;
    }
    return limit;
  }

  static double _faceInsetFromEdge(
    BrepSolid s,
    int fi,
    int edgeId, {
    Set<int> selected = const {},
  }) {
    final e = s.edges[edgeId];
    final a = s.vertices[e.v0];
    final b = s.vertices[e.v1];
    var best = double.infinity;
    for (final loop in s.faces[fi].loops) {
      for (final c in loop.coedges) {
        if (c.edge == edgeId) continue;
        final o = s.edges[c.edge];
        final va = s.vertices[o.v0];
        final vb = s.vertices[o.v1];
        if (_samePt(a, va) ||
            _samePt(a, vb) ||
            _samePt(b, va) ||
            _samePt(b, vb)) {
          continue;
        }
        final d =
            _curveClearance(e.curve, o.curve) *
            (selected.contains(c.edge) ? 0.5 : 1.0);
        if (d > 1e-6 && d < best) best = d;
      }
    }
    return best;
  }

  static bool _samePt(Vector3 p, Vector3 q) => (p - q).length2 < 1e-16;

  static double _curveClearance(BrepCurve a, BrepCurve b) {
    var best = double.infinity;
    void probe(BrepCurve src, BrepCurve dst) {
      for (var i = 0; i <= 8; i++) {
        final t = src.t0 + (src.t1 - src.t0) * (i / 8);
        final p = src.evaluate(t);
        final d = (dst.closestPoint(p) - p).length;
        if (d < best) best = d;
      }
    }

    probe(a, b);
    probe(b, a);
    return best;
  }

  static List<int> facesForEdge(BrepSolid s, int edge) {
    final out = <int>[];
    for (var f = 0; f < s.faces.length; f++) {
      for (final loop in s.faces[f].loops) {
        if (loop.coedges.any((c) => c.edge == edge)) out.add(f);
      }
    }
    return out;
  }

  /// Stable first/second setback directions in the edge's normal section.
  /// Circular cap rims use the cap first and the axial wall second.
  static (Vector3, Vector3, Vector3)? chamferFrame(BrepSolid s, int edgeId) {
    if (edgeId < 0 || edgeId >= s.edges.length) return null;
    final transformed = BrepXform.unwrapAffine(s);
    if (transformed != null &&
        _similarityScale(transformed.transform) != null) {
      final local = chamferFrame(transformed.solid, edgeId);
      if (local == null) return null;
      final m = transformed.transform, origin = m.point(local.$1);
      return (
        origin,
        (m.point(local.$1 + local.$2) - origin).normalized(),
        (m.point(local.$1 + local.$3) - origin).normalized(),
      );
    }
    s = BrepXform.materializeForEditing(s);
    final edge = s.edges[edgeId];
    final adj = facesForEdge(s, edgeId);
    if (adj.length != 2) return null;
    final p = edge.curve.evaluate((edge.curve.t0 + edge.curve.t1) / 2);
    final t = edge.curve
        .tangent((edge.curve.t0 + edge.curve.t1) / 2)
        .normalized();
    Vector3 into(int fi) {
      final f = s.faces[fi];
      final uv = f.surface.parameters(p);
      if (uv == null)
        throw BrepUnsupported('Cannot locate chamfer on its face');
      final n = f.outwardNormal(uv.$1, uv.$2);
      final co = f.loops
          .expand((l) => l.coedges)
          .firstWhere((c) => c.edge == edgeId);
      return n.cross(co.reversed ? -t : t).normalized();
    }

    if (edge.curve is BrepCircle) {
      final plane = adj
          .where((fi) => s.faces[fi].surface is BrepPlane)
          .firstOrNull;
      final cylinder = adj
          .where((fi) => s.faces[fi].surface is BrepCylinder)
          .firstOrNull;
      if (plane != null && cylinder != null) {
        final cap = s.faces[plane], wall = s.faces[cylinder];
        final uv = wall.surface.parameters(p)!;
        final radial = wall.outwardNormal(uv.$1, uv.$2);
        final step = math.max(1e-5, (edge.curve as BrepCircle).radius * 1e-5);
        final minus = BrepClassify.containsPlaneFace(s, cap, p - radial * step);
        final plus = BrepClassify.containsPlaneFace(s, cap, p + radial * step);
        if (minus != plus && minus != null && plus != null) {
          final sign = minus ? -1.0 : 1.0;
          return (p, radial * sign, cap.outwardNormal(0, 0) * sign);
        }
      }
    }
    var first = adj[0], second = adj[1];
    if (edge.curve is BrepCircle && s.faces[second].surface is BrepPlane) {
      first = adj[1];
      second = adj[0];
    }
    return (p, into(first), into(second));
  }

  static BrepSolid _filletOne(
    BrepSolid s,
    int edgeId,
    double r, {
    bool inverse = false,
  }) {
    final e = s.edges[edgeId];
    final adj = facesForEdge(s, edgeId);
    if (adj.length < 2) {
      throw BrepUnsupported('Fillet edge has fewer than two faces');
    }
    final f1 = s.faces[adj[0]];
    final f2 = s.faces[adj[1]];
    if (f1.surface is BrepPlane && f2.surface is BrepPlane) {
      return _filletPlanePlane(s, edgeId, adj[0], adj[1], r, inverse: inverse);
    }
    final plane = f1.surface is BrepPlane
        ? f1
        : f2.surface is BrepPlane
        ? f2
        : null;
    final cyl = f1.surface is BrepCylinder
        ? f1
        : f2.surface is BrepCylinder
        ? f2
        : null;
    if (plane != null && cyl != null && e.curve is BrepCircle) {
      if (!(e.curve as BrepCircle).isClosed && !inverse) {
        final clipped = _blendClippedCylinder(s, edgeId, r, round: true);
        if (clipped != null) return clipped;
      }
      return _filletPlaneCylinder(s, edgeId, plane, cyl, r, inverse: inverse);
    }
    final cone = f1.surface is BrepCone
        ? f1
        : f2.surface is BrepCone
        ? f2
        : null;
    if (plane != null && cone != null && e.curve is BrepCircle && !inverse) {
      return _blendPlaneCone(s, edgeId, plane, cone, r, round: true);
    }
    throw BrepUnsupported(
      'Fillet of ${f1.surface.runtimeType}–${f2.surface.runtimeType} '
      'is not implemented yet',
    );
  }

  static BrepSolid _chamferOne(
    BrepSolid s,
    int edgeId,
    double d, {
    ChamferOptions options = const ChamferOptions(),
  }) {
    final e = s.edges[edgeId];
    final adj = facesForEdge(s, edgeId);
    if (adj.length < 2) {
      throw BrepUnsupported('Chamfer edge has fewer than two faces');
    }
    final f1 = s.faces[adj[0]];
    final f2 = s.faces[adj[1]];
    if (f1.surface is BrepPlane && f2.surface is BrepPlane) {
      return _chamferPlanePlane(s, edgeId, adj[0], adj[1], d, options);
    }
    final plane = f1.surface is BrepPlane
        ? f1
        : f2.surface is BrepPlane
        ? f2
        : null;
    final cyl = f1.surface is BrepCylinder
        ? f1
        : f2.surface is BrepCylinder
        ? f2
        : null;
    if (plane != null && cyl != null && e.curve is BrepCircle) {
      if (!(e.curve as BrepCircle).isClosed) {
        final clipped = _blendClippedCylinder(
          s,
          edgeId,
          d,
          round: false,
          options: options,
        );
        if (clipped != null) return clipped;
      }
      return _chamferPlaneCylinder(s, edgeId, plane, cyl, d, options);
    }
    final cone = f1.surface is BrepCone
        ? f1
        : f2.surface is BrepCone
        ? f2
        : null;
    if (plane != null && cone != null && e.curve is BrepCircle) {
      return _blendPlaneCone(
        s,
        edgeId,
        plane,
        cone,
        d,
        round: false,
        options: options,
      );
    }
    throw BrepUnsupported(
      'Chamfer of ${f1.surface.runtimeType}–${f2.surface.runtimeType} '
      'is not implemented yet',
    );
  }

  static Vector3? _centerOffset(Vector3 n1, Vector3 n2, Vector3 t, double r) {
    final m00 = n1.x, m01 = n1.y, m02 = n1.z;
    final m10 = n2.x, m11 = n2.y, m12 = n2.z;
    final m20 = t.x, m21 = t.y, m22 = t.z;
    final det =
        m00 * (m11 * m22 - m12 * m21) -
        m01 * (m10 * m22 - m12 * m20) +
        m02 * (m10 * m21 - m11 * m20);
    if (det.abs() < 1e-12) return null;
    final inv = 1.0 / det;
    final rhs = Vector3(-r, -r, 0);
    return Vector3(
      inv *
          (rhs.x * (m11 * m22 - m12 * m21) -
              m01 * (rhs.y * m22 - m12 * rhs.z) +
              m02 * (rhs.y * m21 - m11 * rhs.z)),
      inv *
          (m00 * (rhs.y * m22 - m12 * rhs.z) -
              rhs.x * (m10 * m22 - m12 * m20) +
              m02 * (m10 * rhs.z - rhs.y * m20)),
      inv *
          (m00 * (m11 * rhs.z - rhs.y * m21) -
              m01 * (m10 * rhs.z - rhs.y * m20) +
              rhs.x * (m10 * m21 - m11 * m20)),
    );
  }

  static BrepSolid _filletPlanePlane(
    BrepSolid s,
    int edgeId,
    int i1,
    int i2,
    double r, {
    bool inverse = false,
  }) {
    final e = s.edges[edgeId];
    if (e.curve is! BrepLine) {
      throw BrepUnsupported('Plane–plane fillet needs a line edge');
    }
    final f1 = s.faces[i1];
    final f2 = s.faces[i2];
    final p1 = f1.surface as BrepPlane;
    final p2 = f2.surface as BrepPlane;
    var t = (s.vertices[e.v1] - s.vertices[e.v0]);
    if (t.length < 1e-9) {
      throw BrepUnsupported('Degenerate fillet edge');
    }
    t.normalize();
    final n1 = f1.reversed ? -p1.normal : p1.normal;
    final n2 = f2.reversed ? -p2.normal : p2.normal;
    // The signed turn across an outward-oriented edge distinguishes a convex
    // ridge from a reentrant corner, including edges around hole boundaries.
    final loop = f1.loops.firstWhere(
      (l) => l.coedges.any((c) => c.edge == edgeId),
    );
    final use = loop.coedges.firstWhere((c) => c.edge == edgeId);
    var alongFace = use.reversed ? -t : t;
    // Boolean cuts can retain a planar surface with reversed=true and trims
    // wound around its underlying normal. Measure the trim's winding against
    // the outward normal before deciding whether the edge is concave.
    final points = BrepTessellate.sampleLoop(s, loop, 32);
    var area = Vector3.zero();
    for (var i = 0; i < points.length; i++) {
      area += (points[i] - p1.origin).cross(
        points[(i + 1) % points.length] - p1.origin,
      );
    }
    final outer = identical(loop, f1.loops.first);
    if ((area.dot(n1) < 0) == outer) alongFace = -alongFace;
    final concave = n1.cross(n2).dot(alongFace) < 0;
    if (inverse) {
      // Inverse is a circular groove centered on the selected edge. Moving a
      // tangent fillet's center outside the solid creates invalid trim wires.
      final frame = chamferFrame(s, edgeId);
      if (frame == null)
        throw BrepUnsupported('Cannot determine edge interior');
      final a = s.vertices[e.v0], b = s.vertices[e.v1];
      final a1 = a + frame.$2 * r, b1 = b + frame.$2 * r;
      final a2 = a + frame.$3 * r, b2 = b + frame.$3 * r;
      if ([
            a1,
            b1,
          ].any((p) => BrepClassify.containsPlaneFace(s, f1, p) == false) ||
          [
            a2,
            b2,
          ].any((p) => BrepClassify.containsPlaneFace(s, f2, p) == false)) {
        throw BrepUnsupported('Inverse fillet reaches another face boundary');
      }
      return _replaceEdgeWithBlend(
        s,
        edgeId,
        i1,
        i2,
        a1,
        b1,
        a2,
        b2,
        BrepCylinder(origin: a, axis: t, radius: r),
        p1,
        p2,
        blendReversed: true,
        majorArc: concave,
      );
    }
    final reverseCylinder = concave != inverse;
    final off = _centerOffset(n1, n2, t, reverseCylinder ? -r : r);
    if (off == null) {
      throw BrepUnsupported('Could not place fillet cylinder');
    }
    final a = s.vertices[e.v0];
    final b = s.vertices[e.v1];
    final c0 = a + off;
    final cyl = BrepCylinder(origin: c0, axis: t, radius: r);
    // Offset edges on each plane: project cylinder-plane tangent generators.
    final a1 = a + _tangentOffset(n1, t, a, c0);
    final b1 = b + (a1 - a);
    final a2 = a + _tangentOffset(n2, t, a, c0);
    final b2 = b + (a2 - a);

    return _replaceEdgeWithBlend(
      s,
      edgeId,
      i1,
      i2,
      a1,
      b1,
      a2,
      b2,
      cyl,
      p1,
      p2,
      blendReversed: reverseCylinder,
    );
  }

  static Vector3 _tangentOffset(Vector3 n, Vector3 t, Vector3 a, Vector3 c0) {
    // Project the cylinder axis onto the incident plane. Its distance from
    // the original edge is r / tan(interiorAngle / 2), equal to r only for
    // a right angle. Preserve that length so the trim is on both surfaces.
    var dir = (c0 - a);
    dir = dir - t * dir.dot(t);
    dir = dir - n * dir.dot(n);
    return dir;
  }

  static BrepSolid _replaceEdgeWithBlend(
    BrepSolid s,
    int edgeId,
    int i1,
    int i2,
    Vector3 a1,
    Vector3 b1,
    Vector3 a2,
    Vector3 b2,
    BrepCylinder cyl,
    BrepPlane p1,
    BrepPlane p2, {
    BrepSurface? blendSurface,
    bool blendReversed = false,
    bool majorArc = false,
    bool round = true,
  }) {
    final e = s.edges[edgeId];
    final v0 = e.v0;
    final v1 = e.v1;
    final endA = _blendEndPlane(s, edgeId, i1, i2, v0, cyl.axis);
    final endB = _blendEndPlane(s, edgeId, i1, i2, v1, cyl.axis);
    a1 = _clipBlendContact(a1, cyl.axis, endA);
    a2 = _clipBlendContact(a2, cyl.axis, endA);
    b1 = _clipBlendContact(b1, cyl.axis, endB);
    b2 = _clipBlendContact(b2, cyl.axis, endB);
    final bld = BrepBuilder();
    int copyV(Vector3 p) {
      for (var i = 0; i < bld.vertices.length; i++) {
        if ((bld.vertices[i] - p).length < 1e-7) return i;
      }
      return bld.vertex(p);
    }

    final edgeOf = <String, int>{};
    (int, bool) sew(int i, int j, BrepCurve curve) {
      final mid = curve.evaluate(0.5 * (curve.t0 + curve.t1));
      final q =
          '${mid.x.toStringAsFixed(5)}_${mid.y.toStringAsFixed(5)}_${mid.z.toStringAsFixed(5)}';
      final lo = i <= j ? i : j;
      final hi = i <= j ? j : i;
      final k = '$lo-$hi-$q';
      final existing = edgeOf[k];
      if (existing != null) {
        final e = bld.edges[existing];
        return (existing, e.v0 != i);
      }
      final id = bld.edge(i, j, curve);
      edgeOf[k] = id;
      return (id, false);
    }

    final va1 = copyV(a1);
    final vb1 = copyV(b1);
    final va2 = copyV(a2);
    final vb2 = copyV(b2);

    int remap(int fi, int oldV) {
      if (fi == i1) {
        if (oldV == v0) return va1;
        if (oldV == v1) return vb1;
      } else if (fi == i2) {
        if (oldV == v0) return va2;
        if (oldV == v1) return vb2;
      }
      return copyV(s.vertices[oldV]);
    }

    (int, bool) edgeMapped(int fi, int old) {
      final oe = s.edges[old];
      final a = remap(fi, oe.v0);
      final b = remap(fi, oe.v1);
      final moved =
          (bld.vertices[a] - s.vertices[oe.v0]).length > 1e-7 ||
          (bld.vertices[b] - s.vertices[oe.v1]).length > 1e-7;
      final curve = (!moved || oe.curve is BrepCircle)
          ? oe.curve
          : BrepLine(bld.vertices[a], bld.vertices[b]);
      return sew(a, b, curve);
    }

    final eOff1 = sew(va1, vb1, BrepLine(a1, b1)).$1;
    final eOff2 = sew(va2, vb2, BrepLine(a2, b2)).$1;
    var tA = _cylU(cyl, a1);
    var tB = _cylU(cyl, a2);
    if (tB < tA) tB += math.pi * 2;
    if (tB - tA > math.pi) tB -= math.pi * 2;
    if (majorArc) tB += (tB > tA ? -1 : 1) * math.pi * 2;
    final arcA = round
        ? sew(va1, va2, _blendEndCurve(cyl, endA, a1, a2, tA, tB)).$1
        : sew(va1, va2, BrepLine(a1, a2)).$1;
    var tA1 = _cylU(cyl, b1);
    var tB1 = _cylU(cyl, b2);
    if (tB1 < tA1) tB1 += math.pi * 2;
    if (tB1 - tA1 > math.pi) tB1 -= math.pi * 2;
    if (majorArc) tB1 += (tB1 > tA1 ? -1 : 1) * math.pi * 2;
    final arcB = round
        ? sew(vb1, vb2, _blendEndCurve(cyl, endB, b1, b2, tA1, tB1)).$1
        : sew(vb1, vb2, BrepLine(b1, b2)).$1;

    bool usesVertex(BrepFace face, int v) {
      for (final loop in face.loops) {
        for (final c in loop.coedges) {
          final oe = s.edges[c.edge];
          if (oe.v0 == v || oe.v1 == v) return true;
        }
      }
      return false;
    }

    int startOf(BrepCoedge c) {
      final oe = s.edges[c.edge];
      return c.reversed ? oe.v1 : oe.v0;
    }

    int endOf(BrepCoedge c) {
      final oe = s.edges[c.edge];
      return c.reversed ? oe.v0 : oe.v1;
    }

    int vertexOffset(int vertex, BrepCoedge c) {
      if (vertex == v0) return _endOffset(c, i1, i2, va1, va2, s);
      if (vertex == v1) return _endOffset(c, i1, i2, vb1, vb2, s);
      return copyV(s.vertices[vertex]);
    }

    for (var fi = 0; fi < s.faces.length; fi++) {
      final f = s.faces[fi];
      if (fi == i1 || fi == i2) {
        final off = fi == i1 ? eOff1 : eOff2;
        final loops = <BrepLoop>[];
        for (final loop in f.loops) {
          final co = <BrepCoedge>[];
          for (final c in loop.coedges) {
            if (c.edge == edgeId) {
              co.add(BrepCoedge(off, reversed: c.reversed));
            } else {
              final e = edgeMapped(fi, c.edge);
              co.add(BrepCoedge(e.$1, reversed: e.$2 != c.reversed));
            }
          }
          loops.add(BrepLoop(co));
        }
        bld.face(f.surface, loops, reversed: f.reversed);
        continue;
      }
      if (!usesVertex(f, v0) && !usesVertex(f, v1)) {
        final loops = <BrepLoop>[];
        for (final loop in f.loops) {
          final co = <BrepCoedge>[];
          for (final c in loop.coedges) {
            final e = edgeMapped(fi, c.edge);
            co.add(BrepCoedge(e.$1, reversed: e.$2 != c.reversed));
          }
          loops.add(BrepLoop(co));
        }
        bld.face(f.surface, loops, reversed: f.reversed);
        continue;
      }
      final loops = <BrepLoop>[];
      for (final loop in f.loops) {
        final co = <BrepCoedge>[];
        for (var i = 0; i < loop.coedges.length; i++) {
          final c = loop.coedges[i];
          final ns = startOf(c);
          final ne = endOf(c);
          final rs = vertexOffset(ns, c);
          final re = vertexOffset(ne, c);
          final orig = s.edges[c.edge].curve;
          final moved =
              (bld.vertices[rs] - s.vertices[ns]).length > 1e-7 ||
              (bld.vertices[re] - s.vertices[ne]).length > 1e-7;
          if (!moved || orig is BrepCircle) {
            // Reused curves retain their original parameter direction. The
            // surrounding face may traverse them backwards; keep that on the
            // coedge instead of reversing vertex IDs under an unchanged curve.
            final e = c.reversed ? sew(re, rs, orig) : sew(rs, re, orig);
            co.add(BrepCoedge(e.$1, reversed: e.$2 != c.reversed));
          } else {
            final e = sew(rs, re, BrepLine(bld.vertices[rs], bld.vertices[re]));
            co.add(BrepCoedge(e.$1, reversed: e.$2));
          }
          if (ne != v0 && ne != v1) continue;
          final next = loop.coedges[(i + 1) % loop.coedges.length];
          final to = vertexOffset(ne, next);
          if (re == to) continue;
          if (ne == v0) {
            co.add(BrepCoedge(arcA, reversed: re != va1));
          } else {
            co.add(BrepCoedge(arcB, reversed: re != vb1));
          }
        }
        if (co.isNotEmpty) loops.add(BrepLoop(co));
      }
      bld.face(f.surface, loops, reversed: f.reversed);
    }

    final blend = blendSurface ?? cyl;
    final blendLoop = [
      BrepCoedge(eOff1),
      BrepCoedge(arcB),
      BrepCoedge(eOff2, reversed: true),
      BrepCoedge(arcA, reversed: true),
    ];
    final incidentUse = s.faces[i1].loops
        .expand((l) => l.coedges)
        .firstWhere((c) => c.edge == edgeId);
    // Adjacent faces must traverse their shared boundary in opposite directions.
    bld.face(blend, [
      BrepLoop(
        incidentUse.reversed
            ? blendLoop
            : [
                for (final c in blendLoop.reversed)
                  BrepCoedge(c.edge, reversed: !c.reversed),
              ],
      ),
    ], reversed: blendReversed);
    return bld.build();
  }

  static int _endOffset(
    BrepCoedge c,
    int i1,
    int i2,
    int o1,
    int o2,
    BrepSolid s,
  ) {
    final faces = facesForEdge(s, c.edge);
    if (faces.contains(i1)) return o1;
    if (faces.contains(i2)) return o2;
    return o1;
  }

  static double rFor(BrepCylinder c) => c.radius;

  static double _cylU(BrepCylinder c, Vector3 p) {
    final uv = c.parameters(p);
    var u = uv?.$1 ?? 0.0;
    if (u < 0) u += math.pi * 2;
    return u;
  }

  static BrepSolid _filletPlaneCylinder(
    BrepSolid s,
    int edgeId,
    BrepFace planeFace,
    BrepFace cylFace,
    double r, {
    bool inverse = false,
  }) {
    final plane = planeFace.surface as BrepPlane;
    final cyl = cylFace.surface as BrepCylinder;
    final circle = s.edges[edgeId].curve as BrepCircle;
    final nCap = planeFace.reversed ? -plane.normal : plane.normal;
    if (inverse) {
      final radialSign = cylFace.reversed ? -1.0 : 1.0;
      if (r >= cyl.radius - 1e-7 && !cylFace.reversed)
        throw BrepUnsupported('Inverse fillet reaches the cylinder axis');
      return _rebuildCircularRim(
        s,
        edgeId,
        s.faces.indexOf(planeFace),
        s.faces.indexOf(cylFace),
        BrepCircle(
          center: circle.center,
          normal: circle.normal,
          radius: cyl.radius - radialSign * r,
          xDir: circle.xDir,
        ),
        BrepCircle(
          center: circle.center - nCap * r,
          normal: circle.normal,
          radius: cyl.radius,
          xDir: circle.xDir,
        ),
        BrepTorus(
          center: circle.center,
          axis: cyl.axis,
          majorRadius: cyl.radius,
          minorRadius: r,
          xDir: cyl.xDir,
        ),
      );
    }
    final frame = chamferFrame(s, edgeId)!;
    final radial = (frame.$1 - circle.center).normalized();
    final into = frame.$3;
    final major = cyl.radius + frame.$2.dot(radial).sign * r;
    if (major < -1e-7) {
      throw BrepUnsupported('Fillet radius too large for cylinder');
    }
    if (major.abs() <= 1e-7) {
      return _rebuildCircularRim(
        s,
        edgeId,
        s.faces.indexOf(planeFace),
        s.faces.indexOf(cylFace),
        BrepCircle(
          center: circle.center,
          normal: circle.normal,
          radius: 0,
          xDir: circle.xDir,
        ),
        BrepCircle(
          center: circle.center + into * r,
          normal: circle.normal,
          radius: cyl.radius,
          xDir: circle.xDir,
        ),
        BrepSphere(center: circle.center + into * r, radius: r),
      );
    }
    final torus = BrepTorus(
      center: circle.center + into * r,
      axis: cyl.axis,
      majorRadius: major,
      minorRadius: r,
      xDir: cyl.xDir,
    );
    final inner = BrepCircle(
      center: circle.center,
      normal: circle.normal,
      radius: major,
      xDir: circle.xDir,
    );
    return _rebuildCylinderFilletCap(s, cyl, plane, circle, inner, torus, r);
  }

  static BrepSolid _rebuildCylinderFilletCap(
    BrepSolid s,
    BrepCylinder cyl,
    BrepPlane cap,
    BrepCircle outer,
    BrepCircle inner,
    BrepTorus torus,
    double r,
  ) {
    final edgeId = s.edges.indexWhere((e) => identical(e.curve, outer));
    final planeId = s.faces.indexWhere((f) => identical(f.surface, cap));
    final cylinderId = s.faces.indexWhere((f) => identical(f.surface, cyl));
    if (edgeId < 0 || planeId < 0 || cylinderId < 0)
      throw BrepUnsupported('Cannot identify the shared rim');
    final wall = BrepCircle(
      center: torus.center,
      normal: outer.normal,
      radius: cyl.radius,
      xDir: outer.xDir,
    );
    return _rebuildCircularRim(
      s,
      edgeId,
      planeId,
      cylinderId,
      inner,
      wall,
      torus,
    );
  }

  static BrepSolid _rebuildCircularRim(
    BrepSolid s,
    int edgeId,
    int planeId,
    int cylinderId,
    BrepCircle inner,
    BrepCircle wall,
    BrepSurface blend,
  ) {
    final circle = s.edges[edgeId].curve as BrepCircle;
    final cap = s.faces[planeId], wallFace = s.faces[cylinderId];
    final collapsed = inner.radius.abs() < 1e-8;
    if (collapsed &&
        (cap.loops.length != 1 || cap.loops.single.coedges.length != 1)) {
      throw BrepUnsupported('Collapsing circular cap meets another boundary');
    }
    final n = cap.outwardNormal(0, 0);
    final displacement = wall.center - circle.center;
    final depth = displacement.length;
    final into = depth > 1e-9 ? displacement / depth : Vector3.zero();
    if (depth <= 1e-9)
      throw BrepUnsupported('Blend must extend into the selected wall');
    for (var i = 0; i < 32; i++) {
      final t = inner.t0 + (inner.t1 - inner.t0) * i / 32;
      if (BrepClassify.containsPlaneFace(s, cap, inner.evaluate(t)) == false) {
        throw BrepUnsupported('Blend reaches another cap boundary');
      }
    }
    for (final loop in wallFace.loops) {
      for (final co in loop.coedges) {
        if (co.edge == edgeId) continue;
        final curve = s.edges[co.edge].curve;
        if (curve is BrepCircle) {
          final available = (curve.center - circle.center).dot(into);
          if (available > 1e-7 && depth >= available - 1e-7)
            throw BrepUnsupported('Blend reaches the opposite wall boundary');
        }
      }
    }
    final normalized = BrepHeal.sewClosedRims(s, removeSeams: true);
    if (normalized.faces.length != s.faces.length)
      throw BrepUnsupported('Blend input has collapsed faces');
    final edge = s.edges[edgeId];
    if (edge.v0 != edge.v1 || !circle.isClosed) {
      throw BrepUnsupported(
        'A partial circular rim needs endpoint corner blends',
      );
    }
    // Keep every unrelated face and hole. Replacing the whole body with a
    // cylinder here discarded the surrounding block when rounding a bore.
    final b = BrepBuilder();
    for (final v in s.vertices) {
      b.vertex(v);
    }
    b.vertices[edge.v0] = inner.evaluate(inner.t0);
    for (var i = 0; i < s.edges.length; i++) {
      final e = s.edges[i];
      b.edge(e.v0, e.v1, i == edgeId ? inner : e.curve, origin: e.origin);
    }
    final wv = b.vertex(wall.evaluate(wall.t0));
    final we = b.edge(wv, wv, wall);
    BrepCoedge? capUse, wallUse;
    for (var fi = 0; fi < normalized.faces.length; fi++) {
      final f = normalized.faces[fi];
      final loops = <BrepLoop>[];
      for (final loop in f.loops) {
        final uses = <BrepCoedge>[];
        for (final co in loop.coedges) {
          if (co.edge == edgeId && fi == planeId) capUse = co;
          if (co.edge == edgeId && fi == cylinderId) {
            wallUse = BrepCoedge(
              we,
              reversed: collapsed ? wall.normal.dot(n) > 0 : co.reversed,
            );
            uses.add(wallUse);
          } else {
            uses.add(co);
          }
        }
        loops.add(BrepLoop(uses));
      }
      if (!collapsed || fi != planeId)
        b.face(f.surface, loops, reversed: f.reversed, origin: f.origin);
    }
    if (capUse == null || wallUse == null)
      throw BrepUnsupported('Fillet rim has no shared cap and wall');
    final rimPoint = circle.evaluate(circle.t0),
        wallUv = wallFace.surface.parameters(rimPoint)!;
    final expected = n + wallFace.outwardNormal(wallUv.$1, wallUv.$2);
    final middle = (inner.evaluate(inner.t0) + wall.evaluate(wall.t0)) * .5;
    final blendUv = blend.parameters(middle)!;
    b.face(blend, [
      if (!collapsed)
        BrepLoop([BrepCoedge(edgeId, reversed: !capUse.reversed)]),
      BrepLoop([BrepCoedge(we, reversed: !wallUse.reversed)]),
    ], reversed: blend.normalAt(blendUv.$1, blendUv.$2).dot(expected) < 0);
    return b.build();
  }

  static BrepSolid _chamferPlanePlane(
    BrepSolid s,
    int edgeId,
    int i1,
    int i2,
    double d,
    ChamferOptions options,
  ) {
    if (BrepConvexPlanes.supports(s)) {
      return _chamferConvexEdges(s, [edgeId], d, options);
    }
    final e = s.edges[edgeId];
    final f1 = s.faces[i1];
    final f2 = s.faces[i2];
    final p1 = f1.surface as BrepPlane;
    final p2 = f2.surface as BrepPlane;
    var t = s.vertices[e.v1] - s.vertices[e.v0];
    t.normalize();
    final n1 = f1.reversed ? -p1.normal : p1.normal;
    final n2 = f2.reversed ? -p2.normal : p2.normal;
    final a = s.vertices[e.v0];
    final b = s.vertices[e.v1];
    final frame = chamferFrame(s, edgeId)!;
    final dir1 = frame.$2, dir2 = frame.$3;
    final endA = _blendEndPlane(s, edgeId, i1, i2, e.v0, t);
    final endB = _blendEndPlane(s, edgeId, i1, i2, e.v1, t);
    final a1 = _clipBlendContact(a + dir1 * d, t, endA);
    final b1 = _clipBlendContact(b + dir1 * d, t, endB);
    final beta = math.acos(dir1.dot(dir2).clamp(-1.0, 1.0));
    final d2 = options.secondDistance(d, beta);
    final a2 = _clipBlendContact(a + dir2 * d2, t, endA);
    final b2 = _clipBlendContact(b + dir2 * d2, t, endB);
    var cutNormal = (b1 - a1).cross(a2 - a1);
    if (cutNormal.dot(n1 + n2) < 0) cutNormal = -cutNormal;
    if (s.faces.any((f) => f.surface is! BrepPlane) && _isConvexBlendSolid(s)) {
      final pieces = BrepSlice.parts(s, a1, cutNormal);
      final kept = pieces
          .where(
            (part) =>
                part.vertices.isNotEmpty &&
                part.vertices.fold(
                      0.0,
                      (sum, p) => sum + (p - a1).dot(cutNormal),
                    ) <
                    -1e-7,
          )
          .toList();
      if (kept.length == 1) return kept.single;
      throw BrepUnsupported('Chamfer would consume the remaining body');
    }
    if ([
          a1,
          b1,
        ].any((p) => BrepClassify.containsPlaneFace(s, f1, p) == false) ||
        [
          a2,
          b2,
        ].any((p) => BrepClassify.containsPlaneFace(s, f2, p) == false)) {
      throw BrepUnsupported('Chamfer setback reaches another face boundary');
    }
    var nCh = (b1 - a1).cross(a2 - a1);
    if (nCh.dot(n1 + n2) < 0) nCh = -nCh;
    if (nCh.length < 1e-9) {
      throw BrepUnsupported('Chamfer planes are opposite');
    }
    final plane = BrepPlane(origin: a1, normal: nCh.normalized());
    final dummy = BrepCylinder(origin: a, axis: t, radius: d);
    return _replaceEdgeWithBlend(
      s,
      edgeId,
      i1,
      i2,
      a1,
      b1,
      a2,
      b2,
      dummy,
      p1,
      p2,
      blendSurface: plane,
      round: false,
    );
  }

  static BrepSolid _chamferConvexEdges(
    BrepSolid s,
    List<int> edges,
    double distance,
    ChamferOptions options,
  ) {
    final cuts = <BrepPlane>[];
    for (final id in edges.toSet()) {
      final faces = facesForEdge(s, id), frame = chamferFrame(s, id);
      if (faces.length != 2 || frame == null) {
        throw BrepUnsupported('Chamfer needs two adjacent planar faces');
      }
      final edge = s.edges[id];
      final a = s.vertices[edge.v0], tangent = s.vertices[edge.v1] - a;
      final dir1 = frame.$2, dir2 = frame.$3;
      final d2 = options.secondDistance(
        distance,
        math.acos(dir1.dot(dir2).clamp(-1.0, 1.0)),
      );
      // Endpoints may legitimately change adjacency at a bevel corner. The
      // middle contacts must nevertheless lie on both original support faces.
      if (BrepClassify.containsPlaneFace(
                s,
                s.faces[faces[0]],
                frame.$1 + dir1 * distance,
              ) ==
              false ||
          BrepClassify.containsPlaneFace(
                s,
                s.faces[faces[1]],
                frame.$1 + dir2 * d2,
              ) ==
              false) {
        throw BrepUnsupported('Chamfer distance exceeds the adjacent face');
      }
      final p = a + dir1 * distance, q = a + dir2 * d2;
      var normal = tangent.cross(q - p);
      final outward =
          s.faces[faces[0]].outwardNormal(0, 0) +
          s.faces[faces[1]].outwardNormal(0, 0);
      if (normal.dot(outward) < 0) normal = -normal;
      if (normal.length < 1e-12) {
        throw BrepUnsupported('Chamfer planes are opposite');
      }
      cuts.add(BrepPlane(origin: p, normal: normal));
    }
    final result = BrepConvexPlanes.clip(s, cuts);
    return result.solid;
  }

  static BrepSolid _blendPlaneCone(
    BrepSolid s,
    int edgeId,
    BrepFace cap,
    BrepFace face,
    double d, {
    required bool round,
    ChamferOptions options = const ChamferOptions(),
  }) {
    final cone = face.surface as BrepCone,
        circle = s.edges[edgeId].curve as BrepCircle;
    final n = cap.outwardNormal(0, 0), axis = cone.axis, radial = circle.xDir;
    if (n.dot(axis).abs() < 1 - 1e-8)
      throw BrepUnsupported('Cone rim must be perpendicular to its axis');
    final p = circle.center + radial * circle.radius, uv = cone.parameters(p)!;
    final nw = face.outwardNormal(uv.$1, uv.$2);
    final sign = nw.dot(radial) > 0 ? 1.0 : -1.0;
    final first = -radial * sign;
    var along = (axis + radial * ((cone.topRadius - cone.radius) / cone.height))
        .normalized();
    if (along.dot(n) > 0) along = -along;
    final beta = math.acos(first.dot(along).clamp(-1.0, 1.0));
    late Vector3 a, b, center;
    late BrepSurface blend;
    if (round) {
      final depth = -n * d;
      final slope = (cone.topRadius - cone.radius) / cone.height;
      final radius =
          circle.radius +
          slope * depth.dot(axis) -
          sign * d * math.sqrt(1 + slope * slope);
      center = circle.center + depth + radial * radius;
      a = center + n * d;
      b = center + nw * d;
      if (radius <= d * 1e-6)
        throw BrepUnsupported('Fillet radius reaches the cone axis');
      blend = BrepTorus(
        center: circle.center + depth,
        axis: axis,
        majorRadius: radius,
        minorRadius: d,
        xDir: circle.xDir,
      );
    } else {
      final second = options.secondDistance(d, beta);
      a = p + first * d;
      b = p + along * second;
      final origin = circle.center + axis * (b - circle.center).dot(axis);
      final height = (a - b).dot(axis);
      blend = BrepCone(
        origin: origin,
        axis: axis,
        radius: (b - origin).length,
        topRadius: (a - circle.center).length,
        height: height,
        xDir: circle.xDir,
      );
    }
    final ar = (a - circle.center).dot(radial),
        br = (b - circle.center).dot(radial);
    if (ar <= 1e-7 || br <= 1e-7)
      throw BrepUnsupported('Blend reaches the cone axis');
    final z = (b - cone.origin).dot(axis);
    if (z < 1e-7 || z > cone.height - 1e-7)
      throw BrepUnsupported('Blend reaches the opposite end of the cone');
    final inner = BrepCircle(
      center: circle.center,
      normal: circle.normal,
      radius: ar,
      xDir: circle.xDir,
    );
    final wall = BrepCircle(
      center: circle.center + axis * (b - circle.center).dot(axis),
      normal: circle.normal,
      radius: br,
      xDir: circle.xDir,
    );
    return _rebuildCircularRim(
      s,
      edgeId,
      s.faces.indexOf(cap),
      s.faces.indexOf(face),
      inner,
      wall,
      blend,
    );
  }

  static BrepSolid _chamferPlaneCylinder(
    BrepSolid s,
    int edgeId,
    BrepFace planeFace,
    BrepFace cylFace,
    double d,
    ChamferOptions options,
  ) {
    final cyl = cylFace.surface as BrepCylinder;
    final circle = s.edges[edgeId].curve as BrepCircle;
    final depth = options.secondDistance(d, math.pi / 2);
    final frame = chamferFrame(s, edgeId)!;
    final radial = (frame.$1 - circle.center).normalized();
    final major = cyl.radius + frame.$2.dot(radial).sign * d;
    if (major < -1e-7) {
      throw BrepUnsupported('Chamfer too large for cylinder');
    }
    final into = frame.$3;
    final cone = BrepCone(
      origin: circle.center + into * depth,
      axis: -into,
      radius: cyl.radius,
      topRadius: major,
      height: depth,
      xDir: cyl.xDir,
    );
    final inner = BrepCircle(
      center: circle.center,
      normal: circle.normal,
      radius: major,
      xDir: circle.xDir,
    );
    return _rebuildCircularRim(
      s,
      edgeId,
      s.faces.indexOf(planeFace),
      s.faces.indexOf(cylFace),
      inner,
      BrepCircle(
        center: circle.center + into * depth,
        normal: circle.normal,
        radius: cyl.radius,
        xDir: circle.xDir,
      ),
      cone,
    );
  }
}
