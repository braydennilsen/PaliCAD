import 'dart:math' as math;
import 'dart:collection';
import 'dart:isolate';

import 'package:vector_math/vector_math_64.dart';

import 'mesh_data.dart';

/// Spatial broad phase for T-junctions. Segment/box queries retain the exact
/// point-on-segment predicate while avoiding scans of entire coordinate slabs.
class _SegmentPointIndex {
  _SegmentPointIndex(List<(int, Vector3)> points) {
    if (points.isEmpty) {
      leaves = points;
      return;
    }
    low = points.first.$2.clone();
    high = low.clone();
    for (final p in points) {
      for (var axis = 0; axis < 3; axis++) {
        low[axis] = math.min(low[axis], p.$2[axis]);
        high[axis] = math.max(high[axis], p.$2[axis]);
      }
    }
    if (points.length <= 16) {
      leaves = points;
      return;
    }
    final extent = high - low;
    var axis = 0;
    if (extent.y > extent.x) axis = 1;
    if (extent.z > extent[axis]) axis = 2;
    points.sort((a, b) => a.$2[axis].compareTo(b.$2[axis]));
    final mid = points.length ~/ 2;
    left = _SegmentPointIndex(points.sublist(0, mid));
    right = _SegmentPointIndex(points.sublist(mid));
  }
  Vector3 low = Vector3.zero(), high = Vector3.zero();
  List<(int, Vector3)>? leaves;
  _SegmentPointIndex? left, right;

  List<(int, Vector3)> nearSegment(Vector3 a, Vector3 b) {
    final result = <(int, Vector3)>[];
    _collect(a, b - a, result);
    return result;
  }

  void _collect(Vector3 origin, Vector3 direction, List<(int, Vector3)> out) {
    const eps = MeshFaceRegions.planarGrowOffsetEpsilon;
    var start = 0.0, end = 1.0;
    for (var axis = 0; axis < 3; axis++) {
      final d = direction[axis], o = origin[axis];
      if (d == 0) {
        if (o < low[axis] - eps || o > high[axis] + eps) return;
      } else {
        final a = (low[axis] - eps - o) / d, b = (high[axis] + eps - o) / d;
        start = math.max(start, math.min(a, b));
        end = math.min(end, math.max(a, b));
        if (start > end) return;
      }
    }
    final points = leaves;
    if (points != null) {
      out.addAll(points);
    } else {
      left!._collect(origin, direction, out);
      right!._collect(origin, direction, out);
    }
  }
}

/// A face on a plane must have its centroid in that plane's tolerance slab.
/// This broad phase avoids scanning every triangle for every supporting plane.
/// The original per-vertex plane predicate still makes the final decision.
typedef _PlaneCandidate = ({
  int id,
  Vector3 center,
  Vector3 normal,
  double radius,
});

class _FacePlaneIndex {
  _FacePlaneIndex(List<_PlaneCandidate> faces) {
    low = faces.first.center.clone();
    high = low.clone();
    normalLow = Vector3.all(double.infinity);
    normalHigh = Vector3.all(double.negativeInfinity);
    final directionLow = faces.first.normal.clone(),
        directionHigh = directionLow.clone();
    for (final f in faces) {
      for (var axis = 0; axis < 3; axis++) {
        low[axis] = math.min(low[axis], f.center[axis]);
        high[axis] = math.max(high[axis], f.center[axis]);
        normalLow[axis] = math.min(normalLow[axis], f.normal[axis] - f.radius);
        normalHigh[axis] = math.max(
          normalHigh[axis],
          f.normal[axis] + f.radius,
        );
        directionLow[axis] = math.min(directionLow[axis], f.normal[axis]);
        directionHigh[axis] = math.max(directionHigh[axis], f.normal[axis]);
      }
    }
    if (faces.length <= 16) {
      leaves = faces;
    } else {
      final directionSpan = directionHigh - directionLow;
      final radius = faces.fold(0.0, (double r, f) => math.max(r, f.radius));
      final byDirection =
          directionSpan.length2 > math.max(1e-16, 16 * radius * radius);
      final span = byDirection ? directionSpan : high - low;
      var axis = 0;
      if (span.y > span.x) axis = 1;
      if (span.z > span[axis]) axis = 2;
      faces.sort(
        (a, b) => (byDirection ? a.normal[axis] : a.center[axis]).compareTo(
          byDirection ? b.normal[axis] : b.center[axis],
        ),
      );
      final mid = faces.length ~/ 2;
      left = _FacePlaneIndex(faces.sublist(0, mid));
      right = _FacePlaneIndex(faces.sublist(mid));
    }
  }
  late final Vector3 low, high, normalLow, normalHigh;
  List<_PlaneCandidate>? leaves;
  _FacePlaneIndex? left, right;

  bool _exhausted = false;
  void visitNearPlane(
    Vector3 origin,
    Vector3 normal,
    double eps,
    List<int> assigned,
    void Function(int) visit,
  ) {
    if (_exhausted) return;
    bool directionPossible(double sign) {
      for (var axis = 0; axis < 3; axis++) {
        final value = normal[axis] * sign;
        if (value < normalLow[axis] || value > normalHigh[axis]) return false;
      }
      return true;
    }

    if (!directionPossible(1) && !directionPossible(-1)) return;
    var minD = 0.0, maxD = 0.0;
    for (var axis = 0; axis < 3; axis++) {
      final a = (low[axis] - origin[axis]) * normal[axis];
      final b = (high[axis] - origin[axis]) * normal[axis];
      minD += math.min(a, b);
      maxD += math.max(a, b);
    }
    if (minD > eps || maxD < -eps) return;
    final entries = leaves;
    if (entries != null) {
      for (final f in entries) {
        if (assigned[f.id] < 0 && (f.center - origin).dot(normal).abs() <= eps)
          visit(f.id);
      }
      _exhausted = entries.every((f) => assigned[f.id] >= 0);
    } else {
      left!.visitNearPlane(origin, normal, eps, assigned, visit);
      right!.visitNearPlane(origin, normal, eps, assigned, visit);
      _exhausted = left!._exhausted && right!._exhausted;
    }
  }
}

/// Prepared, isolate-sendable mesh data. Expando caches themselves are isolate
/// local: explicitly carry the cache with its mesh, then install it on receipt.
class PreparedMeshRegions {
  PreparedMeshRegions._(this._cache);
  final _MeshRegionCache _cache;
  MeshData get mesh => _cache.mesh;

  MeshData adopt() {
    MeshFaceRegions._caches[mesh] = _cache;
    return mesh;
  }
}

/// Coplanar face-region queries over render meshes, cached per [MeshData].
///
/// Meshes triangulate flat faces into multiple triangles. These helpers group edge-connected, coplanar
/// triangles back into logical faces for selection highlighting and
/// push/pull extrusion. Edges are keyed by quantized vertex positions, not
/// indices, because mesh triangulation duplicates vertices between faces.
class MeshFaceRegions {
  MeshFaceRegions._();

  static final Expando<_MeshRegionCache> _caches = Expando('meshFaceRegions');

  static _MeshRegionCache _cacheFor(MeshData mesh) =>
      _caches[mesh] ??= _MeshRegionCache(mesh);

  /// Run in an import worker before transferring ownership with Isolate.exit.
  /// [onProgress] counts the four display-cache preparation stages.
  static PreparedMeshRegions prepare(
    MeshData mesh, {
    void Function(int completed, int total)? onProgress,
  }) {
    final cache = _cacheFor(mesh);
    onProgress?.call(0, 4);
    cache.planarFaceNormals();
    onProgress?.call(1, 4);
    cache.sharpFeatureEdges();
    onProgress?.call(2, 4);
    cache.cornerNormals();
    onProgress?.call(3, 4);
    cache.silhouetteCandidates();
    cache.displayPrepared = true;
    onProgress?.call(4, 4);
    return PreparedMeshRegions._(cache);
  }

  static final Expando<Future<void>> _preparing = Expando('meshPreparation');

  /// A restored document or a modelling operation may supply a cold mesh to
  /// the viewport. Prepare it off-thread too, deduplicating concurrent draws.
  static Future<void> prepareAsync(MeshData mesh) {
    if (_caches[mesh]?.displayPrepared ?? false) return Future<void>.value();
    return _preparing[mesh] ??= _prepareAsync(mesh);
  }

  static Future<void> _prepareAsync(MeshData mesh) async {
    try {
      final prepared = await Isolate.run(
        () => prepare(mesh),
        debugName: 'pali-mesh-display',
      );
      // The worker saw a copy. All cache indices refer to that same geometry;
      // bind back to the original so selection and drawing share one cache.
      prepared._cache.mesh = mesh;
      prepared.adopt();
    } finally {
      _preparing[mesh] = null;
    }
  }

  static List<List<Vector3>?> cornerNormals(MeshData mesh) =>
      _cacheFor(mesh).cornerNormals();

  /// Coplanar, edge-connected region containing [faceIndex] (sorted face
  /// indices, includes the seed). Empty for degenerate faces.
  static List<int> regionForFace(MeshData mesh, int faceIndex) {
    if (faceIndex < 0 || faceIndex >= mesh.faces.length) return const [];
    return _cacheFor(mesh).regionForFace(faceIndex);
  }

  /// Highlight / pick / extrude region: the CAD face containing [faceIndex].
  ///
  /// A planar original-face group is one CAD face even
  /// when CSG T-junctions leave the triangles edge-disconnected. Flood that
  /// group from the seed across shared edges and T-junctions. A pinch (two
  /// islands that only meet at a vertex — inscribed circle on a square)
  /// stays separate islands, like a through-cylinder splitting a cube face
  /// into four. Tagged curved faces select every triangle of their CAD face,
  /// including fillets and periodic surfaces. Without source face IDs, curved
  /// meshes retain the local planar fallback rather than guessing boundaries.
  static List<int> selectionRegion(MeshData mesh, int faceIndex) {
    if (faceIndex < 0 || faceIndex >= mesh.faces.length) return const [];
    return _cacheFor(mesh).selectionRegion(faceIndex);
  }

  /// Edge-connected coplanar triangles plus vertex-adjacent coplanar slivers.
  ///
  /// Boolean leftovers often sit on the same plane but do not share a full
  /// edge (T-junction). Highlighting used B-Rep tags (the whole plane) while
  /// extrude walked [regionForFace] (one island) — that mismatch extruded
  /// only part of a visually selected face.
  static List<int> growPlanarRegion(MeshData mesh, int faceIndex) {
    if (faceIndex < 0 || faceIndex >= mesh.faces.length) return const [];
    return _cacheFor(mesh).growPlanarRegion(faceIndex);
  }

  /// Expands [faceIndex] into the logical opening region for shell/hollow.
  ///
  /// Flat faces use the coplanar region. Curved caps (for example after
  /// filleting) walk edge-connected faces with a similar outward orientation.
  ///
  /// Box faces are often only two triangles — do not require a larger coplanar
  /// count before accepting them, or orientation walking (dot ≥ 0.85) will
  /// swallow adjacent fillets into the opening and remove them from the shell.
  static List<int> hollowOpeningRegion(MeshData mesh, int faceIndex) {
    if (faceIndex < 0 || faceIndex >= mesh.faces.length) return const [];
    final coplanar = regionForFace(mesh, faceIndex);
    // Two or more coplanar tris ⇒ flat opening (typical tessellated face).
    if (coplanar.length >= 2) return coplanar;

    final seedOutward = _outwardNormal(mesh, faceIndex);
    if (seedOutward == null) return coplanar;

    // Lone axis-aligned triangle: still a flat opening, not a curved cap.
    final dominant = math.max(
      seedOutward.x.abs(),
      math.max(seedOutward.y.abs(), seedOutward.z.abs()),
    );
    if (dominant >= 0.95 && coplanar.isNotEmpty) return coplanar;

    final solidCenter = _meshCentroid(mesh.vertices);
    final cache = _cacheFor(mesh);
    const minDot = 0.85;
    final region = <int>{faceIndex};
    final queue = <int>[faceIndex];

    bool compatible(int f) {
      final outward = _outwardNormal(mesh, f);
      if (outward == null || outward.dot(seedOutward) < minDot) return false;
      final face = mesh.faces[f];
      var sum = Vector3.zero();
      for (final index in face) {
        sum += mesh.vertices[index];
      }
      final center = sum / face.length.toDouble();
      return (center - solidCenter).dot(seedOutward) > 0;
    }

    while (queue.isNotEmpty) {
      final f = queue.removeLast();
      final face = mesh.faces[f];
      for (var i = 0; i < face.length; i++) {
        final key = _edgeKey(
          mesh.vertices[face[i]],
          mesh.vertices[face[(i + 1) % face.length]],
        );
        for (final neighbor in cache.edgeToFaces[key] ?? const <int>[]) {
          if (region.contains(neighbor)) continue;
          if (!compatible(neighbor)) continue;
          region.add(neighbor);
          queue.add(neighbor);
        }
      }
    }

    final result = region.toList()..sort();
    return result.isEmpty ? coplanar : result;
  }

  static Vector3 _meshCentroid(List<Vector3> vertices) {
    if (vertices.isEmpty) return Vector3.zero();
    var sum = Vector3.zero();
    for (final v in vertices) {
      sum += v;
    }
    return sum / vertices.length.toDouble();
  }

  static Vector3? _outwardNormal(MeshData mesh, int faceIndex) {
    final face = mesh.faces[faceIndex];
    if (face.length < 3) return null;
    final a = mesh.vertices[face[0]];
    final b = mesh.vertices[face[1]];
    final c = mesh.vertices[face[2]];
    final n = (b - a).cross(c - a);
    if (n.length < 1e-12) return null;
    n.normalize();
    final region = regionForFace(mesh, faceIndex);
    if (region.isEmpty) return n;
    var sum = Vector3.zero();
    var count = 0;
    for (final f in region) {
      for (final index in mesh.faces[f]) {
        sum += mesh.vertices[index];
        count++;
      }
    }
    if (count == 0) return n;
    final faceCenter = sum / count.toDouble();
    final solidCenter = _meshCentroid(mesh.vertices);
    if ((faceCenter - solidCenter).dot(n) < 0) {
      return -n;
    }
    return n;
  }

  /// Cosine of the max dihedral still treated as a smooth tessellation crease
  /// (~25°). Keeps chamfer / bevel face boundaries (~45°) as hard outlines
  /// while still hiding dense cylinder / fillet tessellation bands.
  static const double smoothCreaseMinDot = 0.90;

  /// Coplanar-region tolerance (normals and plane offset). Boolean meshes
  /// wobble more than 1e-6 on a "flat" face.
  static const double coplanarNormalMinDot = 1.0 - 1e-4;
  static const double coplanarOffsetEpsilon = 1e-4;

  /// Looser plane offset used when growing a face for pick / extrude so
  /// boolean slivers on the same flat still join the parent face.
  static const double planarGrowOffsetEpsilon = 1e-3;

  /// One ID per original CAD face. Uses [MeshData.faceIds] when present;
  /// otherwise a supporting-plane flood plus a smooth-dihedral walk for
  /// curved tessellation.
  static List<int> originalFaceIds(MeshData mesh) =>
      _cacheFor(mesh).originalFaceIds();

  /// Area-weighted normal of each planar original-face group, or null when
  /// that triangle is curved / a singleton. Viewport flat-shading uses this
  /// so CSG fans on a cube face do not pick up lighting seams.
  static List<Vector3?> planarFaceNormals(MeshData mesh) =>
      _cacheFor(mesh).planarFaceNormals();

  /// Real geometric edges of the mesh, as vertex-index pairs: face-region
  /// boundaries and creases. Interior diagonals of flat faces are excluded,
  /// so wireframes show clean face outlines instead of triangle soup.
  ///
  /// Face IDs: triangles that meet across a smooth dihedral (or carry the
  /// same B-Rep [MeshData.faceIds]) are one face. Only the boundary between
  /// those faces is an outline.
  /// T-junctions left by CSG (sphere-on-face, coplanar cap cuts) are split
  /// so matching shorts share a key, then collinear Steiner verts on a
  /// two-face crease collapse back to one CAD edge (a cube's 20 mm rim).
  ///
  /// When [includeSmoothCreases] is false, leftover soft creases on cylinders,
  /// fillets, and thread flanks are omitted (draw-only; picking still uses
  /// the full set). A tagged CAD edge with a sharp section is drawn in full,
  /// including portions that become tangent near an endpoint.
  static List<(int, int)> featureEdges(
    MeshData mesh, {
    bool includeSmoothCreases = true,
  }) {
    final all = _cacheFor(mesh).featureEdges();
    if (includeSmoothCreases) return all;
    return _cacheFor(mesh).sharpFeatureEdges();
  }

  /// View-dependent outline of curved tessellation: triangle edges where
  /// one adjacent face points toward the eye and the other points away.
  ///
  /// Planar CAD faces do not contribute — they already outline via
  /// [featureEdges]. A sphere or cylinder is one face ID, so X-ray would
  /// otherwise show only rims (or nothing). [towardEye] is in mesh space.
  static List<(int, int)> silhouetteEdges(
    MeshData mesh, {
    required Vector3 Function(Vector3 point) towardEye,
  }) {
    return _cacheFor(mesh).silhouetteEdges(towardEye);
  }

  /// Count of edges that belong to only one triangle (open mesh / holes).
  static int openBoundaryEdgeCount(MeshData mesh) {
    if (mesh.faces.isEmpty) return 0;
    final all = [for (var i = 0; i < mesh.faces.length; i++) i];
    return boundaryEdges(mesh, all).length;
  }

  /// Highlight / fillet outline of [region]: CAD feature edges only.
  ///
  /// [boundaryEdges] treats any once-used segment as a border, so CSG
  /// T-junctions and fans on a holed cap show up as extra strokes when a
  /// face is highlighted. Body wireframes already hide those via
  /// [featureEdges] / original-face IDs; the selection overlay must too.
  static List<(int, int)> selectionOutline(MeshData mesh, List<int> region) {
    if (region.isEmpty) return const [];
    final boundary = boundaryEdges(mesh, region);
    if (boundary.isEmpty) return const [];
    final boundaryKeys = {
      for (final e in boundary)
        _edgeKey(mesh.vertices[e.$1], mesh.vertices[e.$2]),
    };
    final out = <(int, int)>[];
    final seen = <String>{};
    for (final e in featureEdges(mesh, includeSmoothCreases: false)) {
      final key = _edgeKey(mesh.vertices[e.$1], mesh.vertices[e.$2]);
      if (!seen.add(key)) continue;
      if (boundaryKeys.contains(key) ||
          _featureCoversBoundary(mesh, e, boundary)) {
        out.add(e);
      }
    }
    return out;
  }

  /// True when [feat] is a collapsed CAD edge whose length is mostly
  /// covered by [boundary] pieces (CSG Steiner splits of a 20 mm cube rim).
  static bool _featureCoversBoundary(
    MeshData mesh,
    (int, int) feat,
    List<(int, int)> boundary,
  ) {
    final pa = mesh.vertices[feat.$1];
    final pb = mesh.vertices[feat.$2];
    final len = (pb - pa).length;
    if (len < 1e-9) return false;
    var covered = 0.0;
    for (final e in boundary) {
      final pu = mesh.vertices[e.$1];
      final pv = mesh.vertices[e.$2];
      if (_pointOnSegment(pu, pa, pb) && _pointOnSegment(pv, pa, pb)) {
        covered += (pu - pv).length;
      }
    }
    return covered > 0.5 * len;
  }

  /// Edges of [region] that lie on its outer boundary, as vertex-index pairs.
  static List<(int, int)> boundaryEdges(MeshData mesh, List<int> region) {
    final counts = <String, (int, int)>{};
    final seen = <String, int>{};
    for (final f in region) {
      final face = mesh.faces[f];
      for (var i = 0; i < face.length; i++) {
        final a = face[i];
        final b = face[(i + 1) % face.length];
        final key = _edgeKey(mesh.vertices[a], mesh.vertices[b]);
        seen[key] = (seen[key] ?? 0) + 1;
        counts.putIfAbsent(key, () => (a, b));
      }
    }
    return [
      for (final entry in seen.entries)
        if (entry.value == 1) counts[entry.key]!,
    ];
  }

  /// Ordered closed boundary loops of [region] as vertex-index rings.
  ///
  /// Walks [boundaryEdges] by quantized position so duplicated mesh vertices
  /// still chain. Each loop lists vertices in order without repeating the
  /// start at the end. Outer and hole loops are returned separately.
  static List<List<int>> orderedBoundaryLoops(MeshData mesh, List<int> region) {
    final edges = boundaryEdges(mesh, region);
    if (edges.isEmpty) return const [];

    String posOf(int i) => quantizePosition(mesh.vertices[i]);

    final remaining = <String, (int, int)>{};
    for (final (a, b) in edges) {
      remaining[_edgeKey(mesh.vertices[a], mesh.vertices[b])] = (a, b);
    }

    final loops = <List<int>>[];
    final safetyLimit = mesh.vertices.length + 4;

    while (remaining.isNotEmpty) {
      final startEntry = remaining.entries.first;
      remaining.remove(startEntry.key);
      final (startA, startB) = startEntry.value;
      final loop = <int>[startA, startB];
      final startPos = posOf(startA);
      var currPos = posOf(startB);

      var closed = false;
      while (remaining.isNotEmpty && loop.length < safetyLimit) {
        String? usedKey;
        int? nextVert;
        String? nextPos;

        for (final entry in remaining.entries) {
          final (x, y) = entry.value;
          final px = posOf(x);
          final py = posOf(y);
          if (px == currPos) {
            usedKey = entry.key;
            nextVert = y;
            nextPos = py;
            break;
          }
          if (py == currPos) {
            usedKey = entry.key;
            nextVert = x;
            nextPos = px;
            break;
          }
        }

        if (usedKey == null || nextVert == null || nextPos == null) break;
        remaining.remove(usedKey);

        if (nextPos == startPos) {
          // Closed — do not append the start vertex again.
          closed = true;
          break;
        }
        loop.add(nextVert);
        currPos = nextPos;
      }

      if (closed && loop.length >= 3) loops.add(loop);
    }

    return loops;
  }

  static String _quantize(Vector3 v) {
    String q(double x) => (x * 1e6).roundToDouble().toString();
    return '${q(v.x)},${q(v.y)},${q(v.z)}';
  }

  static String _edgeKey(Vector3 a, Vector3 b) {
    final ka = _quantize(a);
    final kb = _quantize(b);
    return ka.compareTo(kb) <= 0 ? '$ka|$kb' : '$kb|$ka';
  }

  /// Position key for vertex deduplication (1 µm quantization).
  static String quantizePosition(Vector3 v) => _quantize(v);

  /// Undirected edge key from two positions.
  static String edgeKey(Vector3 a, Vector3 b) => _edgeKey(a, b);

  /// True when [p] lies on the closed segment [a]–[b].
  static bool pointOnSegment(Vector3 p, Vector3 a, Vector3 b) =>
      _pointOnSegment(p, a, b);

  static bool _pointOnSegment(Vector3 p, Vector3 a, Vector3 b) {
    final ab = b - a;
    final len2 = ab.length2;
    if (len2 < 1e-20) return (p - a).length2 < 1e-12;
    final t = (p - a).dot(ab) / len2;
    if (t < -1e-6 || t > 1 + 1e-6) return false;
    return (a + ab * t - p).length2 <
        planarGrowOffsetEpsilon * planarGrowOffsetEpsilon;
  }

  /// Parameter in (0, 1) when [p] lies on the open segment [a]–[b].
  static double? _tOnSegment(Vector3 p, Vector3 a, Vector3 b) {
    final ab = b - a;
    final len2 = ab.length2;
    if (len2 < 1e-20) return null;
    final t = (p - a).dot(ab) / len2;
    if (t <= 1e-6 || t >= 1 - 1e-6) return null;
    if ((a + ab * t - p).length2 >
        planarGrowOffsetEpsilon * planarGrowOffsetEpsilon) {
      return null;
    }
    return t;
  }

  /// Needle leftover: altitude / longest edge is tiny. CSG leaves these on
  /// cube edges and face diagonals; they are not CAD faces.
  static bool isSliverFace(MeshData mesh, int faceIndex) {
    if (faceIndex < 0 || faceIndex >= mesh.faces.length) return true;
    return _cacheFor(mesh)._isSliver(faceIndex);
  }

  /// Triangle indices adjacent to the undirected edge [a]–[b].
  static List<int> facesForEdge(MeshData mesh, int a, int b) {
    if (a < 0 ||
        b < 0 ||
        a >= mesh.vertices.length ||
        b >= mesh.vertices.length) {
      return const [];
    }
    final key = _edgeKey(mesh.vertices[a], mesh.vertices[b]);
    return _cacheFor(mesh).edgeToFaces[key] ?? const [];
  }

  /// Triangles that meet at the CAD crease [a]–[b], including Steiner
  /// fragments of a collapsed feature.
  ///
  /// [facesForEdge] only sees a real mesh edge. After CSG, a cube rim is
  /// often one feature whose endpoints no longer share a triangle — the
  /// two supporting planes still sit on the short pieces in between.
  static List<int> facesAlongCrease(MeshData mesh, int a, int b) {
    if (a < 0 ||
        b < 0 ||
        a >= mesh.vertices.length ||
        b >= mesh.vertices.length) {
      return const [];
    }
    final out = <int>{...facesForEdge(mesh, a, b)};
    if (out.length >= 2) return out.toList();
    final pa = mesh.vertices[a];
    final pb = mesh.vertices[b];
    if ((pb - pa).length2 < 1e-20) return out.toList();
    for (var f = 0; f < mesh.faces.length; f++) {
      if (out.contains(f)) continue;
      final face = mesh.faces[f];
      if (face.length < 3) continue;
      for (var i = 0; i < face.length; i++) {
        final ua = mesh.vertices[face[i]];
        final ub = mesh.vertices[face[(i + 1) % face.length]];
        if ((ua - ub).length2 < 1e-20) continue;
        if (_pointOnSegment(ua, pa, pb) && _pointOnSegment(ub, pa, pb)) {
          out.add(f);
          break;
        }
      }
    }
    return out.toList();
  }
}

class _MeshRegionCache {
  _MeshRegionCache(this.mesh);

  MeshData mesh;
  bool displayPrepared = false;
  final Map<int, List<int>> _regionByFace = {};
  final Map<int, List<int>> _grownByFace = {};
  final Map<int, List<int>> _selectionByFace = {};
  final Map<int, List<int>> _curvedTaggedFaces = {};
  Map<String, List<int>>? _edgeToFaces;
  List<(int, int)>? _featureEdges;
  List<(int, int)>? _sharpFeatureEdges;
  List<int>? _originalFaceIds;
  List<Vector3?>? _planarFaceNormals;

  List<List<Vector3>?>? _cornerNormals;

  List<List<Vector3>?> cornerNormals() {
    final cached = _cornerNormals;
    if (cached != null) return cached;
    final faceCount = mesh.faces.length;
    final faceNormals = List<Vector3?>.filled(faceCount, null);
    final atVertex = List<List<int>>.generate(
      mesh.vertices.length,
      (_) => <int>[],
    );
    for (var f = 0; f < faceCount; f++) {
      final face = mesh.faces[f];
      if (face.length < 3) continue;
      final a = mesh.vertices[face[0]];
      final b = mesh.vertices[face[1]];
      final c = mesh.vertices[face[2]];
      final n = (b - a).cross(c - a);
      if (n.length < 1e-12) continue;
      n.normalize();
      faceNormals[f] = n;
      for (final vi in face) {
        if (vi >= 0 && vi < atVertex.length) atVertex[vi].add(f);
      }
    }

    const minDot = MeshFaceRegions.smoothCreaseMinDot;
    final out = List<List<Vector3>?>.filled(faceCount, null);
    for (var f = 0; f < faceCount; f++) {
      final fn = faceNormals[f];
      final face = mesh.faces[f];
      if (fn == null || face.length < 3) continue;
      final corners = <Vector3>[];
      for (final vi in face) {
        var sum = Vector3.zero();
        for (final g in atVertex[vi]) {
          final gn = faceNormals[g];
          if (gn == null) continue;
          if (gn.dot(fn) >= minDot) sum += gn;
        }
        if (sum.length2 < 1e-20) {
          corners.add(fn.clone());
        } else {
          corners.add(sum.normalized());
        }
      }
      out[f] = corners;
    }
    return _cornerNormals = out;
  }

  Map<String, List<int>> get edgeToFaces {
    final existing = _edgeToFaces;
    if (existing != null) return existing;
    final map = <String, List<int>>{};
    for (var f = 0; f < mesh.faces.length; f++) {
      final face = mesh.faces[f];
      for (var i = 0; i < face.length; i++) {
        final key = MeshFaceRegions._edgeKey(
          mesh.vertices[face[i]],
          mesh.vertices[face[(i + 1) % face.length]],
        );
        (map[key] ??= []).add(f);
      }
    }
    return _edgeToFaces = map;
  }

  Vector3? _faceNormal(int f) {
    final face = mesh.faces[f];
    if (face.length < 3) return null;
    final a = mesh.vertices[face[0]];
    final b = mesh.vertices[face[1]];
    final c = mesh.vertices[face[2]];
    final n = (b - a).cross(c - a);
    if (n.length < 1e-12) return null;
    return n.normalized();
  }

  bool _bothFacesOnOnePlane(int f0, int f1) {
    bool withNormal(Vector3 n, int originFace) {
      final p = mesh.vertices[mesh.faces[originFace][0]];
      bool onPlane(int f) {
        for (final vi in mesh.faces[f]) {
          if ((mesh.vertices[vi] - p).dot(n).abs() >=
              MeshFaceRegions.coplanarOffsetEpsilon) {
            return false;
          }
        }
        return true;
      }

      return onPlane(f0) && onPlane(f1);
    }

    // Prefer a non-degenerate host. A needle sliver's noisy normal would
    // otherwise reject the wall it sits on.
    final n0 = _faceNormal(f0);
    if (n0 != null && withNormal(n0, f0)) return true;
    final n1 = _faceNormal(f1);
    if (n1 != null && withNormal(n1, f1)) return true;
    return false;
  }

  List<int> regionForFace(int seedFace) {
    final cached = _regionByFace[seedFace];
    if (cached != null) return cached;

    final seedNormal = _faceNormal(seedFace);
    if (seedNormal == null) {
      return _regionByFace[seedFace] = const [];
    }
    final seedPoint = mesh.vertices[mesh.faces[seedFace][0]];

    bool coplanar(int f) {
      final n = _faceNormal(f);
      if (n == null) return false;
      // Slightly tolerant: boolean meshes often wobble normals on a
      // "flat" face enough that 1e-6 rejects siblings and leaves one triangle
      // (whose hypotenuse then becomes the sketch X axis).
      if (n.dot(seedNormal).abs() < MeshFaceRegions.coplanarNormalMinDot) {
        return false;
      }
      final d = (mesh.vertices[mesh.faces[f][0]] - seedPoint).dot(seedNormal);
      return d.abs() < MeshFaceRegions.coplanarOffsetEpsilon;
    }

    final region = <int>{seedFace};
    final queue = <int>[seedFace];
    while (queue.isNotEmpty) {
      final f = queue.removeLast();
      final face = mesh.faces[f];
      for (var i = 0; i < face.length; i++) {
        final key = MeshFaceRegions._edgeKey(
          mesh.vertices[face[i]],
          mesh.vertices[face[(i + 1) % face.length]],
        );
        for (final neighbor in edgeToFaces[key] ?? const <int>[]) {
          if (region.contains(neighbor)) continue;
          if (!coplanar(neighbor)) continue;
          region.add(neighbor);
          queue.add(neighbor);
        }
      }
    }

    final result = region.toList()..sort();
    for (final f in result) {
      _regionByFace[f] = result;
    }
    return result;
  }

  List<int> growPlanarRegion(int seedFace) {
    final cached = _grownByFace[seedFace];
    if (cached != null) return cached;

    final seed = regionForFace(seedFace);
    final region = <int>{if (seed.isNotEmpty) ...seed else seedFace};

    final seedNormal = _faceNormal(seedFace);
    final seedFaceVerts = mesh.faces[seedFace];
    if (seedNormal == null || seedFaceVerts.length < 3) {
      final out = region.toList()..sort();
      for (final f in out) {
        _grownByFace[f] = out;
      }
      return out;
    }
    final seedPoint = mesh.vertices[seedFaceVerts[0]];

    bool coplanar(int f) {
      final n = _faceNormal(f);
      if (n == null) return false;
      if (n.dot(seedNormal).abs() < MeshFaceRegions.coplanarNormalMinDot) {
        return false;
      }
      for (final vi in mesh.faces[f]) {
        if ((mesh.vertices[vi] - seedPoint).dot(seedNormal).abs() >=
            MeshFaceRegions.planarGrowOffsetEpsilon) {
          return false;
        }
      }
      return true;
    }

    String posOf(int i) => MeshFaceRegions.quantizePosition(mesh.vertices[i]);
    final verts = <String>{};
    for (final f in region) {
      for (final vi in mesh.faces[f]) {
        verts.add(posOf(vi));
      }
    }

    bool onRegion(int f) {
      var shared = 0;
      for (final vi in mesh.faces[f]) {
        if (verts.contains(posOf(vi))) {
          shared++;
          continue;
        }
        final p = mesh.vertices[vi];
        for (final rf in region) {
          final rfFace = mesh.faces[rf];
          for (var i = 0; i < rfFace.length; i++) {
            if (MeshFaceRegions._pointOnSegment(
              p,
              mesh.vertices[rfFace[i]],
              mesh.vertices[rfFace[(i + 1) % rfFace.length]],
            )) {
              shared++;
              break;
            }
          }
        }
      }
      final face = mesh.faces[f];
      if (face.length < 3) return false;
      final a = mesh.vertices[face[1]] - mesh.vertices[face[0]];
      final b = mesh.vertices[face[2]] - mesh.vertices[face[0]];
      final sliver = 0.5 * a.cross(b).length < 1e-4;
      return shared >= (sliver ? 1 : 2);
    }

    var changed = true;
    while (changed) {
      changed = false;
      for (var f = 0; f < mesh.faces.length; f++) {
        if (region.contains(f)) continue;
        if (!coplanar(f)) continue;
        if (!onRegion(f)) continue;
        region.add(f);
        for (final vi in mesh.faces[f]) {
          verts.add(posOf(vi));
        }
        changed = true;
      }
    }

    final result = region.toList()..sort();
    for (final f in result) {
      _grownByFace[f] = result;
    }
    return result;
  }

  List<int> selectionRegion(int faceIndex) {
    final cached = _selectionByFace[faceIndex];
    if (cached != null) return cached;

    final planar = planarFaceNormals();
    final stored = mesh.faceIds;
    if (stored != null && stored.length == mesh.faces.length) {
      final curved = _curvedTaggedFaces[stored[faceIndex]];
      if (curved != null) return curved;
    }

    var seed = faceIndex;
    if (_isSliver(seed)) {
      final host = _hostFaceForSliver(seed);
      if (host == null) return growPlanarRegion(faceIndex);
      seed = host;
    }

    if (planar[seed] == null) {
      return growPlanarRegion(seed);
    }

    final ids = originalFaceIds();
    final id = ids[seed];
    if (id < 0) return growPlanarRegion(seed);

    final seedN = planar[seed]!;
    final tagged =
        stored != null &&
        stored.length == mesh.faces.length &&
        stored[seed] >= 0;
    final group = <int>[];
    for (var f = 0; f < mesh.faces.length; f++) {
      if (ids[f] != id) continue;
      // Native cap tessellation needs thin triangles around text and curved
      // holes. Their CAD face tag is authoritative; dropping them cuts gaps
      // through the selection fill and can disconnect the rest of the face.
      final sliver = _isSliver(f);
      if (sliver && !tagged) continue;
      final n = _faceNormal(f);
      if (!sliver && n != null && n.dot(seedN) < 0) continue;
      group.add(f);
    }
    if (group.length <= 1) return growPlanarRegion(seed);

    final groupSet = group.toSet();
    final seen = <int>{seed};
    final queue = [seed];
    while (queue.isNotEmpty) {
      final f = queue.removeLast();
      final face = mesh.faces[f];
      for (var i = 0; i < face.length; i++) {
        final key = MeshFaceRegions._edgeKey(
          mesh.vertices[face[i]],
          mesh.vertices[face[(i + 1) % face.length]],
        );
        for (final n in edgeToFaces[key] ?? const <int>[]) {
          if (!groupSet.contains(n) || !seen.add(n)) continue;
          queue.add(n);
        }
      }
      for (final n in group) {
        if (seen.contains(n)) continue;
        if (!_tJunctionLinked(f, n)) continue;
        seen.add(n);
        queue.add(n);
      }
    }

    final result = seen.toList()..sort();
    for (final f in result) {
      _selectionByFace[f] = result;
    }
    if (faceIndex != seed) {
      _selectionByFace[faceIndex] = result;
    }
    return result;
  }

  /// Needle leftover on a CAD face: the largest real triangle that owns
  /// this sliver's original-face ID, else a unique coplanar neighbor.
  int? _hostFaceForSliver(int sliver) {
    final ids = originalFaceIds();
    final id = ids[sliver];
    int? best;
    var bestArea = -1.0;
    void consider(int f) {
      if (f == sliver || _isSliver(f)) return;
      final a = _area2(f);
      if (a > bestArea) {
        bestArea = a;
        best = f;
      }
    }

    if (id >= 0) {
      for (var f = 0; f < mesh.faces.length; f++) {
        if (ids[f] == id) consider(f);
      }
      if (best != null) return best;
    }
    int? host;
    var conflict = false;
    final face = mesh.faces[sliver];
    for (var i = 0; i < face.length; i++) {
      final key = MeshFaceRegions._edgeKey(
        mesh.vertices[face[i]],
        mesh.vertices[face[(i + 1) % face.length]],
      );
      for (final nbr in edgeToFaces[key] ?? const <int>[]) {
        if (_isSliver(nbr)) continue;
        if (!_bothFacesOnOnePlane(sliver, nbr)) continue;
        if (host == null) {
          host = nbr;
        } else if (host != nbr && ids[host] != ids[nbr]) {
          conflict = true;
        }
      }
    }
    return conflict ? null : host;
  }

  /// Vertex of one triangle lies in the open interior of an edge of the
  /// other. Reunites CSG Steiner islands. A pinch that only shares an
  /// endpoint (inscribed circle on a square) is not a T-junction.
  bool _tJunctionLinked(int a, int b) {
    bool vertOnOpenEdge(int vertFace, int edgeFace) {
      for (final vi in mesh.faces[vertFace]) {
        final p = mesh.vertices[vi];
        final ef = mesh.faces[edgeFace];
        for (var i = 0; i < ef.length; i++) {
          if (MeshFaceRegions._tOnSegment(
                p,
                mesh.vertices[ef[i]],
                mesh.vertices[ef[(i + 1) % ef.length]],
              ) !=
              null) {
            return true;
          }
        }
      }
      return false;
    }

    return vertOnOpenEdge(a, b) || vertOnOpenEdge(b, a);
  }

  double _area2(int f) {
    final face = mesh.faces[f];
    if (face.length < 3) return 0;
    final a = mesh.vertices[face[0]];
    final b = mesh.vertices[face[1]];
    final c = mesh.vertices[face[2]];
    return (b - a).cross(c - a).length2;
  }

  bool _vertsOnPlane(int f, Vector3 origin, Vector3 normal, double eps) {
    for (final vi in mesh.faces[f]) {
      if ((mesh.vertices[vi] - origin).dot(normal).abs() >= eps) {
        return false;
      }
    }
    return true;
  }

  bool _isSliver(int f) {
    final face = mesh.faces[f];
    if (face.length < 3) return true;
    final ab = mesh.vertices[face[1]] - mesh.vertices[face[0]];
    final ac = mesh.vertices[face[2]] - mesh.vertices[face[0]];
    final longest2 = math.max(ab.length2, ac.length2);
    final bc = mesh.vertices[face[2]] - mesh.vertices[face[1]];
    final long2 = math.max(longest2, bc.length2);
    if (long2 < 1e-20) return true;
    // Altitude / longest edge < 1e-4 → needle leftover, not a CAD face.
    return ab.cross(ac).length2 < long2 * 1e-8;
  }

  /// One ID per original face: stored [MeshData.faceIds],
  /// or a supporting-plane flood plus a walk across
  /// edges whose dihedral is smoother than [smoothCreaseMinDot].
  ///
  /// Quantized plane keys fragment a CSG-cut cube face: intersection verts
  /// sit slightly off-plane, fan triangles get noisy / degenerate normals,
  /// and every radius draws as a feature. Split children of one input face
  /// stay on one ID, then coplanar tris group by supporting plane — every
  /// vertex on the plane, not matching quantized normals.
  /// A sphere / cylinder triangle has a vertex off the plane and stays out.
  /// Crease slivers on an edge lie on two planes; the largest host claims
  /// them. The rim still draws because [featureEdges] keeps edges where a
  /// second CAD face meets them. Do not union two hosts
  /// through a sliver.
  List<int> originalFaceIds() {
    final stored = mesh.faceIds;
    if (stored != null && stored.length == mesh.faces.length) {
      var tagged = false;
      for (final id in stored) {
        if (id >= 0) {
          tagged = true;
          break;
        }
      }
      // All-negative IDs mean "untagged", not one CAD face. Feature
      // edges skip id < 0, so treating -1 as stored would hide every
      // crease on imported / legacy meshes.
      if (tagged) return stored;
    }
    final cached = _originalFaceIds;
    if (cached != null) return cached;
    final n = mesh.faces.length;
    final ids = List<int>.filled(n, -1);
    var next = 0;

    // Largest reliable-normal triangles first so a cube face claims its
    // fans and needles before a tiny leftover can become its own plane.
    const planeEps = MeshFaceRegions.planarGrowOffsetEpsilon;
    final seeds = <int>[
      for (var f = 0; f < n; f++)
        if (_faceNormal(f) != null && _area2(f) > 1e-16) f,
    ]..sort((a, b) => _area2(b).compareTo(_area2(a)));

    final candidates = <_PlaneCandidate>[];
    for (final f in seeds) {
      if (_isSliver(f)) continue;
      final face = mesh.faces[f];
      final a = mesh.vertices[face[0]],
          b = mesh.vertices[face[1]],
          c = mesh.vertices[face[2]];
      final ab = b - a, ac = c - a, bc = c - b;
      final cross = ab.cross(ac), length = cross.length;
      final longest = math.sqrt(
        math.max(ab.length2, math.max(ac.length2, bc.length2)),
      );
      // If every vertex is within epsilon of a plane, its normal lies within
      // 8*epsilon/altitude of one of the triangle's two unit normals. This
      // conservative bound admits small/oblique triangles too; the original
      // per-vertex predicate below still decides membership exactly.
      final radius = math.min(2.0, 8 * planeEps * longest / length + 1e-10);
      candidates.add((
        id: f,
        center:
            face.fold(Vector3.zero(), (sum, vi) => sum + mesh.vertices[vi]) /
            face.length.toDouble(),
        normal: cross / length,
        radius: radius,
      ));
    }
    // Tiny triangles have wide normal bounds. Keep them out of the same
    // hierarchy as well-conditioned faces, otherwise a few needles prevent
    // pruning almost every node in a large curved mesh.
    final buckets = <int, List<_PlaneCandidate>>{};
    for (final candidate in candidates) {
      final bucket = (math.log(candidate.radius) / math.ln2).floor();
      (buckets[bucket] ??= []).add(candidate);
    }
    final planeIndices = [
      for (final bucket in buckets.values) _FacePlaneIndex(bucket),
    ];

    for (final seed in seeds) {
      if (ids[seed] >= 0) continue;
      if (_isSliver(seed)) continue;
      final nrm = _faceNormal(seed);
      if (nrm == null) continue;
      final origin = mesh.vertices[mesh.faces[seed][0]];
      final id = next++;
      void assign(int f) {
        if (ids[f] >= 0 || _isSliver(f)) return;
        if (_vertsOnPlane(f, origin, nrm, planeEps)) ids[f] = id;
      }

      for (final index in planeIndices) {
        index.visitNearPlane(origin, nrm, planeEps * 1.000001, ids, assign);
      }
    }

    for (var seed = 0; seed < n; seed++) {
      if (ids[seed] >= 0) continue;
      if (_isSliver(seed)) continue;
      final id = next++;
      ids[seed] = id;
      final stack = [seed];
      while (stack.isNotEmpty) {
        final f = stack.removeLast();
        final nF = _faceNormal(f);
        final face = mesh.faces[f];
        if (nF == null || face.length < 3) continue;
        for (var i = 0; i < face.length; i++) {
          final key = MeshFaceRegions._edgeKey(
            mesh.vertices[face[i]],
            mesh.vertices[face[(i + 1) % face.length]],
          );
          for (final nbr in edgeToFaces[key] ?? const <int>[]) {
            if (ids[nbr] >= 0) continue;
            if (_isSliver(nbr)) continue;
            final nn = _faceNormal(nbr);
            if (nn == null) continue;
            if (nn.dot(nF) < MeshFaceRegions.smoothCreaseMinDot) continue;
            ids[nbr] = id;
            stack.add(nbr);
          }
        }
      }
    }

    var changed = true;
    while (changed) {
      changed = false;
      for (var f = 0; f < n; f++) {
        if (!_isSliver(f)) continue;
        int? hostId;
        var conflict = false;
        final face = mesh.faces[f];
        for (var i = 0; i < face.length; i++) {
          final key = MeshFaceRegions._edgeKey(
            mesh.vertices[face[i]],
            mesh.vertices[face[(i + 1) % face.length]],
          );
          for (final nbr in edgeToFaces[key] ?? const <int>[]) {
            if (_isSliver(nbr)) continue;
            if (!_bothFacesOnOnePlane(f, nbr)) continue;
            final id = ids[nbr];
            if (id < 0) continue;
            if (hostId == null) {
              hostId = id;
            } else if (id != hostId) {
              conflict = true;
            }
          }
        }
        if (!conflict && hostId != null && ids[f] != hostId) {
          ids[f] = hostId;
          changed = true;
        }
      }
    }
    return _originalFaceIds = ids;
  }

  List<Vector3?> planarFaceNormals() {
    final cached = _planarFaceNormals;
    if (cached != null) return cached;
    final ids = originalFaceIds();
    final n = mesh.faces.length;
    final out = List<Vector3?>.filled(n, null);
    final byId = <int, List<int>>{};
    for (var f = 0; f < n; f++) {
      (byId[ids[f]] ??= []).add(f);
    }
    for (final faces in byId.values) {
      final reliable = [
        for (final f in faces)
          if (!_isSliver(f)) f,
      ];
      if (reliable.length < 2) continue;
      var nSum = Vector3.zero();
      Vector3? origin;
      for (final f in reliable) {
        final face = mesh.faces[f];
        if (face.length < 3) continue;
        nSum += (mesh.vertices[face[1]] - mesh.vertices[face[0]]).cross(
          mesh.vertices[face[2]] - mesh.vertices[face[0]],
        );
        origin ??= mesh.vertices[face[0]];
      }
      if (nSum.length2 < 1e-20 || origin == null) continue;
      final nHat = nSum.normalized();
      var planar = true;
      for (final f in reliable) {
        final fn = _faceNormal(f);
        if (fn != null &&
            fn.dot(nHat).abs() < MeshFaceRegions.smoothCreaseMinDot) {
          planar = false;
          break;
        }
        if (!_vertsOnPlane(
          f,
          origin,
          nHat,
          MeshFaceRegions.planarGrowOffsetEpsilon,
        )) {
          planar = false;
          break;
        }
      }
      if (!planar) continue;
      for (final f in faces) {
        out[f] = nHat;
      }
    }
    // Reuse the groups already built during worker display preparation. A
    // curved face's topology tag is authoritative even across periodic seams
    // or tiny display slivers; a planar flood would select tessellation bands
    // and can even leak into a tangent neighbor. Untagged inferred groups are
    // deliberately excluded. Subsequent hover/pick lookup is constant-time.
    final stored = mesh.faceIds;
    if (stored != null && stored.length == n) {
      for (final entry in byId.entries) {
        if (entry.key >= 0 &&
            stored[entry.value.first] == entry.key &&
            out[entry.value.first] == null) {
          _curvedTaggedFaces[entry.key] = entry.value;
        }
      }
    }
    return _planarFaceNormals = out;
  }

  List<(int, int)> featureEdges() {
    final cached = _featureEdges;
    if (cached != null) return cached;

    final ids = originalFaceIds();
    // T-split every triangle edge at any *used* mesh vertex on it.
    // Unify / weld leave orphan mid-edge leftovers; treating those as
    // Steiner points splits a cube's 20mm side into two 10mm strokes.
    final usedVert = List<bool>.filled(mesh.vertices.length, false);
    for (final face in mesh.faces) {
      for (final i in face) {
        if (i >= 0 && i < usedVert.length) usedVert[i] = true;
      }
    }
    final points = <(int, Vector3)>[];
    final seenPos = <String>{};
    for (var vi = 0; vi < mesh.vertices.length; vi++) {
      if (!usedVert[vi]) continue;
      if (seenPos.add(MeshFaceRegions.quantizePosition(mesh.vertices[vi]))) {
        points.add((vi, mesh.vertices[vi]));
      }
    }

    final pointIndex = _SegmentPointIndex(points);
    final idSets = <String, Set<int>>{};
    final reps = <String, (int, int)>{};
    final uses = <String, int>{};

    void addSeg(int a, int b, int id) {
      if (a == b) return;
      final key = MeshFaceRegions._edgeKey(mesh.vertices[a], mesh.vertices[b]);
      (idSets[key] ??= {}).add(id);
      uses[key] = (uses[key] ?? 0) + 1;
      reps[key] = (a, b);
    }

    for (var f = 0; f < mesh.faces.length; f++) {
      final face = mesh.faces[f];
      if (face.length < 3) continue;
      final id = ids[f];
      if (id < 0) continue;
      for (var i = 0; i < face.length; i++) {
        final a = face[i];
        final b = face[(i + 1) % face.length];
        final pa = mesh.vertices[a];
        final pb = mesh.vertices[b];
        final mids = <(double, int)>[];
        for (final (vi, p) in pointIndex.nearSegment(pa, pb)) {
          if (vi == a || vi == b) continue;
          final t = MeshFaceRegions._tOnSegment(p, pa, pb);
          if (t != null) mids.add((t, vi));
        }
        mids.sort((x, y) => x.$1.compareTo(y.$1));
        var prev = a;
        var prevT = 0.0;
        for (final m in mids) {
          if (m.$1 - prevT < 1e-6) continue;
          addSeg(prev, m.$2, id);
          prev = m.$2;
          prevT = m.$1;
        }
        addSeg(prev, b, id);
      }
    }

    final result = <(int, int)>[];
    for (final key in idSets.keys) {
      final nIds = idSets[key]!.length;
      final nUses = uses[key]!;
      if (nIds == 1 && nUses != 1) continue;
      result.add(reps[key]!);
    }

    final vertFaceIds = List<Set<int>>.generate(
      mesh.vertices.length,
      (_) => <int>{},
    );
    for (var f = 0; f < mesh.faces.length; f++) {
      final id = ids[f];
      if (id < 0) continue;
      for (final vi in mesh.faces[f]) {
        if (vi >= 0 && vi < vertFaceIds.length) vertFaceIds[vi].add(id);
      }
    }
    return _featureEdges = _collapseColinearFeatures(result, vertFaceIds);
  }

  /// A vertex with only two CAD faces
  /// that sits between two collinear feature strokes is Steiner (infinite-
  /// plane CSG nicking a cube rim). Cube corners and a sphere hitting a
  /// rim have ≥3 face IDs and stay split.
  List<(int, int)> _collapseColinearFeatures(
    List<(int, int)> edges,
    List<Set<int>> vertFaceIds,
  ) {
    final segments = <int, (int, int)>{};
    final byPair = <(int, int), int>{};
    final adjacent = <int, Set<int>>{};
    (int, int) pair(int a, int b) => a < b ? (a, b) : (b, a);
    var nextId = 0;
    void insert(int a, int b) {
      final id = nextId++;
      segments[id] = (a, b);
      byPair[pair(a, b)] = id;
      (adjacent[a] ??= {}).add(id);
      (adjacent[b] ??= {}).add(id);
    }

    void remove(int id) {
      final e = segments.remove(id)!;
      byPair.remove(pair(e.$1, e.$2));
      adjacent[e.$1]!.remove(id);
      adjacent[e.$2]!.remove(id);
    }

    for (final e in edges) insert(e.$1, e.$2);
    final ready = SplayTreeMap<int, int>(), priorities = <int, int>{};
    void refresh(int v) {
      final old = priorities.remove(v);
      if (old != null) ready.remove(old);
      final uses = adjacent[v]!;
      if (uses.length != 2 ||
          v < 0 ||
          v >= vertFaceIds.length ||
          vertFaceIds[v].length >= 3)
        return;
      final first = uses.reduce(math.min), last = uses.reduce(math.max);
      final e = segments[first]!, f = segments[last]!;
      final a = e.$1 == v ? e.$2 : e.$1, b = f.$1 == v ? f.$2 : f.$1;
      if (a == b ||
          !_collinearTriple(
            mesh.vertices[a],
            mesh.vertices[v],
            mesh.vertices[b],
          ))
        return;
      // Match the former scan order: the first surviving segment containing
      // this vertex, with its first endpoint visited before its second.
      final priority = first * 2 + (e.$1 == v ? 0 : 1);
      priorities[v] = priority;
      ready[priority] = v;
    }

    for (final v in adjacent.keys) refresh(v);
    while (ready.isNotEmpty) {
      final priority = ready.firstKey()!, v = ready.remove(priority)!;
      priorities.remove(v);
      final ids = adjacent[v]!.toList();
      final e = segments[ids[0]]!, f = segments[ids[1]]!;
      final a = e.$1 == v ? e.$2 : e.$1, b = f.$1 == v ? f.$2 : f.$1;
      remove(ids[0]);
      remove(ids[1]);
      if (!byPair.containsKey(pair(a, b))) insert(a, b);
      refresh(v);
      refresh(a);
      refresh(b);
    }
    return segments.values.toList();
  }

  bool _collinearTriple(Vector3 a, Vector3 m, Vector3 b) {
    final ab = b - a;
    final len2 = ab.length2;
    if (len2 < 1e-20) return true;
    final t = (m - a).dot(ab) / len2;
    if (t <= 1e-6 || t >= 1 - 1e-6) return false;
    final dist2 = (a + ab * t - m).length2;
    final eps = MeshFaceRegions.planarGrowOffsetEpsilon;
    return dist2 < eps * eps;
  }

  /// [featureEdges] minus soft tessellation creases (cylinder / fillet / screw).
  List<(int, int)> sharpFeatureEdges() {
    final cached = _sharpFeatureEdges;
    if (cached != null) return cached;

    final result = <(int, int)>[];
    for (final (a, b) in featureEdges()) {
      final key = MeshFaceRegions._edgeKey(mesh.vertices[a], mesh.vertices[b]);
      final adjacent = edgeToFaces[key] ?? const <int>[];
      if (adjacent.length != 2) {
        result.add((a, b));
        continue;
      }
      final n1 = _faceNormal(adjacent[0]);
      final n2 = _faceNormal(adjacent[1]);
      if (n1 == null || n2 == null) {
        result.add((a, b));
        continue;
      }
      // Keep only sharp creases and face-region boundaries.
      if (n1.dot(n2) < MeshFaceRegions.smoothCreaseMinDot) {
        result.add((a, b));
      }
    }
    // A single CAD edge may approach tangency at its ends, as the elliptical
    // seams between inverse fillets do. Classifying each display chord alone
    // cuts those outlines short. Keep all chords of an edge once any part is
    // sharp; fully smooth joins and untagged mesh edges retain their filtering.
    final sharpIds = <int>{
      for (final (a, b) in result)
        if (mesh.brepEdgeOfVertices(a, b) case final int id) id,
    };
    if (sharpIds.isEmpty) return _sharpFeatureEdges = result;
    final sharpSegments = result.toSet();
    return _sharpFeatureEdges = [
      for (final edge in featureEdges())
        if (sharpSegments.contains(edge) ||
            sharpIds.contains(mesh.brepEdgeOfVertices(edge.$1, edge.$2)))
          edge,
    ];
  }

  List<(int, int, Vector3, Vector3, Vector3)>? _silhouetteCandidates;

  List<(int, int, Vector3, Vector3, Vector3)> silhouetteCandidates() {
    final cached = _silhouetteCandidates;
    if (cached != null) return cached;
    final features = <String>{
      for (final e in featureEdges())
        MeshFaceRegions._edgeKey(mesh.vertices[e.$1], mesh.vertices[e.$2]),
    };
    final seen = <String>{};
    final out = <(int, int, Vector3, Vector3, Vector3)>[];
    for (var f = 0; f < mesh.faces.length; f++) {
      final face = mesh.faces[f];
      if (face.length < 3) continue;
      for (var i = 0; i < face.length; i++) {
        final a = face[i];
        final b = face[(i + 1) % face.length];
        if (a == b) continue;
        final key = MeshFaceRegions._edgeKey(
          mesh.vertices[a],
          mesh.vertices[b],
        );
        if (!seen.add(key) || features.contains(key)) continue;
        final adj = edgeToFaces[key];
        if (adj == null || adj.length != 2) continue;
        final n1 = _faceNormal(adj[0]);
        final n2 = _faceNormal(adj[1]);
        if (n1 == null || n2 == null) continue;
        final mid = (mesh.vertices[a] + mesh.vertices[b]) * 0.5;
        out.add((a, b, mid, n1, n2));
      }
    }
    return _silhouetteCandidates = out;
  }

  List<(int, int)> silhouetteEdges(Vector3 Function(Vector3 point) towardEye) {
    final out = <(int, int)>[];
    for (final (a, b, mid, n1, n2) in silhouetteCandidates()) {
      final view = towardEye(mid);
      if (view.length2 >= 1e-20 && n1.dot(view) * n2.dot(view) <= 0) {
        out.add((a, b));
      }
    }
    return out;
  }
}
