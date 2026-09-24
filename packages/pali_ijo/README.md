# pali ijo
PLACEHOLDER DIRECTORY TO BE UPLAODED SOON

A Dart BREP CAD kernel.

**The BREP is the solid.** Curves, surfaces and their shared boundaries are stored independently of the triangles used for display. Supported modelling operations preserve that representation. An unsupported BREP operation throws `BrepUnsupported`; it must not silently convert the document to a mesh.

Imported meshes and explicitly converted bodies have a separate mesh representation. Threads preserve native BREP inputs; legacy mesh APIs continue to produce mesh solids. BREP support is substantial but is not yet general for every surface intersection or offset; see the limits below.

```dart
import 'package:pali_ijo/pali_ijo.dart';
import 'package:vector_math/vector_math_64.dart';

final block = BrepPrimitives.box(30, 20, 10);
final hole = BrepPrimitives.cylinder(
  origin: Vector3(15, 10, -1),
  axis: Vector3(0, 0, 1),
  radius: 4,
  height: 12,
);
final part = IjoOps.booleanSolid(block, hole, IjoOps.subtractOp);
final display = IjoOps.tessellateClosed(part, op: 'Bore');

print(IjoMeasure.ofSolid(part).volume);
print(MeshKit.isClosedSolid(display));
final saved = IjoDocument.solid(part).encode();
final restored = IjoDocument.decode(saved).solid;
```

## Use from Dart, the CLI, or the browser

Import `package:pali_ijo/pali_ijo.dart` for the Dart library. The kernel has no Flutter or FFI dependency. PaliCAD uses this API directly.

The CLI reads a document from stdin and writes a document to stdout. A solid document retains BREP geometry through supported commands:

```bash
tool/build_cli.sh

pali_ijo box --size 40,40,20 \
  | pali_ijo hollow --thickness 2 --faces top \
  | pali_ijo export --out part.3mf
```

Use `pali_ijo -v boolean --op subtract --with tool.json`, or `IJO_DEBUG=1`, for stderr breadcrumbs. `IJO_DEBUG_FILE=path` also appends them to a file. These messages identify the last recorded operation, rather than guaranteeing the cause of a crash.

`tool/build_wasm.sh` builds the WasmGC command surface. See [`wasm/pali_ijo_wasm.dart`](wasm/pali_ijo_wasm.dart).

## Persistence and export

Use `IjoDocument.solid(solid).encode()` and `IjoDocument.decode(text)` to preserve topology, surface parameters and provenance. `BrepSolid.toJson` / `fromJson` are also available. A document's `.solid` accesses its BREP; `.mesh` generates or accesses a display mesh. Saving only that mesh discards the analytic representation.

Affine geometry, unwrapped circular/conical helices, ruled surfaces, rational B-splines, implicit plane/quadric sections and general surface-intersection branches have native JSON representations. Readers predating these geometry types cannot read documents containing them.

STL, OBJ and 3MF exports are tessellations by definition. Use `IjoStep.writeSolid` and `IjoStep.readSolid` / `readSolids` for STEP BREP interchange. `StepPart.solid` supports multiple BREP parts. PaliCAD exports the cached BREP with its world placement and imports both the BREP and a derived display mesh.

STEP supports native planes, cylinders, cones, spheres and tori, plus rational B-spline curves and surfaces. Conics and affine quadrics have exact rational representations. Compatible rational ruled walls and their plane sections now convert algebraically, without fitting. This includes differently weighted spline boundaries and matching angular conics. STEP also retains cylinder/quadric and general spline-intersection parent-surface references and recovers their native implicit curves on supported imports. General branch recovery verifies geometric agreement within bounded parameter intervals. Other bounded intersection curves and incompatible ruled parameterizations still use adaptive cubic fitting with a checked target of `1e-7` mm; this is an interchange approximation, not an exact symbolic conversion or a certified global error bound. The original in-memory BREP is unchanged. Enclosed cavities retain their shell orientation. Missing or unsupported geometry rejects the import instead of dropping faces or cavities.

An imported curved STEP surface remains a BREP, but that does not imply every modelling operation supports it. Recognized ellipsoids and ruled spline walls regain their corresponding native operation paths. Bounded rational spline surfaces support regular transverse intersections and the corresponding supported face trimming. Singular contacts, coincident patches and discontinuous domains still have more limited editing support. STEP assembly occurrence graphs are not yet supported. Mesh-only STEP export remains faceted; importing it creates planar BREP faces and cannot recover the original smooth surfaces.

## Modelling APIs

- **Primitives:** `BrepPrimitives` creates native solids.
- **Booleans:** `IjoOps.booleanSolid` operates on BREPs. `booleanPreferSolid` also returns a validated display mesh. Offset cylinder/sphere union, subtraction and intersection retain native spherical and cylindrical surfaces with exact quadric intersection curves. Intersected spherical faces are divided into bounded analytic patches, preserving existing circular cap/band trims and subsequent cuts. Mixed BREP/mesh operands require an explicit conversion decision from the host.
- **Extrude:** `BrepExtrude` constructs profile and planar-face extrusions. Face extrusion retains line, circular and elliptical boundaries, including holes. Oblique extrusion can use an exact affine shear. `tryMovePlanarFace` and `attachFromFace` support direct editing where applicable.
- **Hollow:** `IjoOps.hollowSolid` offsets supported analytic faces and connects opening rims. An empty opening list produces a sealed cavity. Collapsed or unsupported offsets throw.
- **Pipe:** `BrepPipe.run` creates cylindrical runs and toroidal elbows, including hollow tubes. `IjoOps.pipePreferSolid` returns the BREP and its display mesh. Invalid reversals and singular bends throw.
- **Loft:** `BrepLoft.sections` builds through ordered, parallel planar sections. Native prisms, cylinders and cones are retained where exact; other supported correspondences use ruled surfaces. `IjoOps.loftPreferSolid` returns the BREP and display mesh. Two-section lofts can include paired holes.
- **Slice:** `IjoOps.sliceBrep` returns BREP bodies. Full spheres and supported capped cylinders have direct analytic plane cuts; other supported solids use trimmed booleans. Both sides must validate before a split succeeds.
- **Mirror and scale:** `BrepXform` preserves curves, surfaces and provenance. General nonsingular affine scale retains exact transformed geometry, including ellipsoids, rather than replacing it with triangles.
- **Fillet and chamfer:** `IjoOps.filletSolid` and `chamferSolid` support planar corners, full circular cylinder/cone rims, and smooth extrusion rims containing lines, arcs, ellipses and cubic/rational spline profiles. Chamfers support equal setbacks, two distances, and distance plus angle. See [blend geometry and limits](doc/blends.md). These remain bounded blend constructors, not support for every arbitrary surface junction.
- **Threads:** `BrepThreads.apply` retains native helix curves and ruled faces for shafts and bores. It supports partial lengths, both ends, handedness, multiple starts, tapered threads and gradual runouts. Complete cylindrical walls use direct boundary construction; interrupted walls use BREP intersections. See [native threads](doc/brep_threads.md) for dimension conventions, usage and limits.
- **Body decomposition:** `IjoOps.splitBrep` keeps inward cavity shells attached to their enclosing exterior.

The `MeshData` overloads and `Mesh*` APIs remain available for mesh-only workflows. A method returning only `MeshData` cannot carry a BREP into the caller's next operation. Integrations should retain the solid and use its mesh as derived display data.

## Validation, measurement and selection

`BrepHeal.topologyIsClosed` checks structural closure. `IjoOps.tessellateClosed` requires that closure and a closed display tessellation; unsupported results throw. Display seam repair is still possible, but does not replace the stored BREP. Structural closure alone does not prove absence of self-intersections or correct geometric incidence.

`IjoMeasure.ofSolid` integrates over supported surface domains. Its `isExact` flag describes whether analytic domains were used; numerical quadrature still has finite error. If a trim domain is unavailable, measurement uses the complete tessellation and reports `isExact: false`, preserving the original solid. Some bounds are conservative envelopes, especially for general affine solids and implicit intersection curves.

Use `BrepFaceSelect` and the BREP edge selectors when operating on solids. `FaceSelect` and `EdgeSelect` serve mesh regions and crease selections.

See the [tool hardening analysis](../../debug/tool-hardening/README.md), [STEP follow-up](../../debug/step-and-surface-hardening/README.md), [curved-intersection work](../../debug/curved-intersection-hardening/README.md) and [general spline intersections](../../debug/general-spline-intersections/README.md) for root causes, regression evidence and remaining work.

## Tests

Run `dart test` from this package. PaliCAD's host integration tests live in the main application's `test` directory and run with `flutter test`. Neither requires GUI interaction.

## License

MIT. See [LICENSE](LICENSE). Third-party notices are in [NOTICE](NOTICE).
