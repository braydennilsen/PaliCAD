# Changelog

## 0.1.3 — Unreleased

 - Expanded Fillet and Chamfer support for circular rims, hollow bodies, inside corners, successive blends, and mitred shapes.
- Improved Hollow for chamfered bodies and solids containing extruded text, splines, and elliptical curves.
- Improved threading of multiple holes, rounded parts, partial threads, and custom sizes.
- Fixed further modeling operations on supported scaled and mirrored geometry.
- Improved extrusion and Boolean operations involving text, projected curves, touching outlines, and complex SVG profiles.
- Improved STEP import, assembly placement, curve precision, and repair of malformed boundaries.
- Expanded 3MF support for compressed packages, component instances, units, and mirrored geometry.
- Improved tessellation of curved surfaces, threaded geometry, thin features, and faces containing complex holes.
- Added native cuts through toroidal surfaces

## 0.1.2

- Fixed incomplete arcs in extruded profiles that could cause unions against curved walls to fail.

- Fixed STEP import crashes caused by collapsed triangles on trimmed cylindrical surfaces.
- Fixed polygon extrusions being mistaken for curves, preserving straight sides across side counts, holes and combined sketch regions.
- Removed the fixed 10,000 mm fillet radius cap; limits now follow the selected geometry.
- Fixed hollowing pointed cones while preserving wall thickness and a closed solid.
- Fixed face extrusions retaining overlapping walls when joining another part of the same body, and improved triangulation of concave faces.
- Fixed leftover surface boundaries in unions of stretched bodies, including old outlines visible in X-ray.

## 0.1.1

Pali Ijo's version now matches PaliCAD's desktop release version.

- Improved Boolean operations and slicing at shared faces, edges and curved boundaries, including unions of mirrored bodies.
- Fixed several fillet, inverse fillet and chamfer failures, including subsequent cuts through modified surfaces.
- Improved hollowing of sliced and spherical bodies.
- Improved internal-thread performance and support for nominal-size holes, partial lengths, multiple starts and tapered threads.
- Expanded STEP curve and surface support, improved handling of face holes, and added lighter display settings for large assemblies.
- Improved triangulation of planar, spherical, lofted and spline surfaces. Fixed incomplete curved-edge outlines.
- Improved mesh Boolean closure, boundary repair and separation of connected parts.
- Improved handling of saved curved geometry and expanded regression and previous-release compatibility coverage.

## 0.1.0 — Initial release

Pali Ijo is a Dart BREP CAD kernel for creating and editing 3D objects. It stores their faces, edges and surfaces, generates triangle meshes for display, and provides import and export tools. It can be used independently of Flutter.

### Creating and editing models

- BREP primitives: boxes, cylinders, cones, truncated cones, spheres, tori, regular prisms and pyramids.
- Extrusion of planar profiles, profiles with holes and existing planar faces.
- Boolean union, subtraction and intersection.
- Edge fillets and chamfers, including two-distance and angle-based chamfers.
- Hollowing with selected openings or an enclosed cavity.
- Pipes and tubes along paths, with rounded bends.
- Lofts between parallel planar sections, including paired holes in two-section lofts.
- Internal and external threads with partial lengths, both-end threading, multiple starts, left- or right-hand winding, taper and gradual runouts.
- Plane cuts and separation of disconnected bodies.
- Translation, rotation, mirroring and uniform or per-axis scaling.
- Curved and spline geometry, including supported surface intersections and trimmed faces.

### Display and model inspection

- BREP tessellation with draft, standard and fine display quality through the Dart API.
- Face and edge identifiers connecting the display mesh to the model geometry.
- Face and edge selection by geometry, position, direction, size and modelling history.
- Volume, surface area, centre of mass, bounding-box and sphericity measurements.
- Closure checks and topology cleanup.

### Mesh tools

- Revolutions, lathes, profile sweeps, helical sweeps and skins between rings.
- Mesh booleans, offsets, Minkowski sums and differences, fillets, chamfers, hollowing, threads and plane cuts.
- Vertex welding, triangulation, coplanar-face merging, removal of degenerate faces, boundary repair and separation of connected components.
- Convex-hull construction and convexity checks.

### Files and integration

- JSON saving and loading of BREP geometry, topology and source associations, or mesh documents.
- STEP import and export for supported BREP geometry and multiple parts.
- STL, OBJ and 3MF mesh import and export.
- Dart library API and data packing for worker isolates.
- Native command-line tool with 30 commands, file input and output, shell pipelines and JSON reports.
- WebAssembly entry point with partial browser support.
- Operation logging for diagnostics.

### Availability and limits

Modelling support depends on the surfaces and intersections involved. Some complex blends, offsets and surface contacts remain unsupported. General profile sweeps and revolutions currently produce meshes. Browser execution and some CLI options also have limitations described in the [README](README.md#current-limits).

Pali Ijo is available under the [MIT license](LICENSE). Third-party attribution is documented in [THIRD_PARTY_NOTICES.txt](THIRD_PARTY_NOTICES.txt).
