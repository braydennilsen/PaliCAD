# PaliCAD
PLACEHOLDER DIRECTORY TO BE UPLOADED SOON

Direct-modeling CAD for **macOS / Windows** and **iPad**, built with **Flutter / Dart**. Modelling is **[pali ijo](packages/pali_ijo/)**, a standalone pure-Dart kernel with no Flutter dependency.

Sketch on construction planes, extrude and edit solids, and save projects as `.pali` files.

**Version:** 0.1.0 · **License:** [PolyForm Internal Use 1.0.0](LICENSE)

---

## Features

- Sketch tools: line, polyline, circle, rectangle (corner / center / 3-point), polygon, spline, fillet, chamfer, mirror, project, and more
- 3D tools: extrude, fillet, chamfer, hollow, threads, pipe, slice, mirror, scale, primitives, boolean union / subtract / intersect
- Selection modes: object, face, edge, plane
- Construction planes and face-based plane creation
- Import / export STL, STEP, 3MF, DXF (2D), .pali custom project file
- Undo / redo, autosave, keyboard shortcuts, light / dark themes

---

## Requirements

| Tool | Notes |
|------|--------|
| [Flutter](https://docs.flutter.dev/get-started/install) | Stable channel (3.35+ recommended) |
| Xcode | macOS builds |
| Visual Studio (Desktop C++) | Windows builds |
---

## Build
```bash
git clone https://github.com/braydennilsen/palicad.git
cd palicad
flutter pub get

# macOS
flutter run -d macos

# Windows
flutter run -d windows

# Linux
TBD

# Android
TBD
```

---

## License

PaliCAD is licensed under **[PolyForm Internal Use 1.0.0](LICENSE)**, © 2026 Brayden Nilsen. The pali ijo kernel is **[MIT](packages/pali_ijo/LICENSE)** licensed, © 2026 Brayden Nilsen. Other Flutter dependencies keep their own licenses (Settings → About → All licenses).[
