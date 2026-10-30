/// The `pali_ijo` executable.
///
/// Build it with `dart compile exe bin/pali_ijo.dart -o pali_ijo`, which is
/// what makes the kernel reachable from any language: a Python or Rust or C
/// program that can run a process and pipe JSON can drive it, with no FFI and
/// no bindings to maintain.
///
/// Everything except the `dart:io` wiring lives in `lib/src/cli/`, so the
/// commands can be tested in-process and can compile for a WasmGC target that
/// has no `dart:io` at all.
library;

import 'dart:io';
import 'dart:typed_data';

import 'package:pali_ijo/src/cli/cli.dart';

Future<void> main(List<String> args) async {
  final host = ProcessHost();
  final code = await IjoCli(host).run(args);
  await host.flush();
  exitCode = code;
}

/// An [IjoHost] backed by the real filesystem and the real process streams.
class ProcessHost implements IjoHost {
  @override
  Future<Uint8List> readFile(String path) async {
    final file = File(path);
    if (!await file.exists()) {
      throw CliFileException('No such file: $path');
    }
    try {
      return await file.readAsBytes();
    } on FileSystemException catch (e) {
      throw CliFileException('Could not read $path: ${e.message}');
    }
  }

  @override
  Future<void> writeFile(String path, List<int> bytes) async {
    try {
      final file = File(path);
      // Writing into a directory that does not exist yet is a normal thing to
      // ask for in a build script, so it is created rather than refused.
      await file.parent.create(recursive: true);
      await file.writeAsBytes(bytes);
    } on FileSystemException catch (e) {
      throw CliFileException('Could not write $path: ${e.message}');
    }
  }

  @override
  Future<Uint8List> readStdin() async {
    final chunks = <List<int>>[];
    await for (final chunk in stdin) {
      chunks.add(chunk);
    }
    final builder = BytesBuilder();
    for (final chunk in chunks) {
      builder.add(chunk);
    }
    return builder.toBytes();
  }

  @override
  void writeStdout(String text) => stdout.write('$text\n');

  @override
  void writeStdoutBytes(List<int> bytes) => stdout.add(bytes);

  @override
  void writeStderr(String text) => stderr.writeln(text);

  /// Waits for the output streams to drain.
  ///
  /// Without this, a large model written to stdout can be truncated when the
  /// process exits before the pipe has been emptied — which shows up as a
  /// corrupt file only for big parts, and only sometimes.
  Future<void> flush() async {
    try {
      await stdout.flush();
      await stderr.flush();
    } on Object {
      // A closed pipe — the far side of `| head` having exited, say — is not
      // an error worth reporting on the way out.
    }
  }
}
