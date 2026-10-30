/// The WasmGC entry point: the kernel, callable from JavaScript.
///
/// Built by `tool/build_wasm.sh`. This is the third of the three ways to
/// consume pali ijo, alongside importing the Dart package and running the
/// native executable — and the one a browser-based designer would use.
///
/// What it exposes is deliberately the *same* surface as the command line
/// rather than a hand-written binding per operation. One function that takes
/// an argument list and a document and returns a document means there is one
/// API to document, one to test, and nothing to fall out of step when an
/// operation gains a flag.
///
/// This file lives outside `lib/` because it is the only part of the package
/// that imports `dart:js_interop`, which has no meaning on the VM. Everything
/// it calls is plain Dart.
library;

import 'dart:convert';
import 'dart:js_interop';
import 'dart:js_interop_unsafe';
import 'dart:typed_data';

import 'package:pali_ijo/pali_ijo.dart';
import 'package:pali_ijo/src/cli/cli.dart';

void main() {
  final api = JSObject()
    ..setProperty('run'.toJS, _run.toJS)
    ..setProperty('formats'.toJS, _formats.toJS);
  globalContext.setProperty('paliIjo'.toJS, api);
}

/// Runs one command.
///
/// From JavaScript:
///
/// ```js
/// const box = await paliIjo.run(['box', '--size', '40,40,20'], null);
/// const part = await paliIjo.run(
///   ['fillet', '--radius', '4', '--edges', 'convex'],
///   box.stdout,
/// );
/// const info = await paliIjo.run(['measure', '--json'], part.stdout);
/// console.log(JSON.parse(info.stdout).volume);
/// ```
///
/// Returns `{code, stdout, stderr, files}`. A non-zero `code` means the
/// command failed and `stderr` says why; nothing throws, matching how the
/// same command behaves in a shell.
///
/// `files` carries anything written with `--out`, since a browser has no
/// filesystem to write it to. Its values are `Uint8Array`, so a binary export
/// survives — `paliIjo.run(['export', '--out', 'part.3mf'], solid)` then
/// `result.files['part.3mf']` is bytes ready for a `Blob`.
JSPromise<JSObject> _run(JSArray<JSString> argv, JSString? stdin) {
  final args = [
    for (final a in argv.toDart) a.toDart,
  ];
  final host = MemoryHost(
    stdin: stdin == null ? null : utf8.encode(stdin.toDart),
  );

  return (() async {
    final code = await IjoCli(host).run(args);
    final files = JSObject();
    for (final entry in host.files.entries) {
      files.setProperty(entry.key.toJS, entry.value.toJS);
    }
    return JSObject()
      ..setProperty('code'.toJS, code.toJS)
      ..setProperty('stdout'.toJS, host.stdout.toJS)
      ..setProperty('stderr'.toJS, host.stderr.toJS)
      ..setProperty('files'.toJS, files);
  })()
      .toJS;
}

/// The file extensions `--out` and `export --format` understand.
JSArray<JSString> _formats() => [
      for (final f in IjoFormat.values) f.extension.toJS,
    ].toJS;

extension on Uint8List {
  JSUint8Array get toJS => JSUint8Array.withLength(length)
    ..toDart.setAll(0, this);
}
