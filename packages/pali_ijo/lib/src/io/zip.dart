import 'dart:convert';
import 'dart:typed_data';

import 'package:archive/archive.dart' as zip;

/// ZIP container used by 3MF. Stored output remains deterministic; input accepts
/// ordinary DEFLATE, data descriptors and ZIP64 through the portable decoder.
/// Entries are inflated and CRC-checked only when requested, so slicer G-code,
/// thumbnails and unrelated metadata do not consume uncompressed memory.
class ZipArchive {
  ZipArchive();

  final Map<String, Uint8List> _entries = {};
  final Map<String, zip.ArchiveFile> _encoded = {};

  Iterable<String> get names => {..._encoded.keys, ..._entries.keys};

  int get length => names.length;

  Uint8List? operator [](String name) {
    final cached = _entries[name];
    if (cached != null) return cached;
    final entry = _encoded[name];
    if (entry == null) return null;
    try {
      final content = entry.content as List<int>;
      if (content.length != entry.size ||
          zip.getCrc32(content) != entry.crc32) {
        throw const FormatException('Size or CRC checksum mismatch');
      }
      return _entries[name] = content is Uint8List
          ? content
          : Uint8List.fromList(content);
    } catch (error) {
      throw FormatException('Invalid ZIP entry "$name": $error');
    }
  }

  void addBytes(String name, List<int> bytes) {
    _entries[name] = Uint8List.fromList(bytes);
  }

  void addText(String name, String text) => addBytes(name, utf8.encode(text));

  String? readText(String name) {
    final bytes = this[name];
    return bytes == null ? null : utf8.decode(bytes);
  }

  /// The archive as ZIP bytes.
  Uint8List toBytes() {
    final out = BytesBuilder();
    final directory = BytesBuilder();
    var count = 0;

    for (final key in names) {
      final name = utf8.encode(key);
      final data = this[key]!;
      final crc = _crc32(data);
      final offset = out.length;

      final local = ByteData(30);
      local.setUint32(0, 0x04034b50, Endian.little); // local file header
      local.setUint16(4, 20, Endian.little); // version needed
      local.setUint16(6, 0, Endian.little); // flags
      local.setUint16(8, 0, Endian.little); // stored
      local.setUint16(10, 0, Endian.little); // time
      local.setUint16(12, 0x21, Endian.little); // date: 1980-01-01
      local.setUint32(14, crc, Endian.little);
      local.setUint32(18, data.length, Endian.little);
      local.setUint32(22, data.length, Endian.little);
      local.setUint16(26, name.length, Endian.little);
      local.setUint16(28, 0, Endian.little); // extra
      out.add(local.buffer.asUint8List());
      out.add(name);
      out.add(data);

      final central = ByteData(46);
      central.setUint32(0, 0x02014b50, Endian.little); // central directory
      central.setUint16(4, 20, Endian.little); // version made by
      central.setUint16(6, 20, Endian.little); // version needed
      central.setUint16(8, 0, Endian.little); // flags
      central.setUint16(10, 0, Endian.little); // stored
      central.setUint16(12, 0, Endian.little); // time
      central.setUint16(14, 0x21, Endian.little); // date
      central.setUint32(16, crc, Endian.little);
      central.setUint32(20, data.length, Endian.little);
      central.setUint32(24, data.length, Endian.little);
      central.setUint16(28, name.length, Endian.little);
      central.setUint16(30, 0, Endian.little); // extra
      central.setUint16(32, 0, Endian.little); // comment
      central.setUint16(34, 0, Endian.little); // disk
      central.setUint16(36, 0, Endian.little); // internal attrs
      central.setUint32(38, 0, Endian.little); // external attrs
      central.setUint32(42, offset, Endian.little);
      directory.add(central.buffer.asUint8List());
      directory.add(name);
      count++;
    }

    final directoryOffset = out.length;
    final directoryBytes = directory.toBytes();
    out.add(directoryBytes);

    final end = ByteData(22);
    end.setUint32(0, 0x06054b50, Endian.little); // end of central directory
    end.setUint16(4, 0, Endian.little); // this disk
    end.setUint16(6, 0, Endian.little); // disk with directory
    end.setUint16(8, count, Endian.little);
    end.setUint16(10, count, Endian.little);
    end.setUint32(12, directoryBytes.length, Endian.little);
    end.setUint32(16, directoryOffset, Endian.little);
    end.setUint16(20, 0, Endian.little); // comment
    out.add(end.buffer.asUint8List());

    return out.toBytes();
  }

  /// Reads a ZIP, walking the central directory rather than scanning for local
  /// headers, so an entry's own name and size come from the authoritative
  /// record.
  static ZipArchive fromBytes(Uint8List bytes) {
    if (_findEndOfCentralDirectory(bytes) < 0) {
      throw const FormatException(
        'Not a ZIP archive: no end-of-directory record',
      );
    }
    try {
      final decoder = zip.ZipDecoder();
      final decoded = decoder.decodeBytes(bytes);
      final headers = decoder.directory.fileHeaders;
      if (headers.length != decoder.directory.totalCentralDirectoryEntries ||
          decoder.directory.numberOfThisDisk != 0 ||
          decoder.directory.diskWithTheStartOfTheCentralDirectory != 0) {
        throw const FormatException('Incomplete or split ZIP directory');
      }
      // The central directory is authoritative, including for 64-bit data
      // descriptors. Older decoder versions read those trailing sizes as
      // 32-bit values, so do not use their local-header size for validation.
      final metadata = {for (final h in headers) h.filename: h};
      if (metadata.length != headers.length) {
        throw const FormatException('Duplicate ZIP entry');
      }
      final archive = ZipArchive();
      for (final file in decoded.files) {
        if (!file.isFile) continue;
        final header = metadata[file.name]!;
        file.size = header.uncompressedSize!;
        file.crc32 = header.crc32;
        if (archive._encoded.containsKey(file.name)) {
          throw FormatException('Duplicate ZIP entry "${file.name}"');
        }
        archive._encoded[file.name] = file;
      }
      return archive;
    } catch (error) {
      throw FormatException('Invalid ZIP archive: $error');
    }
  }

  /// The end-of-central-directory record sits at the end, but a trailing
  /// comment can push it back by up to 64 KiB, so it is searched for.
  static int _findEndOfCentralDirectory(Uint8List bytes) {
    if (bytes.length < 22) return -1;
    final earliest = bytes.length - 22 - 65535;
    final lowest = earliest < 0 ? 0 : earliest;
    for (var i = bytes.length - 22; i >= lowest; i--) {
      if (bytes[i] == 0x50 &&
          bytes[i + 1] == 0x4b &&
          bytes[i + 2] == 0x05 &&
          bytes[i + 3] == 0x06) {
        return i;
      }
    }
    return -1;
  }

  static Uint32List? _crcTable;

  static int _crc32(List<int> bytes) {
    final table = _crcTable ??= _buildCrcTable();
    var crc = 0xffffffff;
    for (final b in bytes) {
      crc = table[(crc ^ b) & 0xff] ^ (crc >> 8);
    }
    return (crc ^ 0xffffffff) & 0xffffffff;
  }

  static Uint32List _buildCrcTable() {
    final table = Uint32List(256);
    for (var i = 0; i < 256; i++) {
      var c = i;
      for (var k = 0; k < 8; k++) {
        c = (c & 1) != 0 ? 0xedb88320 ^ (c >> 1) : c >> 1;
      }
      table[i] = c;
    }
    return table;
  }
}
