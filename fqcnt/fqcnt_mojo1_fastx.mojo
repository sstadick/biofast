"""Count uncompressed FASTA/multiline FASTQ using Mojoverse's mojo-io.

Keep reusable sequence and quality buffers, like kseq and Goose fastx;
this is a parser benchmark, not just a line-length counter.
"""

from std.sys import argv
from mojo_io import BufferedReader


def main() raises:
    var args = argv()
    if len(args) < 2:
        print("Usage: fqcnt <in.fq>")
        return
    var reader = BufferedReader(open(args[1], "r"), buffer_capacity=16384)
    var name = List[UInt8]()
    var sequence = List[UInt8]()
    var quality = List[UInt8]()
    var pending_header = False
    var count = 0
    var bases = 0
    var qualities = 0
    var malformed = False
    while True:
        if not pending_header:
            name.clear()
            while True:
                var line = reader.read_line()
                if not line:
                    break
                var bytes = line.value()
                if len(bytes) and (bytes[0] == 62 or bytes[0] == 64):
                    name.extend(bytes[1:])
                    pending_header = True
                    break
            if not pending_header:
                break
        pending_header = False
        sequence.clear()
        quality.clear()
        var fastq = False
        while True:
            var line = reader.read_line()
            if not line:
                break
            var bytes = line.value()
            if len(bytes) == 0:
                continue
            if bytes[0] == 62 or bytes[0] == 64:
                name.clear()
                name.extend(bytes[1:])
                pending_header = True
                break
            if bytes[0] == 43:
                fastq = True
                break
            sequence.extend(bytes)
        if fastq:
            # Even an empty sequence consumes a quality line.
            while True:
                var line = reader.read_line()
                if not line:
                    malformed = len(quality) != len(sequence)
                    break
                quality.extend(line.value())
                if len(quality) >= len(sequence):
                    malformed = len(quality) != len(sequence)
                    break
            if malformed:
                break
        count += 1
        bases += len(sequence)
        qualities += len(quality)
    reader.close()
    print(String(count, "\t", bases, "\t", qualities))
    if malformed:
        var errors = open("/dev/stderr", "w")
        errors.write("ERROR: malformated FASTX\n")
        errors.close()
