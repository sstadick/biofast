"""Uncompressed BED coverage using Mojoverse's mojo-lapper and mojo-io.

Like the Goose/Rust programs, keep one index per contig and buffer output.
mojo-lapper uses a maximum-interval-length bound, not their block-max index.
"""

from std.sys import argv
from lapper.lapper import Interval, Lapper
from mojo_io import BufferedReader, BufferedWriter


def atoi(text: Span[UInt8, _]) -> Int32:
    var i = 0
    while i < len(text) and (text[i] == 32 or (text[i] >= 9 and text[i] <= 13)):
        i += 1
    var negative = False
    if i < len(text) and (text[i] == 45 or text[i] == 43):
        negative = text[i] == 45
        i += 1
    var value = Int64(0)
    while i < len(text) and text[i] >= 48 and text[i] <= 57:
        value = value * 10 + Int64(text[i] - 48)
        i += 1
    return Int32(-value if negative else value)


def coordinate(value: Int32) -> UInt32:
    # Preserve signed BED coordinate order in mojo-lapper's unsigned index.
    return UInt32(value) ^ 0x80000000


@fieldwise_init
struct Record(Copyable):
    var contig: String
    var start: Int32
    var stop: Int32


def parse_bed(line: Span[UInt8, _]) raises -> Optional[Record]:
    var first = -1
    var second = -1
    var third = len(line)
    for i in range(len(line)):
        if line[i] == 9:
            if first < 0:
                first = i
            elif second < 0:
                second = i
            else:
                third = i
                break
    if second < 0:
        return None
    return Record(
        String(from_utf8=line[:first]),
        atoi(line[first + 1 : second]),
        atoi(line[second + 1 : third]),
    )


def main() raises:
    var args = argv()
    if len(args) != 3:
        print("Usage: bedcov <loaded.bed> <streamed.bed>")
        return

    var ids = Dict[String, Int]()
    var intervals = List[List[Interval]]()
    var loaded = BufferedReader(open(args[1], "r"), buffer_capacity=65536)
    while True:
        var line = loaded.read_line()
        if not line:
            break
        var parsed = parse_bed(line.value())
        if not parsed:
            continue
        var record = parsed.take()
        if record.start > record.stop:
            continue
        if record.contig not in ids:
            ids[record.contig] = len(intervals)
            intervals.append(List[Interval]())
        intervals[ids[record.contig]].append(
            Interval(coordinate(record.start), coordinate(record.stop), 0)
        )
    loaded.close()

    var lappers = List[Lapper[]]()
    while len(intervals):
        var index = Lapper(intervals.pop(0))
        lappers.append(index^)

    var results = List[Interval]()
    var queries = BufferedReader(open(args[2], "r"), buffer_capacity=65536)
    var output = BufferedWriter(open("/dev/stdout", "w"), buffer_capacity=65536)
    while True:
        var line = queries.read_line()
        if not line:
            break
        var parsed = parse_bed(line.value())
        if not parsed:
            continue
        var record = parsed.take()
        results.clear()
        if record.contig in ids:
            lappers[ids[record.contig]].find(
                coordinate(record.start), coordinate(record.stop), results
            )
        var coverage = 0
        var covered_start = 0
        var covered_end = 0
        for interval in results:
            var start = max(
                Int(Int32(interval.start ^ 0x80000000)), Int(record.start)
            )
            var stop = min(
                Int(Int32(interval.stop ^ 0x80000000)), Int(record.stop)
            )
            if start > covered_end:
                coverage += covered_end - covered_start
                covered_start = start
                covered_end = stop
            elif stop > covered_end:
                covered_end = stop
        coverage += covered_end - covered_start
        output.write_string(
            String(
                record.contig,
                "\t",
                record.start,
                "\t",
                record.stop,
                "\t",
                len(results),
                "\t",
                coverage,
                "\n",
            )
        )
    queries.close()
    output.close()
