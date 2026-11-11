"""Uncompressed BED coverage using Mojoverse's mojo-lapper and mojo-io.

Like the Goose/Rust programs, keep one index per contig and buffer output.
"""

from std.sys import argv
from lapper.lapper import Interval, IntervalVisitor, Lapper
from mojo_io import BufferedReader, BufferedWriter


@fieldwise_init
struct Coverage(IntervalVisitor):
    """Accumulate ordered overlaps without storing and rereading each match."""

    var query_start: Int32
    var query_stop: Int32
    var count: Int
    var coverage: Int
    var covered_start: Int
    var covered_end: Int

    @always_inline
    def visit(mut self, interval: Interval):
        var start = max(Int(Int32(interval.start ^ 0x80000000)), Int(self.query_start))
        var stop = min(Int(Int32(interval.stop ^ 0x80000000)), Int(self.query_stop))
        if start > self.covered_end:
            self.coverage += self.covered_end - self.covered_start
            self.covered_start = start
            self.covered_end = stop
        elif stop > self.covered_end:
            self.covered_end = stop
        self.count += 1


@always_inline
def atoi(text: Span[UInt8, _], mut offset: Int) -> Int32:
    var i = offset
    # Tabs delimit fields, including empty fields; never skip into the next one.
    while i < len(text) and (
        text[i] == 32 or (text[i] >= 10 and text[i] <= 13)
    ):
        i += 1
    var negative = False
    if i < len(text) and (text[i] == 45 or text[i] == 43):
        negative = text[i] == 45
        i += 1
    var value = Int64(0)
    while i < len(text) and text[i] >= 48 and text[i] <= 57:
        value = value * 10 + Int64(text[i] - 48)
        i += 1
    offset = i
    return Int32(-value if negative else value)


def coordinate(value: Int32) -> UInt32:
    # Preserve signed BED coordinate order in mojo-lapper's unsigned index.
    return UInt32(value) ^ 0x80000000


@fieldwise_init
struct Record(Copyable):
    var contig_end: Int
    var start: Int32
    var stop: Int32


@always_inline
def parse_bed(line: Span[UInt8, _]) -> Optional[Record]:
    var first = 0
    while first < len(line) and line[first] != 9:
        first += 1
    var offset = first + 1
    var start = atoi(line, offset)
    # Usually atoi stops at the tab. Preserve atol's acceptance of suffixes
    # without rescanning every start coordinate to find its delimiter first.
    while offset < len(line) and line[offset] != 9:
        offset += 1
    if offset >= len(line):
        return None
    offset += 1
    var stop = atoi(line, offset)
    return Record(first, start, stop)


struct Contigs:
    """Own contig keys and cache the last lookup, without retaining reader views.
    """

    var ids: Dict[String, Int]
    var name: String
    var id: Int

    def __init__(out self):
        self.ids = Dict[String, Int]()
        self.name = String()
        self.id = -2  # Distinguish an uninitialized cache from a missing key.

    def resolve(
        mut self, name: Span[UInt8, _], create: Bool = False
    ) raises -> Int:
        if self.id == -2 or name != self.name.as_bytes():
            # Validate before retaining a new name. Equal bytes were already
            # validated, so runs on the same contig need neither a copy nor hash.
            self.name = String(from_utf8=name)
            self.id = self.ids.get(self.name, -1)
        if create and self.id < 0:
            self.id = len(self.ids)
            self.ids[self.name] = self.id
        return self.id


def main() raises:
    var args = argv()
    if len(args) != 3:
        print("Usage: bedcov <loaded.bed> <streamed.bed>")
        return

    var contigs = Contigs()
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
        var id = contigs.resolve(
            line.value()[: record.contig_end], record.start <= record.stop
        )
        if record.start > record.stop:
            continue
        if id == len(intervals):
            intervals.append(List[Interval]())
        intervals[id].append(
            Interval(coordinate(record.start), coordinate(record.stop), 0)
        )
    loaded.close()

    var lappers = List[Lapper[]]()
    while len(intervals):
        var index = Lapper(intervals.pop(0))
        lappers.append(index^)

    var text = String(capacity_bytes=65536)
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
        var accumulator = Coverage(record.start, record.stop, 0, 0, 0, 0)
        var id = contigs.resolve(line.value()[: record.contig_end])
        if id >= 0:
            lappers[id].visit(
                coordinate(record.start), coordinate(record.stop), accumulator
            )
        var coverage = accumulator.coverage + accumulator.covered_end - accumulator.covered_start
        # Write each number once into a reusable batch, rather than sizing and
        # allocating a new variadic String for every record.
        contigs.name.write_to(text)
        text.write_string("\t")
        record.start.write_to(text)
        text.write_string("\t")
        record.stop.write_to(text)
        text.write_string("\t")
        accumulator.count.write_to(text)
        text.write_string("\t")
        coverage.write_to(text)
        text.write_string("\n")
        if text.byte_length() >= 65536:
            output.write_string(text)
            text.resize(0)
    output.write_string(text)
    queries.close()
    output.close()
