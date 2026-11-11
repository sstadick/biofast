// `bedcov_c1_cgr.c` in Goose with a different index: load one BED file into
// one goose-bio `lapper` per contig, then print, for every interval of a
// second BED file, how many loaded intervals overlap it and how many of its
// bases they cover.
//
// A lapper indexes the intervals of one sequence, so this program keeps the
// contig names itself, where cgranges and `iitree` keep them in the index.
// `../lib/bed.goose` reads the files. Input must not be compressed.
import std;
import os;
import dictionary;

import io from "github.com/sstadick/goose-io@v0.2.2";
import bstr from "github.com/sstadick/goose-bstr@v0.3.1";
import simd_bytes from "github.com/sstadick/goose-simd@v0.2.0";
import bio from "github.com/sstadick/goose-bio@v0.2.0";
import bio.lapper;
import bio.radix_sort;
import bed from "../lib";

// The contig names seen so far, numbered from 0 in order of first use.
root struct Contigs {
    names: u8[varint][>..]              = [],
    ids:   dictionary<const u8[:], i32> = dictionary<const u8[:], i32> {},
}

// The number of the contig called name, which is added if it is new.
fn id_of(contigs: Contigs&, name: const u8[:]) -> i32 {
    let known = contigs.ids.get(name);
    if known { return known; }
    let id = contigs.ids.count as! i32;
    contigs.ids.insert(contigs.names.push(name), id);
    id
}

checkfree fn contig_of(interval: lapper::Interval<i32, i32>) -> i32 { interval.val }

fn main() {
    let arguments = args();
    if arguments.len != 3 {
        print("Usage: bedcov <loaded.bed> <streamed.bed>");
        return;
    }

    // Every loaded interval carries its contig's number as its value. Sorting
    // by that number gives each contig one run of intervals to index.
    let contigs = Contigs { .. };
    var intervals: lapper::Interval<i32, i32>[>..] = [];
    let loaded = bed::each_bed_record(arguments[1]) { contig, start, end =>
        // cgranges drops intervals that end before they start.
        if start <= end { intervals.push(lapper::Interval { start: start, stop: end, val: contigs.id_of(contig) }); }
    };
    if !loaded { bed::fail(arguments[1]); }
    radix_sort::radix_sort(intervals, contig_of);

    var index: i32[>..] = [];
    var lappers: lapper::Lapper<i32, i32>[>..] = [];
    var first = 0;
    for id in 0..contigs.ids.count {
        var last = first;
        while last < intervals.len && intervals[last].val == id { last += 1; }
        lappers.push(lapper::new_lapper(intervals[first..last], index));
        first = last;
    }

    var out: u8[>..] = [];
    let streamed = bed::each_bed_record(arguments[2]) { contig, start1, end1 =>
        var count = 0;
        var coverage = 0;
        var covered_start = 0;
        var covered_end = 0;
        let id = contigs.ids.get(contig);
        if id {
            lappers[id].find(start1, end1) { interval =>
                // Overlaps arrive by start, so each one extends the covered
                // run or starts a new one.
                let start0 = max(interval.start, start1);
                let end0 = min(interval.stop, end1);
                if start0 > covered_end {
                    coverage += covered_end - covered_start;
                    covered_start = start0;
                    covered_end = end0;
                } else if end0 > covered_end {
                    covered_end = end0;
                }
                count += 1;
            };
        }
        coverage += covered_end - covered_start;

        format(out, contig, "\t", start1, "\t", end1, "\t", count, "\t", coverage, "\n");
        if out.len >= bed::BUFFER_SIZE {
            write_stdout(out);
            out.clear();
        }
    };
    write_stdout(out);
    if !streamed { bed::fail(arguments[2]); }
}
