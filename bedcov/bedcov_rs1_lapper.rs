//! `bedcov_c1_cgr.c` in Rust with a different index: load one BED file into
//! one [rust-lapper](https://github.com/sstadick/rust-lapper) `Lapper` per
//! contig, then print, for every interval of a second BED file, how many
//! loaded intervals overlap it and how many of its bases they cover.
//!
//! Input must not be compressed.

use rust_lapper::{Interval, Lapper};
use std::collections::HashMap;
use std::env;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::process;

const BUFFER_SIZE: usize = 0x10000;

/// C's `atol` narrowed to 32 bits: optional blanks and sign, then digits.
fn atoi(text: &[u8]) -> i32 {
    let mut rest = text;
    while let [b' ' | b'\t'..=b'\r', tail @ ..] = rest {
        rest = tail;
    }
    let negative = match rest {
        [sign @ (b'+' | b'-'), tail @ ..] => {
            rest = tail;
            *sign == b'-'
        }
        _ => false,
    };
    let mut value = 0i64;
    for &byte in rest.iter().take_while(|byte| byte.is_ascii_digit()) {
        value = value.wrapping_mul(10).wrapping_add((byte - b'0') as i64);
    }
    (if negative { value.wrapping_neg() } else { value }) as i32
}

/// `parse_bed`: the contig, start, and end of a line with at least three
/// tab-separated fields.
fn parse_bed(line: &[u8]) -> Option<(&[u8], i32, i32)> {
    let mut fields = line.split(|&byte| byte == b'\t');
    let contig = fields.next()?;
    let start = atoi(fields.next()?);
    let end = atoi(fields.next()?);
    Some((contig, start, end))
}

/// Call `visit` with each line of the file, without its line terminator.
fn for_each_line(path: &str, mut visit: impl FnMut(&[u8]) -> io::Result<()>) -> io::Result<()> {
    let mut reader = BufReader::with_capacity(BUFFER_SIZE, File::open(path)?);
    let mut line = Vec::new();
    loop {
        line.clear();
        if reader.read_until(b'\n', &mut line)? == 0 {
            return Ok(());
        }
        if line.last() == Some(&b'\n') {
            line.pop();
        }
        if line.len() > 1 && line.last() == Some(&b'\r') {
            line.pop();
        }
        visit(&line)?;
    }
}

/// `read_bed`: the intervals of the file, indexed by contig. Each carries its
/// line's number among the file's records, as in C.
fn read_bed(path: &str) -> io::Result<HashMap<Vec<u8>, Lapper<i32, i32>>> {
    let mut intervals: HashMap<Vec<u8>, Vec<Interval<i32, i32>>> = HashMap::new();
    let mut label = 0;
    for_each_line(path, |line| {
        if let Some((contig, start, stop)) = parse_bed(line) {
            // cgranges drops intervals that end before they start.
            if start <= stop {
                let interval = Interval { start, stop, val: label };
                match intervals.get_mut(contig) {
                    Some(known) => known.push(interval),
                    None => {
                        intervals.insert(contig.to_vec(), vec![interval]);
                    }
                }
            }
            label += 1;
        }
        Ok(())
    })?;
    Ok(intervals.into_iter().map(|(contig, intervals)| (contig, Lapper::new(intervals))).collect())
}

fn run(loaded: &str, streamed: &str) -> io::Result<()> {
    let lappers = read_bed(loaded)?;

    let mut out = BufWriter::with_capacity(BUFFER_SIZE, io::stdout().lock());
    for_each_line(streamed, |line| {
        let (contig, start1, end1) = match parse_bed(line) {
            Some(record) => record,
            None => return Ok(()),
        };
        let (mut count, mut coverage, mut covered_start, mut covered_end) = (0i64, 0i64, 0i64, 0i64);
        if let Some(lapper) = lappers.get(contig) {
            // Overlaps arrive by start, so each one extends the covered run
            // or starts a new one.
            for interval in lapper.find(start1, end1) {
                let start0 = interval.start.max(start1) as i64;
                let end0 = interval.stop.min(end1) as i64;
                if start0 > covered_end {
                    coverage += covered_end - covered_start;
                    covered_start = start0;
                    covered_end = end0;
                } else {
                    covered_end = covered_end.max(end0);
                }
                count += 1;
            }
        }
        coverage += covered_end - covered_start;
        out.write_all(contig)?;
        writeln!(out, "\t{}\t{}\t{}\t{}", start1, end1, count, coverage)
    })?;
    out.flush()
}

fn main() {
    let arguments: Vec<String> = env::args().collect();
    if arguments.len() != 3 {
        println!("Usage: bedcov <loaded.bed> <streamed.bed>");
        return;
    }
    if let Err(error) = run(&arguments[1], &arguments[2]) {
        eprintln!("bedcov: {}", error);
        process::exit(1);
    }
}
