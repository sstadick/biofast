#!/usr/bin/env python3
"""Build the C, Goose, and Rust programs with the Makefiles and compare them.

The flags are the Makefiles': -O3 -march=native for C and for the C that Goose
generates, and target-cpu=native for Rust. Each program's output is first
checked against the C program's. Then every program is run several times with
its output discarded, and the median CPU time (user plus system) and peak
memory are printed as Markdown tables.
"""

import argparse
from dataclasses import dataclass
import hashlib
import os
from pathlib import Path
import shlex
import shutil
import statistics
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parent
DATA_URL = "https://github.com/lh3/biofast/releases/tag/biofast-data-v1"
TIME = "/usr/bin/time"


@dataclass(frozen=True)
class Program:
    source: str  # In the benchmark's directory.
    language: str
    library: str
    target: str  # The Makefile target that builds it.
    executable: str  # What the target builds, in the benchmark's directory.


@dataclass(frozen=True)
class Benchmark:
    directory: str
    # Column name and the data files of each command line.
    cases: dict[str, tuple[str, ...]]
    # The first program is the reference for the others' output.
    programs: tuple[Program, ...]


BENCHMARKS = (
    Benchmark(
        "bedcov",
        {"g2r": ("ex-rna.bed", "ex-anno.bed"), "r2g": ("ex-anno.bed", "ex-rna.bed")},
        (
            Program("bedcov_c1_cgr.c", "C", "cgranges", "bedcov_c1_cgr", "bedcov_c1_cgr"),
            Program(
                "bedcov_goose1_iitree.goose",
                "Goose",
                "goose-bio iitree",
                "bedcov_goose1_iitree",
                "bedcov_goose1_iitree",
            ),
            Program(
                "bedcov_goose2_lapper.goose",
                "Goose",
                "goose-bio lapper",
                "bedcov_goose2_lapper",
                "bedcov_goose2_lapper",
            ),
            Program("bedcov_rs1_lapper.rs", "Rust", "rust-lapper", "bedcov_rs", "bin/bedcov_lapper"),
        ),
    ),
    Benchmark(
        "fqcnt",
        {"plain": ("M_abscessus_HiSeq.fq",)},
        (
            Program("fqcnt_c1_kseq.c", "C", "kseq.h", "fqcnt_c1_kseq", "fqcnt_c1_kseq"),
            Program(
                "fqcnt_goose1_fastx.goose",
                "Goose",
                "goose-bio fastx",
                "fqcnt_goose1_fastx",
                "fqcnt_goose1_fastx",
            ),
            Program("fqcnt_rs1_rustbio.rs", "Rust", "rust-bio", "fqcnt_rs", "bin/fqcnt_rustbio"),
            Program("fqcnt_rs2_needletail.rs", "Rust", "needletail", "fqcnt_rs", "bin/fqcnt_needletail"),
        ),
    ),
)


@dataclass(frozen=True)
class Sample:
    cpu: float  # User plus system seconds.
    peak: float  # Peak resident set size in MB.


def run(command: list[str], cwd: Path) -> None:
    result = subprocess.run(command, cwd=cwd, text=True, capture_output=True, check=False)
    if result.returncode != 0:
        sys.exit(f"{shlex.join(command)} failed:\n{result.stdout}{result.stderr}")


def build(cc: str) -> None:
    """Build every program with the flags in the Makefiles."""
    for benchmark in BENCHMARKS:
        targets = list(dict.fromkeys(program.target for program in benchmark.programs))
        print(f"building {benchmark.directory}: {' '.join(targets)}", file=sys.stderr)
        run(["make", "-B", f"CC={cc}", *targets], ROOT / benchmark.directory)


def command(benchmark: Benchmark, program: Program, case: str, data: Path) -> list[str]:
    executable = ROOT / benchmark.directory / program.executable
    return [str(executable), *(str(data / file) for file in benchmark.cases[case])]


def output_digest(command: list[str]) -> str:
    """Run command and return a digest of what it wrote to standard output."""
    digest = hashlib.md5()
    with subprocess.Popen(command, stdout=subprocess.PIPE) as process:
        while chunk := process.stdout.read(1 << 20):
            digest.update(chunk)
    if process.returncode != 0:
        sys.exit(f"{shlex.join(command)} failed with status {process.returncode}")
    return digest.hexdigest()


def check(data: Path) -> None:
    """Exit unless every program prints what its benchmark's first program prints."""
    print("checking every program's output against C", file=sys.stderr)
    wrong = []
    for benchmark in BENCHMARKS:
        for case in benchmark.cases:
            digests = [
                output_digest(command(benchmark, program, case, data))
                for program in benchmark.programs
            ]
            wrong += [
                f"{program.source} {case}"
                for program, digest in zip(benchmark.programs, digests)
                if digest != digests[0]
            ]
    if wrong:
        sys.exit(f"output differs from C: {', '.join(wrong)}")


def measure(prefix: list[str], command: list[str], report: Path) -> Sample:
    """Run command with its output discarded and return what it cost."""
    # GNU time is a small parent. Linux carries a process's peak resident set
    # size across exec, so a child of this Python process would never report
    # less than Python's own footprint.
    command = [*prefix, TIME, "-f", "%U %S %M", "-o", str(report), *command]
    status = subprocess.run(command, stdout=subprocess.DEVNULL, check=False).returncode
    if status != 0:
        sys.exit(f"{shlex.join(command)} failed with status {status}")
    user, system, peak_kb = report.read_text().split()
    return Sample(float(user) + float(system), int(peak_kb) / 1024)


def table(header: list[str], rows: list[list[str]]) -> str:
    """Format a Markdown table whose columns also line up as plain text."""
    left = 3  # Program, language, and library read left to right.
    widths = [max(len(row[i]) for row in [header, *rows]) for i in range(len(header))]
    lines = [
        "|"
        + "|".join(
            cell.ljust(width) if i < left else cell.rjust(width)
            for i, (cell, width) in enumerate(zip(row, widths))
        )
        + "|"
        for row in [header, *rows]
    ]
    rule = "|" + "|".join(
        ":" + "-" * (width - 1) if i < left else "-" * (width - 1) + ":"
        for i, width in enumerate(widths)
    ) + "|"
    return "\n".join([lines[0], rule, *lines[1:]])


def report(benchmark: Benchmark, samples: dict[tuple[str, str], list[Sample]]) -> None:
    header = ["Program", "Language", "Library"]
    for case in benchmark.cases:
        header += [f"t<sub>{case}</sub> (s)", f"M<sub>{case}</sub> (Mb)"]
    rows = []
    for program in benchmark.programs:
        name = program.source.replace("_", "\\_")
        row = [f"[{name}]({program.source})", program.language, program.library]
        for case in benchmark.cases:
            cpu = statistics.median(sample.cpu for sample in samples[program.source, case])
            peak = statistics.median(sample.peak for sample in samples[program.source, case])
            row += [f"{cpu:.2f}", f"{peak:.1f}"]
        rows.append(row)
    print(table(header, rows))
    print()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "data", type=Path, help=f"directory of the unpacked biofast-data-v1 files, from {DATA_URL}"
    )
    parser.add_argument("--rounds", type=int, default=5, help="timed runs of each program (default: 5)")
    parser.add_argument("--cpu", type=int, help="pin each timed run to this CPU with taskset")
    parser.add_argument(
        "--cc",
        default=os.environ.get("CC", "cc"),
        help="C compiler for the C programs and for Goose's C (default: $CC, or cc)",
    )
    args = parser.parse_args()
    if args.rounds < 1:
        parser.error("--rounds must be positive")
    for tool in [TIME, args.cc, "make", "goose", "cargo", *(["taskset"] if args.cpu is not None else [])]:
        if shutil.which(tool) is None:
            parser.error(f"{tool} not found")
    files = sorted({file for benchmark in BENCHMARKS for case in benchmark.cases.values() for file in case})
    missing = [file for file in files if not (args.data / file).is_file()]
    if missing:
        parser.error(f"{args.data} lacks {', '.join(missing)}; unpack the data from {DATA_URL}")

    samples: dict[tuple[str, str], list[Sample]] = {}
    with tempfile.TemporaryDirectory(prefix="biofast-bench-") as tmp:
        scratch = Path(tmp)
        build(args.cc)
        check(args.data)
        prefix = [] if args.cpu is None else ["taskset", "-c", str(args.cpu)]
        # Interleave the programs so slow drift affects every program alike.
        for round_number in range(1, args.rounds + 1):
            print(f"round {round_number}/{args.rounds}", file=sys.stderr)
            for benchmark in BENCHMARKS:
                for case in benchmark.cases:
                    for program in benchmark.programs:
                        # Byte-identical executables showed different timings
                        # in local runs. Sample fresh copies for every program;
                        # copying itself is outside the timed region.
                        executable, *files = command(benchmark, program, case, args.data)
                        copy = scratch / f"{Path(executable).name}-{round_number}"
                        shutil.copy2(executable, copy)
                        sample = measure(prefix, [str(copy), *files], scratch / "report")
                        copy.unlink()
                        samples.setdefault((program.source, case), []).append(sample)

    for benchmark in BENCHMARKS:
        report(benchmark, samples)


if __name__ == "__main__":
    main()
