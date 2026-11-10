"""Check the Mojo command-line programs, including cases outside biofast-data-v1."""

from pathlib import Path
import random
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


class MojoPrograms(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        for directory, targets in (
            ("bedcov", ["bedcov_c1_cgr", "bedcov_mojo1_lapper"]),
            ("fqcnt", ["fqcnt_c1_kseq", "fqcnt_mojo1_fastx"]),
        ):
            subprocess.run(["make", "-C", str(ROOT / directory), "CC=clang", *targets], check=True)

    def compare(self, directory, programs, inputs, expected, stderr=b""):
        with tempfile.TemporaryDirectory() as tmp:
            files = []
            for i, data in enumerate(inputs):
                path = Path(tmp) / str(i)
                path.write_bytes(data)
                files.append(str(path))
            for program in programs:
                with self.subTest(program=program):
                    result = subprocess.run(
                        [str(ROOT / directory / program), *files], capture_output=True, check=True
                    )
                    self.assertEqual(result.stdout, expected)
                    self.assertEqual(result.stderr, stderr)

    def bedcov(self, loaded, queries, expected):
        self.compare("bedcov", ["bedcov_c1_cgr", "bedcov_mojo1_lapper"],
                     [loaded, queries], expected)

    def test_bed_boundaries_and_contigs(self):
        loaded = (b"chr2\t4\t8\nchr1\t10\t20\nchr1\t5\t12\nchr1\t10\t20\n"
                  b"chr1\t30\t30\nchr1\t80\t40\nnot a BED record\n")
        queries = (b"chr1\t0\t5\nchr1\t5\t10\nchr1\t11\t18\nchr1\t20\t30\n"
                   b"chr1\t29\t31\nchr2\t +3\t7tail\textra\r\nmissing\t1\t99")
        expected = (b"chr1\t0\t5\t0\t0\nchr1\t5\t10\t1\t5\n"
                    b"chr1\t11\t18\t3\t7\nchr1\t20\t30\t0\t0\n"
                    b"chr1\t29\t31\t1\t0\nchr2\t3\t7\t1\t3\n"
                    b"missing\t1\t99\t0\t0\n")
        self.bedcov(loaded, queries, expected)
        self.bedcov(b"", queries, b"".join(
            b"\t".join(line.split(b"\t")[:3]) + b"\t0\t0\n"
            for line in expected.splitlines()
        ))
        self.bedcov(loaded, b"", b"")
        self.bedcov(
            b"chr1\t0\t8\nchr1\t2147483640\t2147483647\n",
            b"chr1\t2\t5\nchr1\t2147483645\t2147483647\n",
            b"chr1\t2\t5\t1\t3\nchr1\t2147483645\t2147483647\t1\t2\n",
        )
        # A UTF-8 contig name spanning multiple reader buffers is retained.
        name = "chr-é".encode() * 12000
        self.bedcov(name + b"\t1\t9\n", name + b"\t4\t11", name + b"\t4\t11\t1\t5\n")

    def test_bed_union_against_brute_force(self):
        rng = random.Random(71)
        intervals = [(rng.randrange(150), rng.randrange(150, 350)) for _ in range(30)]
        intervals += [(rng.randrange(350), rng.randrange(351, 400)) for _ in range(30)]
        rng.shuffle(intervals)
        queries = [(i, i + 13) for i in range(0, 420, 7)]
        loaded = b"".join(f"chr1\t{s}\t{e}\n".encode() for s, e in intervals)
        data = b"".join(f"chr1\t{s}\t{e}\n".encode() for s, e in queries)
        expected = b""
        for start, stop in queries:
            overlaps = [(s, e) for s, e in intervals if s < stop and e > start]
            covered = {base for s, e in overlaps for base in range(max(s, start), min(e, stop))}
            expected += f"chr1\t{start}\t{stop}\t{len(overlaps)}\t{len(covered)}\n".encode()
        self.bedcov(loaded, data, expected)

    def test_fastx(self):
        cases = [
            (b"", b"0\t0\t0\n"),
            (b">a comment\nAC\nGT\n\n>b\nN", b"2\t5\t0\n"),
            (b"@a\r\nAC\r\nG\r\n+comment\r\n!\r\n@@\r\n"
             b"@b\nTA\n+\n++", b"2\t5\t5\n"),
            (b"@empty\n\n+\n\n>x\nAGT\n", b"2\t3\t0\n"),
            (b">empty\n>next\nAC\n", b"2\t2\t0\n"),
            (b"@long " + b"h" * 20000 + b"\n" + b"A" * 65539 + b"\n+\n"
             + b"!" * 65539 + b"\n", b"1\t65539\t65539\n"),
        ]
        for data, expected in cases:
            with self.subTest(length=len(data)):
                self.compare("fqcnt", ["fqcnt_c1_kseq", "fqcnt_mojo1_fastx"], [data], expected)

    def test_truncated_and_overlong_quality(self):
        for quality in [b"", b"!", b"!!!\n", b"!\n"]:
            with self.subTest(quality=quality):
                self.compare(
                    "fqcnt", ["fqcnt_c1_kseq", "fqcnt_mojo1_fastx"],
                    [b"@good\nG\n+\n!\n@bad\nAC\n+\n" + quality], b"1\t1\t1\n",
                    stderr=b"ERROR: malformated FASTX\n",
                )


if __name__ == "__main__":
    unittest.main()
