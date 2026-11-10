|Program | Language | Library | t<sub>g2r</sub> (s) | M<sub>g2r</sub> (Mb) | t<sub>r2g</sub> (s) | M<sub>r2g</sub> (Mb) |
|:-------|:---------|:--------|--------------------:|---------------------:|--------------------:|---------------------:|
|[bedcov\_c1\_cgr.c](bedcov_c1_cgr.c)          |C         |[cgranges.h](../lib/cgranges.h)|  5.2|  138.4 | 10.7|  19.1 |
|[bedcov\_cr1\_klib.cr](bedcov_cr1_klib.cr)    |Crystal   |[klib.cr](../lib/klib.cr)      |  8.8|  319.6 | 14.8|  40.7 |
|[bedcov\_nim1\_klib.nim](bedcov_nim1_klib.nim)|Nim       |[klib.nim](../lib/klib.nim)    | 16.6|  248.4 | 26.0|  34.1 |
|[bedcov\_jl1\_klib.jl](bedcov_jl1_klib.jl)    |Julia     |[Klib.jl](../lib/Klib.jl)      | 25.9|  428.1 | 63.0| 257.0 |
|[bedcov\_go1.go](bedcov_go1.go)               |Go        |                               | 34.0|  318.9 | 21.8|  47.3 |
|[bedcov\_js1\_cgr.js](bedcov_js1_cgr.jl)      |Javascript|                               | 76.4| 2219.9 | 80.0| 316.8 |
|[bedcov\_lua1\_cgr.lua](bedcov_lua1_cgr.lua)  |LuaJIT    |                               |174.1| 2668.0 |218.9| 364.6 |

* Crystal, Nim and Julia use a long contiguous memory block for loaded
  intervals on each chromosome. Javascript and LuaJIT instead use a list of
  separate objects.

* Bounds check is completely disabled for Nim and partially disabled for
  Crystal. Bounds check might be partially disabled in Julia, too. Not sure.

* For "g2r", sorting takes significant time. C uses a radix sort. Other
  languages rely on their standard libraries.

## Goose, Rust and Mojo

These programs were measured on a different machine from the table above, a
Ryzen 9 3950X running Ubuntu 24.04, so compare them with the C row here and
not with the rows above. The compilers were clang 18.1.3, Goose
v2026.3-269-g74c641c5, rustc 1.97.1 and Mojo 1.2.0 (2016d0cd).

|Program | Language | Library | t<sub>g2r</sub> (s) | M<sub>g2r</sub> (Mb) | t<sub>r2g</sub> (s) | M<sub>r2g</sub> (Mb) |
|:-------|:---------|:--------|--------------------:|---------------------:|--------------------:|---------------------:|
|[bedcov\_goose2\_lapper.goose](bedcov_goose2_lapper.goose)|Goose|[goose-bio][goose-bio] lapper| 1.44| 178.8 | 2.46| 25.5 |
|[bedcov\_goose1\_iitree.goose](bedcov_goose1_iitree.goose)|Goose|[goose-bio][goose-bio] iitree| 1.95| 138.2 | 3.24| 20.2 |
|[bedcov\_rs1\_lapper.rs](bedcov_rs1_lapper.rs)            |Rust |[rust-lapper][rust-lapper]   | 2.48| 212.0 | 3.63| 30.0 |
|[bedcov\_c1\_cgr.c](bedcov_c1_cgr.c)                      |C    |[cgranges.h](../lib/cgranges.h)| 2.96| 139.2 | 5.66| 20.2 |
|[bedcov\_mojo1\_lapper.mojo](bedcov_mojo1_lapper.mojo)    |Mojo |[mojo-lapper][mojo-lapper]   | 6.68| 184.1 | 6.06| 33.5 |

* "t" is the median CPU time (user plus system) of seven runs pinned to one
  core, with the output sent to `/dev/null`. "M" is the peak resident set
  size reported by GNU time, also the median of seven runs. Each timed run
  uses a fresh copy of the executable. All five programs print the same bytes
  on both Biofast BED datasets.
  [run\_bench.py](../run_bench.py) builds, checks and times them:
  `./run_bench.py <data directory> --cc clang --cpu 6 --rounds 7`.

* Everything was built as the Makefiles build it, with `CC=clang`:
  `-O3 -march=native` for C and for the C that Goose generates, and
  `cargo install` with `RUSTFLAGS="-C target-cpu=native"` for Rust.
  Mojo uses `-O3 --march native` with the compiler and Mojoverse dependencies
  pinned by the root `pixi.toml` and `pixi.lock`.
  Goose's `-O2` selects its inlining level; `lib/goose_cc.sh` appends the
  native C flags after Goose's default C flags. Native builds are the default.

* The iitree program is a port of the C. goose-bio's iitree is cgranges'
  radix sort and implicit interval tree, with one difference: it calls a
  block with each overlap, where `cr_overlap` fills an array of indices that
  the caller reads the intervals back through.

* The Goose and Rust lapper programs use a different index. A lapper holds
  the sorted intervals of one contig with an index over blocks of them, so these
  programs keep one lapper per contig and look the contig up themselves. The
  index is the extra memory.

* The Mojo program uses [mojo-lapper][mojo-lapper] with one index per contig
  and [mojo-io][mojo-io] for buffered input/output. This version of mojo-lapper
  uses comparison sorting and a maximum-interval-length search bound, not
  the Goose/Rust block-max index. Its `find` appends overlaps to a reusable
  list. These are different implementations, not a controlled language-only
  comparison. Mojo BED contig names must be valid UTF-8.

* These are whole-program times, and much of the C's time is not spent in
  cgranges. The Goose and Rust programs parse numbers with a digit loop and
  format output into a 64 KiB buffer, where the C calls `atol` and `printf`.
  [Measured in goose-bio][why] with clang, doing the same in the C took it
  from 2.93 s to 2.38 s on g2r and from 5.59 s to 3.87 s on r2g.

* The Goose, Rust and Mojo programs read uncompressed files only and have
  no `-c` or `-C` option.

[goose-bio]: https://github.com/sstadick/goose-bio
[rust-lapper]: https://github.com/sstadick/rust-lapper
[mojo-lapper]: https://github.com/sstadick/mojo-lapper
[mojo-io]: https://github.com/sstadick/mojo-io
[why]: https://github.com/sstadick/goose-bio/blob/v0.2.0/bench/BIOFAST.md#why-goose-is-ahead-of-c-on-bedcov
