CFLAGS=-g -Wall -O3 -march=native
LIBS=-lz
PROG=fqcnt_c1_kseq fqcnt_nim1_klib fqcnt_go1 fqcnt_rs fqcnt_goose1_fastx fqcnt_mojo1_fastx
NIM=nim
GO=go
CARGO=cargo
GOOSE=goose
# Goose compiles to C, which this script builds with -O3 -march=native.
GOOSE_CC=$(abspath ../lib/goose_cc.sh)
RUSTFLAGS=-C target-cpu=native
AMM=amm
MILL=mill
PROJECT_ROOT := $(shell realpath ../)
MOJO=pixi run --manifest-path $(PROJECT_ROOT)/pixi.toml --locked mojo
MOJOFLAGS=-O3 --march native

.PHONY:all clean

all:$(PROG)

fqcnt_c1_kseq:fqcnt_c1_kseq.c ../lib/kseq.h
	$(CC) $(CFLAGS) -I../lib -o $@ $< $(LIBS)

fqcnt_nim1_klib:fqcnt_nim1_klib.nim ../lib/klib.nim
	$(NIM) c -d:release -d:nimEmulateOverflowChecks --bound_checks:off -p:../lib -o:$@ $<

fqcnt_go1:fqcnt_go1.go
	$(GO) build $<

fqcnt_goose1_fastx:fqcnt_goose1_fastx.goose ../lib/goose_cc.sh
	CC=$(CC) $(GOOSE) build $< --cc $(GOOSE_CC) -O2 -o $@

fqcnt_mojo1_fastx:fqcnt_mojo1_fastx.mojo ../pixi.toml ../pixi.lock
	$(MOJO) build $(MOJOFLAGS) $(abspath $<) -o $(abspath $@)

.PHONY: fqcnt_rs
fqcnt_rs:
	RUSTFLAGS="$(RUSTFLAGS)" $(CARGO) install --path $(PROJECT_ROOT) --force --bin fqcnt_rustbio --bin fqcnt_needletail --root .

.PHONY: rust
rust: fqcnt_rs
	@echo "Rust binaries are in ./bin/"

.PHONY: fqcnt_fsharp
fqcnt_fsharp:fsharp/fqcnt_fsharp.fs
	dotnet build -c Release fsharp/

# Scala - amm
# Ammnoite is required to compile and run the scala script. Please obtain from ammonite.io:
# $ sudo sh -c '(echo "#!/usr/bin/env sh" \
#     && curl -L https://github.com/lihaoyi/Ammonite/releases/download/2.1.4/2.13-2.1.4) > /usr/local/bin/amm \
#     && chmod +x /usr/local/bin/amm' \
#     && amm
fqcnt_scala_amm_fgbio: fqcnt_scala_fgbio.sc
	echo "java -Xmx4g -jar \$$(which $(AMM)) --no-home-predef $< \$$@" > $@ 
	chmod 755 $@

# Scala - mill: http://www.lihaoyi.com/mill/
fqcnt_scala_jar_fgbio: scala/tools/src/com/github/biofast/FgBio.scala
	cd scala && $(MILL) tools.localJar

clean:
	rm -fr *.dSYM $(PROG)
