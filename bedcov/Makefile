CFLAGS=-g -Wall -O3 -march=native
LIBS=-lz
PROG=bedcov_c1_cgr bedcov_nim1_klib bedcov_cr1_klib bedcov_go1 bedcov_d_dgr bedcov_fsharp1 \
	bedcov_goose1_iitree bedcov_goose2_lapper
NIM=nim
CRYSTAL=crystal
GO=go
LDC=ldc2
GOOSE=goose
# Goose compiles to C, which this script builds with -O3 -march=native.
GOOSE_CC=$(abspath ../lib/goose_cc.sh)
RUSTFLAGS=-C target-cpu=native
CARGO=cargo
PROJECT_ROOT := $(shell realpath ../)

.PHONY:all clean

all:$(PROG)

bedcov_c1_cgr:bedcov_c1_cgr.c ../lib/cgranges.h ../lib/cgranges.c
	$(CC) $(CFLAGS) -I../lib -o $@ $< ../lib/cgranges.c $(LIBS)

bedcov_nim1_klib:bedcov_nim1_klib.nim ../lib/klib.nim
	$(NIM) c -d:release -d:nimEmulateOverflowChecks --bound_checks:off -p:../lib -o:$@ $<

bedcov_cr1_klib:bedcov_cr1_klib.cr ../lib/klib.cr
	$(CRYSTAL) build --release $<

bedcov_go1:bedcov_go1.go
	$(GO) build $<

bedcov_d_dgr:bedcov_d_dgr.d
	$(LDC) -O3 -release -flto=full $<

bedcov_fsharp1:fsharp/IITree.fs fsharp/Bedcov.fs
	dotnet build -c Release fsharp/

bedcov_goose1_iitree:bedcov_goose1_iitree.goose ../lib/bed.goose ../lib/goose_cc.sh
	CC=$(CC) $(GOOSE) build $< --cc $(GOOSE_CC) -O2 -o $@

bedcov_goose2_lapper:bedcov_goose2_lapper.goose ../lib/bed.goose ../lib/goose_cc.sh
	CC=$(CC) $(GOOSE) build $< --cc $(GOOSE_CC) -O2 -o $@

# Installs bin/bedcov_lapper, built from bedcov_rs1_lapper.rs.
.PHONY: bedcov_rs
bedcov_rs:
	RUSTFLAGS="$(RUSTFLAGS)" $(CARGO) install --path $(PROJECT_ROOT) --force --bin bedcov_lapper --root .

clean:
	rm -fr *.dSYM $(PROG) bin
