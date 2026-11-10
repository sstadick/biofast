#!/bin/sh
# The C compiler for Goose programs. `goose build` runs the C compiler itself
# and has no option for its flags, so the Makefiles pass this script as the
# compiler to build Goose's C with the flags they build the C programs with.
exec "${CC:-cc}" "$@" -O3 -march=native
