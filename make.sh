#!/bin/bash
#
#  This file is part of varr-numerics, an experimental variable resolution
#  primitive numerics library.
#
#  Copyright (C) 2019 Cardiff University
#
#  varr-numerics is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  varr-numerics is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with varr-numerics.  If not, see <http://www.gnu.org/licenses/>.
#

set -x

C_COMPILER=$CC

BUILD_OUTPUT_DIRECTORY=./build

mkdir -p $BUILD_OUTPUT_DIRECTORY

#
# gcc:
#

if [ ! -n "$VARR_CFLAGS" ]; then
   VARR_CFLAGS=$(pkg-config --cflags gsl)
fi

VARR_CFLAGS+=" -O3 -mtune=native -march=native -ffast-math -ffinite-math-only -mavx "

#
# icc 2018/update-2:
#

# VARR_CFLAGS+="-O3 -xCOMMON-AVX512 -mtune=native -fp-model fast=2 -qopt-report=5 -qopt-report-phase=vec"

if [ ! -n "$VARR_LDFLAGS" ]; then
   VARR_LDFLAGS=$(pkg-config --libs gsl)
fi

COMMON_VARR_CFLAGS="-D_POSIX_C_SOURCE=200112L -std=c99"
VARR_CFLAGS=$VARR_CFLAGS" "$COMMON_VARR_CFLAGS

rm -f $BUILD_OUTPUT_DIRECTORY/*.o
rm -f $BUILD_OUTPUT_DIRECTORY/libvarr.so

$C_COMPILER $VARR_CFLAGS -fPIC ./src/varr_sin.c -I./varr/ -c -o $BUILD_OUTPUT_DIRECTORY/varr_sin.o -Werror
$C_COMPILER $VARR_CFLAGS -fPIC ./src/varr_cos.c -I./varr/ -c -o $BUILD_OUTPUT_DIRECTORY/varr_cos.o -Werror
$C_COMPILER $VARR_CFLAGS -fPIC ./src/varr_phasor.c -I./varr/ -c -o $BUILD_OUTPUT_DIRECTORY/varr_phasor.o -Werror
$C_COMPILER $VARR_CFLAGS -fPIC ./src/varr_sixth_root.c -I./varr/ -c -o $BUILD_OUTPUT_DIRECTORY/varr_sixth_root.o -Werror
$C_COMPILER $VARR_CFLAGS -fPIC ./src/varr_atan.c -I./varr/ -c -o $BUILD_OUTPUT_DIRECTORY/varr_atan.o -Werror
$C_COMPILER $VARR_CFLAGS -fPIC ./src/varr_log.c -I./varr/ -c -o $BUILD_OUTPUT_DIRECTORY/varr_log.o -Werror
$C_COMPILER $VARR_CFLAGS -fPIC ./src/varr_exp.c -I./varr/ -c -o $BUILD_OUTPUT_DIRECTORY/varr_exp.o -Werror
$C_COMPILER $VARR_CFLAGS -fPIC ./src/varr_3_over_4.c -I./varr/ -c -o $BUILD_OUTPUT_DIRECTORY/varr_3_over_4.o -Werror
$C_COMPILER $VARR_CFLAGS -fPIC ./src/varr_extimer.c -I./varr/ -c -o $BUILD_OUTPUT_DIRECTORY/varr_extimer.o -Werror
$C_COMPILER $VARR_CFLAGS -fPIC ./src/varr_general_bound_linbuf.c -I./varr/ -c -o $BUILD_OUTPUT_DIRECTORY/varr_general_bound_linbuf.o -Werror

$C_COMPILER $COMMON_VARR_CFLAGS -O0 -g3 -fPIC ./test/entry.c -I./varr/ -c -o $BUILD_OUTPUT_DIRECTORY/entry.o -Werror
$C_COMPILER $COMMON_VARR_CFLAGS -O0 -g3 -fPIC ./test/varr_sixthroot.test.c -I./varr/ -c -o $BUILD_OUTPUT_DIRECTORY/varr_sixthroot.test.o -Werror
$C_COMPILER $COMMON_VARR_CFLAGS -O0 -g3 -fPIC ./test/varr_3_over_4.test.c -I./varr/ -c -o $BUILD_OUTPUT_DIRECTORY/varr_3_over_4.test.o -Werror
$C_COMPILER $COMMON_VARR_CFLAGS -O0 -g3 -fPIC ./test/varr_phasor.test.c -I./varr/ -c -o $BUILD_OUTPUT_DIRECTORY/varr_phasor.test.o -Werror
$C_COMPILER $COMMON_VARR_CFLAGS -O0 -g3 -fPIC ./test/varr_atan.test.c -I./varr/ -c -o $BUILD_OUTPUT_DIRECTORY/varr_atan.test.o -Werror
$C_COMPILER $COMMON_VARR_CFLAGS -O0 -g3 -fPIC ./test/varr_log.test.c -I./varr/ -c -o $BUILD_OUTPUT_DIRECTORY/varr_log.test.o -Werror
$C_COMPILER $COMMON_VARR_CFLAGS -O0 -g3 -fPIC ./test/varr_sin.test.c -I./varr/ -c -o $BUILD_OUTPUT_DIRECTORY/varr_sin.test.o -Werror
$C_COMPILER $COMMON_VARR_CFLAGS -O0 -g3 -fPIC ./test/varr_exp.test.c -I./varr/ -c -o $BUILD_OUTPUT_DIRECTORY/varr_exp.test.o -Werror
$C_COMPILER $COMMON_VARR_CFLAGS -O0 -g3 -fPIC ./test/varr_general_bound_linbuf.test.c -I./varr/ -c -o $BUILD_OUTPUT_DIRECTORY/varr_general_bound_linbuf.test.o -Werror
$C_COMPILER $COMMON_VARR_CFLAGS -O0 -g3 -fPIC ./test/varr_sequence_analysis.test.c -I./varr/ -c -o $BUILD_OUTPUT_DIRECTORY/varr_sequence_analysis.test.o -Werror
$C_COMPILER $COMMON_VARR_CFLAGS -O0 -g3 -fPIC ./test/test_results.c -I./varr/ -c -o $BUILD_OUTPUT_DIRECTORY/test_results.o -Werror

$C_COMPILER $COMMON_VARR_CFLAGS -O0 -g3 -fPIC ./test/timings_double.c -I./varr/ -c -o $BUILD_OUTPUT_DIRECTORY/timings_double.o -Werror
$C_COMPILER $COMMON_VARR_CFLAGS -O0 -g3 -fPIC ./test/timings_complex.c -I./varr/ -c -o $BUILD_OUTPUT_DIRECTORY/timings_complex.o -Werror
$C_COMPILER $COMMON_VARR_CFLAGS -O0 -g3 -fPIC ./test/sequence_generation.c -I./varr/ -c -o $BUILD_OUTPUT_DIRECTORY/sequence_generation.o -Werror
$C_COMPILER $COMMON_VARR_CFLAGS -O0 -g3 -fPIC ./test/sequence_comparison.c -I./varr/ -c -o $BUILD_OUTPUT_DIRECTORY/sequence_comparison.o -Werror

$C_COMPILER $BUILD_OUTPUT_DIRECTORY/varr_extimer.o $BUILD_OUTPUT_DIRECTORY/varr_exp.o $BUILD_OUTPUT_DIRECTORY/varr_log.o $BUILD_OUTPUT_DIRECTORY/varr_sin.o $BUILD_OUTPUT_DIRECTORY/varr_cos.o $BUILD_OUTPUT_DIRECTORY/varr_phasor.o $BUILD_OUTPUT_DIRECTORY/varr_sixth_root.o $BUILD_OUTPUT_DIRECTORY/varr_3_over_4.o $BUILD_OUTPUT_DIRECTORY/varr_atan.o $BUILD_OUTPUT_DIRECTORY/varr_general_bound_linbuf.o -Werror --shared -o $BUILD_OUTPUT_DIRECTORY/libvarr.so $VARR_LDFLAGS -lm

$C_COMPILER $COMMON_VARR_CFLAGS -O0 -g3 ./test/varr_utils.c -I./varr/ -c -o $BUILD_OUTPUT_DIRECTORY/varr_utils.o -Werror
$C_COMPILER -O0 -g3 $BUILD_OUTPUT_DIRECTORY/entry.o $BUILD_OUTPUT_DIRECTORY/varr_utils.o $BUILD_OUTPUT_DIRECTORY/varr_sixthroot.test.o $BUILD_OUTPUT_DIRECTORY/varr_3_over_4.test.o $BUILD_OUTPUT_DIRECTORY/varr_exp.test.o $BUILD_OUTPUT_DIRECTORY/varr_phasor.test.o $BUILD_OUTPUT_DIRECTORY/varr_atan.test.o $BUILD_OUTPUT_DIRECTORY/varr_sin.test.o $BUILD_OUTPUT_DIRECTORY/varr_log.test.o $BUILD_OUTPUT_DIRECTORY/test_results.o $BUILD_OUTPUT_DIRECTORY/timings_double.o $BUILD_OUTPUT_DIRECTORY/timings_complex.o $BUILD_OUTPUT_DIRECTORY/sequence_comparison.o $BUILD_OUTPUT_DIRECTORY/sequence_generation.o $BUILD_OUTPUT_DIRECTORY/varr_sequence_analysis.test.o $BUILD_OUTPUT_DIRECTORY/varr_general_bound_linbuf.test.o -L$BUILD_OUTPUT_DIRECTORY/ -lvarr -o $BUILD_OUTPUT_DIRECTORY/test $VARR_LDFLAGS -lrt -lm
