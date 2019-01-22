# varr-numerics
Experimental simple variable-resolution primitive numerics.  This repository should remain private.

# Compilation

This project includes a basic build script named `make.sh`, located in its base directory. 

To build, run `make.sh` using `bash` in its own directory:

```shell
$ ./make.sh

+ C_COMPILER=gcc
+ BUILD_OUTPUT_DIRECTORY=./build
+ mkdir -p ./build
++ pkg-config --cflags gsl
+ CFLAGS='-O3 -mtune=native -march=native -ffast-math -ffinite-math-only -mavx '
+ COMMON_CFLAGS='-D_POSIX_C_SOURCE=200112L -std=c99'
+ CFLAGS='-O3 -mtune=native -march=native -ffast-math -ffinite-math-only -mavx  -D_POSIX_C_SOURCE=200112L -std=c99'
+ rm -f ./build/entry.o ./build/sequence_comparison.o ./build/sequence_generation.o ./build/test_results.o ./build/timings_complex.o ./build/timings_double.o ./build/varr_3_over_4.o ./build/varr_3_over_4.test.o ./build/varr_atan.o ./build/varr_atan.test.o ./build/varr_cos.o ./build/varr_exp.o ./build/varr_exp.test.o ./build/varr_extimer.o ./build/varr_general_bound_linbuf.o ./build/varr_general_bound_linbuf.test.o ./build/varr_log.o ./build/varr_log.test.o ./build/varr_phasor.o ./build/varr_phasor.test.o ./build/varr_sequence_analysis.test.o ./build/varr_sin.o ./build/varr_sin.test.o ./build/varr_sixth_root.o ./build/varr_sixthroot.test.o ./build/varr_utils.o
+ rm -f ./build/libvarr.so
+ gcc -O3 -mtune=native -march=native -ffast-math -ffinite-math-only -mavx -D_POSIX_C_SOURCE=200112L -std=c99 -fPIC ./src/varr_sin.c -I./include/ -c -o ./build/varr_sin.o -Werror
+ gcc -O3 -mtune=native -march=native -ffast-math -ffinite-math-only -mavx -D_POSIX_C_SOURCE=200112L -std=c99 -fPIC ./src/varr_cos.c -I./include/ -c -o ./build/varr_cos.o -Werror
+ gcc -O3 -mtune=native -march=native -ffast-math -ffinite-math-only -mavx -D_POSIX_C_SOURCE=200112L -std=c99 -fPIC ./src/varr_phasor.c -I./include/ -c -o ./build/varr_phasor.o -Werror
+ gcc -O3 -mtune=native -march=native -ffast-math -ffinite-math-only -mavx -D_POSIX_C_SOURCE=200112L -std=c99 -fPIC ./src/varr_sixth_root.c -I./include/ -c -o ./build/varr_sixth_root.o -Werror
+ gcc -O3 -mtune=native -march=native -ffast-math -ffinite-math-only -mavx -D_POSIX_C_SOURCE=200112L -std=c99 -fPIC ./src/varr_atan.c -I./include/ -c -o ./build/varr_atan.o -Werror
+ gcc -O3 -mtune=native -march=native -ffast-math -ffinite-math-only -mavx -D_POSIX_C_SOURCE=200112L -std=c99 -fPIC ./src/varr_log.c -I./include/ -c -o ./build/varr_log.o -Werror
+ gcc -O3 -mtune=native -march=native -ffast-math -ffinite-math-only -mavx -D_POSIX_C_SOURCE=200112L -std=c99 -fPIC ./src/varr_exp.c -I./include/ -c -o ./build/varr_exp.o -Werror
+ gcc -O3 -mtune=native -march=native -ffast-math -ffinite-math-only -mavx -D_POSIX_C_SOURCE=200112L -std=c99 -fPIC ./src/varr_3_over_4.c -I./include/ -c -o ./build/varr_3_over_4.o -Werror
+ gcc -O3 -mtune=native -march=native -ffast-math -ffinite-math-only -mavx -D_POSIX_C_SOURCE=200112L -std=c99 -fPIC ./src/varr_extimer.c -I./include/ -c -o ./build/varr_extimer.o -Werror
+ gcc -O3 -mtune=native -march=native -ffast-math -ffinite-math-only -mavx -D_POSIX_C_SOURCE=200112L -std=c99 -fPIC ./src/varr_general_bound_linbuf.c -I./include/ -c -o ./build/varr_general_bound_linbuf.o -Werror
+ gcc -D_POSIX_C_SOURCE=200112L -std=c99 -O0 -g3 -fPIC ./test/entry.c -I./include/ -c -o ./build/entry.o -Werror
+ gcc -D_POSIX_C_SOURCE=200112L -std=c99 -O0 -g3 -fPIC ./test/varr_sixthroot.test.c -I./include/ -c -o ./build/varr_sixthroot.test.o -Werror
+ gcc -D_POSIX_C_SOURCE=200112L -std=c99 -O0 -g3 -fPIC ./test/varr_3_over_4.test.c -I./include/ -c -o ./build/varr_3_over_4.test.o -Werror
+ gcc -D_POSIX_C_SOURCE=200112L -std=c99 -O0 -g3 -fPIC ./test/varr_phasor.test.c -I./include/ -c -o ./build/varr_phasor.test.o -Werror
+ gcc -D_POSIX_C_SOURCE=200112L -std=c99 -O0 -g3 -fPIC ./test/varr_atan.test.c -I./include/ -c -o ./build/varr_atan.test.o -Werror
+ gcc -D_POSIX_C_SOURCE=200112L -std=c99 -O0 -g3 -fPIC ./test/varr_log.test.c -I./include/ -c -o ./build/varr_log.test.o -Werror
+ gcc -D_POSIX_C_SOURCE=200112L -std=c99 -O0 -g3 -fPIC ./test/varr_sin.test.c -I./include/ -c -o ./build/varr_sin.test.o -Werror
+ gcc -D_POSIX_C_SOURCE=200112L -std=c99 -O0 -g3 -fPIC ./test/varr_exp.test.c -I./include/ -c -o ./build/varr_exp.test.o -Werror
+ gcc -D_POSIX_C_SOURCE=200112L -std=c99 -O0 -g3 -fPIC ./test/varr_general_bound_linbuf.test.c -I./include/ -c -o ./build/varr_general_bound_linbuf.test.o -Werror
+ gcc -D_POSIX_C_SOURCE=200112L -std=c99 -O0 -g3 -fPIC ./test/varr_sequence_analysis.test.c -I./include/ -c -o ./build/varr_sequence_analysis.test.o -Werror
+ gcc -D_POSIX_C_SOURCE=200112L -std=c99 -O0 -g3 -fPIC ./test/test_results.c -I./include/ -c -o ./build/test_results.o -Werror
+ gcc -D_POSIX_C_SOURCE=200112L -std=c99 -O0 -g3 -fPIC ./test/timings_double.c -I./include/ -c -o ./build/timings_double.o -Werror
+ gcc -D_POSIX_C_SOURCE=200112L -std=c99 -O0 -g3 -fPIC ./test/timings_complex.c -I./include/ -c -o ./build/timings_complex.o -Werror
+ gcc -D_POSIX_C_SOURCE=200112L -std=c99 -O0 -g3 -fPIC ./test/sequence_generation.c -I./include/ -c -o ./build/sequence_generation.o -Werror
+ gcc -D_POSIX_C_SOURCE=200112L -std=c99 -O0 -g3 -fPIC ./test/sequence_comparison.c -I./include/ -c -o ./build/sequence_comparison.o -Werror
++ pkg-config --libs gsl
+ gcc ./build/varr_extimer.o ./build/varr_exp.o ./build/varr_log.o ./build/varr_sin.o ./build/varr_cos.o ./build/varr_phasor.o ./build/varr_sixth_root.o ./build/varr_3_over_4.o ./build/varr_atan.o ./build/varr_general_bound_linbuf.o -Werror --shared -o ./build/libvarr.so -lgsl -lgslcblas -lm
+ gcc -D_POSIX_C_SOURCE=200112L -std=c99 -O0 -g3 ./test/varr_utils.c -I./include/ -c -o ./build/varr_utils.o -Werror
++ pkg-config --libs gsl
+ gcc -O0 -g3 ./build/entry.o ./build/varr_utils.o ./build/varr_sixthroot.test.o ./build/varr_3_over_4.test.o ./build/varr_exp.test.o ./build/varr_phasor.test.o ./build/varr_atan.test.o ./build/varr_sin.test.o ./build/varr_log.test.o ./build/test_results.o ./build/timings_double.o ./build/timings_complex.o ./build/sequence_comparison.o ./build/sequence_generation.o ./build/varr_sequence_analysis.test.o ./build/varr_general_bound_linbuf.test.o -L./build/ -lvarr -o ./build/test -lgsl -lgslcblas -lm -lrt
```

Your compiler (eg. `gcc`, `icc`) and build flags (eg. `-O3`) can be customized using simple options in the above build script.

## Compatible Compilers

`Intel` `icc` (`2017-2`, `17.0.2 20170213`) and `gcc` `6.3.0` (`Debian 6.3.0-18+deb9u1 6.3.0 20170516`) compilers are both known to be compatible with this library.

## Build Dependencies

The `GSL` library is required (version `2.2.0` or higher is recommended) to compile.

Include (`-I`) and link (`-L/l`) a local `GSL` installation or module if GSL configuration is not provided by pkg-config.

## `AVX` Extensions

This library includes a number of optional features, such as some `AVX-256` and `AVX-512` vectorized `VARR` functions.

These features are typically enabled for compilation by preprocessor directives in the header file `varr.h` (of `include/`).

To compile this library with `AVX-256` features enabled, `#define __VARR_HAS_AVX__` in `varr.h` or as compiler arguments.  Add the appropriate vector compiler instruction flags to your build script (eg. `-mavx`).

Similarly to compile with `AVX-512` features enabled, `#define __VARR_USE_AVX512__` in `varr.h` or as compiler arguments.  Add the appropriate vector instruction flags to your build script.

# Tests

`varr-numerics` includes a number of self tests, which are compiled by `make.sh`.  

Upon successful compilation execute `./test` in the `build/` subdirectory.  This will perform a number of self tests, report the success or failure of each, and will report on the overall numerical precision of the library on your system:

```
/build$ ./test
general bound linbuf numerical tests:
-------- Test 1 Begins --------
Scalar evaluation:
Timing: machine stock evaluation: 1910.9
Timing: custom evaluation: 441.707
Worst (relative) numerical difference: 5.65787e-12 (50000000 tests in the range [0, 6.28319], linear sampling)
Vector evaluation:
Timing: machine stock evaluation: 1879.25
Timing: custom evaluation: 163.919
Worst (relative) numerical difference: 5.65787e-12 (50000000 tests in the range [0, 6.28319], linear sampling)
-------- Test 1 Ends --------
-------- Test 2 Begins --------
Scalar evaluation:
Timing: machine stock evaluation: 1907.91
Timing: custom evaluation: 441.383
Worst (relative) numerical difference: 5.43016e-12 (50000003 tests in the range [0, 6.28319], linear sampling)
Vector evaluation:
Timing: machine stock evaluation: 1875.01
Timing: custom evaluation: 162.373
Worst (relative) numerical difference: 5.43016e-12 (50000003 tests in the range [0, 6.28319], linear sampling)
-------- Test 2 Ends --------
-------- Test 3 Begins --------
Scalar evaluation:
Timing: machine stock evaluation: 1906.05
Timing: custom evaluation: 441.621
Worst (relative) numerical difference: 5.43016e-12 (50000003 tests in the range [0, 6.28319], linear sampling)
Vector evaluation:
Timing: machine stock evaluation: 1790.85
Timing: custom evaluation: 83.863
Worst (relative) numerical difference: 5.43016e-12 (50000003 tests in the range [0, 6.28319], linear sampling)
-------- Test 3 Ends --------
-------- Test 4 Begins --------
pow(x, 3/4) numerical tests:
Timing: machine stock evaluation: 5762.5
Timing: custom evaluation: 1138.24
Worst (relative) numerical difference: 1.08802e-14 (100000000 tests in the range [1e-18, 1e+18], logarithmic sampling)
-------- Test 4 Ends --------
pow(x, 1/6) numerical tests:
-------- Test 5 Begins --------
Timing: machine stock evaluation: 2807.64
Timing: custom evaluation: 583.439
Worst (relative) numerical difference: 8.10463e-15 (50000000 tests in the range [1e-18, 1e+18], logarithmic sampling)
-------- Test 5 Ends --------
-------- Test 6 Begins --------
Vector evaluation:
Timing: machine stock evaluation: 2552.67
Timing: custom evaluation: 436.972
Worst (relative) numerical difference: 8.10463e-15 (50000000 tests in the range [1e-18, 1e+18], logarithmic sampling)
-------- Test 6 Ends --------
-------- Test 7 Begins --------
Vector evaluation:
Timing: machine stock evaluation: 2556.09
Timing: custom evaluation: 437.974
Worst (relative) numerical difference: 8.10463e-15 (50000003 tests in the range [1e-18, 1e+18], logarithmic sampling)
-------- Test 7 Ends --------
-------- Test 8 Begins --------
Vector evaluation:
Timing: machine stock evaluation: 2468.68
Timing: custom evaluation: 353.476
Worst (relative) numerical difference: 8.10463e-15 (50000003 tests in the range [1e-18, 1e+18], logarithmic sampling)
-------- Test 8 Ends --------
-------- Test 9 Begins --------
Timing: machine stock evaluation: 2807.52
Timing: custom evaluation: 517.087
Worst (relative) numerical difference: 1.11087e-07 (50000000 tests in the range [1e-18, 1e+18], logarithmic sampling)
-------- Test 9 Ends --------
-------- Test 10 Begins --------
Vector evaluation:
Timing: machine stock evaluation: 2475.37
Timing: custom evaluation: 197.924
Worst (relative) numerical difference: 1.11077e-07 (50000003 tests in the range [1e-18, 1e+18], logarithmic sampling)
-------- Test 10 Ends --------
-------- Test 11 Begins --------
sequence analysis numerical tests (linear generation):
Numerical discrepancy: observed: 0, maximum allowed: 1e-14
-------- Test 11 Ends --------
-------- Test 12 Begins --------
Sequence analysis numerical tests (logarithmic generation):
Numerical discrepancy: observed: 5.32459e-11, maximum allowed: 1e-10
-------- Test 12 Ends --------
-------- Test 13 Begins --------
Sequence comparison numerical tests (exponential weight):
Numerical discrepancy: observed: 2.22045e-16, maximum allowed: 5e-16
-------- Test 13 Ends --------
-------- Test 14 Begins --------
log(x) numerical tests:
Case 1: normalizing sampling log:
Timing: machine stock evaluation: 212.716
Timing: custom evaluation: 153.464
Worst (relative) numerical difference: 9.62256e-10 (10000000 tests in the range [1e-10, 1e+10], logarithmic sampling)
Case 2: quad series convergent:
Timing: machine stock evaluation: 212.817
Timing: custom evaluation: 1527.6
Worst (relative) numerical difference: 1.11562e-12 (10000000 tests in the range [0.001, 1000], logarithmic sampling)
Case 3: batch normalizing sampling log:
Timing: machine stock evaluation: 989.111
Timing: custom evaluation: 342.671
Worst (relative) numerical difference: 9.62258e-10 (50000000 tests in the range [1e-10, 1e+10], logarithmic sampling)
-------- Test 14 Ends --------
-------- Test 15 Begins --------
log(x) (sublinear) numerical tests:
Case 4: batch sublinear sampling log:
Timing: machine stock evaluation: 206.722
Timing: custom evaluation: 259.45
Worst (relative) numerical difference: 1.55169e-05 (10000000 tests in the range [1e-10, 1e+10], logarithmic sampling)
Case 5: batch sublinear sampling log:
Timing: machine stock evaluation: 935.304
Timing: custom evaluation: 408.376
Worst (relative) numerical difference: 1.5524e-05 (50000000 tests in the range [1e-10, 1e+10], logarithmic sampling)
-------- Test 15 Ends --------
arctangent numerical tests:
-------- Test 16 Begins --------
Scalar evaluation:
Timing: machine stock evaluation: 973.441
Timing: custom evaluation: 467.018
Worst (relative) numerical difference: 2.80512e-10 (50000000 tests in the range [-50, 50], linear sampling)
Vector evaluation:
Timing: machine stock evaluation: 953.465
Timing: custom evaluation: 174.418
Worst (relative) numerical difference: 2.80512e-10 (50000000 tests in the range [-50, 50], linear sampling)
-------- Test 16 Ends --------
-------- Test 17 Begins --------
Scalar evaluation:
Timing: machine stock evaluation: 984.815
Timing: custom evaluation: 466.065
Worst (relative) numerical difference: 2.80631e-10 (50000003 tests in the range [-50, 50], linear sampling)
Vector evaluation:
Timing: machine stock evaluation: 953.199
Timing: custom evaluation: 175.339
Worst (relative) numerical difference: 2.80631e-10 (50000003 tests in the range [-50, 50], linear sampling)
-------- Test 17 Ends --------
-------- Test 18 Begins --------
Vector evaluation:
Timing: machine stock evaluation: 869.761
Timing: custom evaluation: 95.69
Worst (relative) numerical difference: 2.80631e-10 (50000003 tests in the range [-50, 50], linear sampling)
-------- Test 18 Ends --------
phasor(phi) = exp(i*phi) numerical tests:
-------- Test 19 Begins --------
Scalar evaluation:
Timing: machine stock evaluation: 618.991
Timing: custom evaluation: 202.398
Worst (relative) numerical difference: 5.48674e-13 (10000000 tests in the range [-18.8496, 18.8496], linear sampling)
Vector evaluation:
Timing: machine stock evaluation: 593.46
Timing: custom evaluation: 127.986
Worst (relative) numerical difference: 5.48673e-13 (10000000 tests in the range [-18.8496, 18.8496], linear sampling)
-------- Test 19 Ends --------
-------- Test 20 Begins --------
Scalar evaluation:
Timing: machine stock evaluation: 617.906
Timing: custom evaluation: 203.372
Worst (relative) numerical difference: 5.48673e-13 (10000003 tests in the range [-18.8496, 18.8496], linear sampling)
Vector evaluation:
Timing: machine stock evaluation: 594.333
Timing: custom evaluation: 128.5
Worst (relative) numerical difference: 5.48679e-13 (10000003 tests in the range [-18.8496, 18.8496], linear sampling)
-------- Test 20 Ends --------
-------- Test 21 Begins --------
sin(x) numerical tests:
Timing: machine stock evaluation: 270.738
Timing: custom evaluation: 94.69
Worst (relative) numerical difference: 2.01588e-11 (10000000 tests in the range [0, 6.28319], linear sampling)
-------- Test 21 Ends --------
exp(x) numerical tests:
-------- Test 22 Begins --------
Scalar evaluation:
Timing: machine stock evaluation: 863.034
Timing: custom evaluation: 517.708
Worst (relative) numerical difference: 1.43219e-14 (50000000 tests in the range [-10, 10], linear sampling)
Vector evaluation:
Timing: machine stock evaluation: 845.31
Timing: custom evaluation: 193.939
Worst (relative) numerical difference: 1.43219e-14 (50000000 tests in the range [-10, 10], linear sampling)
-------- Test 22 Ends --------
-------- Test 23 Begins --------
Scalar evaluation:
Timing: machine stock evaluation: 862.228
Timing: custom evaluation: 518.017
Worst (relative) numerical difference: 1.43219e-14 (50000003 tests in the range [-10, 10], linear sampling)
Vector evaluation:
Timing: machine stock evaluation: 845.945
Timing: custom evaluation: 192.577
Worst (relative) numerical difference: 1.43219e-14 (50000003 tests in the range [-10, 10], linear sampling)
-------- Test 23 Ends --------
-------- Test 24 Begins --------
Vector evaluation:
Timing: machine stock evaluation: 764.505
Timing: custom evaluation: 112.637
Worst (relative) numerical difference: 1.43219e-14 (50000003 tests in the range [-10, 10], linear sampling)
-------- Test 24 Ends --------
-------- Result --------
Number of tests evaluated: 24
Number of failed tests: 0
Notes and messages:
varr general bound linbuf (cos specialization): Sampling evaluation (no remainder loop): Pass: observed numerical error: 5.65787e-12; maximum allowed: 5.66e-12. Threshold missed by 0.0376991%.
varr general bound linbuf (cos specialization): Sampling evaluation (with remainder loop): Pass: observed numerical error: 5.43016e-12; maximum allowed: 5.66e-12. Threshold missed by 4.06073%.
varr general bound linbuf (cos specialization): Sampling evaluation (with remainder loop, in-place): Pass: observed numerical error: 5.43016e-12; maximum allowed: 5.66e-12. Threshold missed by 4.06073%.
varr-3_over_4/pow(x, 3/4): Sampling evaluation: Pass: observed numerical error: 1.08802e-14; maximum allowed: 1.12e-14. Threshold missed by 2.85549%.
varr-sixthroot/pow(x, 1/6): Sampling evaluation: Pass: observed numerical error: 8.10463e-15; maximum allowed: 8.66e-15. Threshold missed by 6.41307%.
varr-sixthroot/pow(x, 1/6) (batch evaluation): Sampling evaluation: Pass: observed numerical error: 8.10463e-15; maximum allowed: 8.66e-15. Threshold missed by 6.41307%.
varr-sixthroot/pow(x, 1/6) (batch evaluation, with remainder loop): Sampling evaluation: Pass: observed numerical error: 8.10463e-15; maximum allowed: 8.66e-15. Threshold missed by 6.41307%.
varr-sixthroot/pow(x, 1/6) (batch evaluation, with remainder loop, in-place): Sampling evaluation: Pass: observed numerical error: 8.10463e-15; maximum allowed: 8.66e-15. Threshold missed by 6.41307%.
varr-sixthroot/pow(x, 1/6) (scalar evaluation, with remainder loop, in-place, sublinear): Sampling evaluation: Pass: observed numerical error: 1.11087e-07; maximum allowed: 1.12e-07. Threshold missed by 0.815336%.
varr-sixthroot/pow(x, 1/6) (batch evaluation, with remainder loop, in-place, sublinear): Sampling evaluation: Pass: observed numerical error: 1.11077e-07; maximum allowed: 1.12e-07. Threshold missed by 0.823726%.
varr-sequence-analysis: sequence generation (linear): Pass: observed numerical error: 0; maximum allowed: 1e-14. Threshold missed by 100%.
varr-sequence-analysis: sequence generation (logscale): Pass: observed numerical error: 5.32459e-11; maximum allowed: 1e-10. Threshold missed by 46.7541%.
varr-sequence-analysis: sequence comparison (exponential): Pass: observed numerical error: 2.22045e-16; maximum allowed: 5e-16. Threshold missed by 55.5911%.
varr-log/natural-logarithm: Sublinear sampling evaluation: Pass: observed numerical error: 1.5524e-05; maximum allowed: 1.56e-05. Threshold missed by 0.487075%.
varr-atan/arctangent: Sampling evaluation (no remainder loop): Pass: observed numerical error: 2.80512e-10; maximum allowed: 2.82e-10. Threshold missed by 0.527781%.
varr-atan/arctangent: Sampling evaluation (with remainder loop): Pass: observed numerical error: 2.80631e-10; maximum allowed: 2.82e-10. Threshold missed by 0.485562%.
varr-atan/arctangent: Sampling evaluation (with remainder loop, in-place): Pass: observed numerical error: 2.80631e-10; maximum allowed: 2.82e-10. Threshold missed by 0.485562%.
varr-phasor/phasor(phi) = exp(i*phi): Sampling evaluation (no remainder loop): Pass: observed numerical error: 5.48674e-13; maximum allowed: 5.5e-13. Threshold missed by 0.241162%.
varr-phasor/phasor(phi) = exp(i*phi): Sampling evaluation (with remainder loop): Pass: observed numerical error: 5.48679e-13; maximum allowed: 5.5e-13. Threshold missed by 0.240111%.
varr-sin/sin(x): Sampling evaluation: Pass: observed numerical error: 2.01588e-11; maximum allowed: 2.1e-11. Threshold missed by 4.00557%.
varr-exp/exp(x): Sampling evaluation (no remainder loop): Pass: observed numerical error: 1.43219e-14; maximum allowed: 1.47e-14. Threshold missed by 2.57227%.
varr-exp/exp(x): Sampling evaluation (with remainder loop): Pass: observed numerical error: 1.43219e-14; maximum allowed: 1.47e-14. Threshold missed by 2.57227%.
varr-exp/exp(x): Sampling evaluation (with remainder loop, in-place): Pass: observed numerical error: 1.43219e-14; maximum allowed: 1.47e-14. Threshold missed by 2.57227%.
```

# `VARR` Functions

This library provides a number of primitive `VARR` functions of various kinds.

These functions are grouped according to their mathematical classes.  For example, `VARR` implementations of the real arctangent function are declared by the header file `varr_atan.h`.

## `VARR` Clamped Real Arctangent

The header file `varr_atan.h` declares a `VARR` implementation of the real arctangent function (hereafter `atan(x)`).

This implementation uses sampling techniques to compute the arctangent.  It allocates (and maintains references to) memory for this purpose.

The number of sampling points requested of this implementation is customizable, and the amount of memory allocated by it is approximately proportional to the number of sampling points requested.

The numerical accuracy of this implementation generally increases with the number of sampling points requested.

### Allocation

This `VARR` function is allocated and disallocated as per the following example:

```c++
#include "varr_atan.h"

// Allocate a VARR sampling arctangent function with 100,000 sampling points:
size_t const
   number_of_samples = 100000u;
VARRAtanDEvaluator
   evaluator = clamping_linear_interpolating_atand(number_of_samples);
// Ultimately, disallocate the same:
evaluator.disallocate(&evaluator);
```

The above method `clamping_linear_interpolating_atand` allocates and populates a grid of real (double precision) arctangent values.  The size of this sampling grid is indicated by its argument, `number_of_samples`, which must be nonzero.

The `disallocate` function above must be called once (once per object as indicated) when this `VARR` function is to be disallocated, never to be used by the application again.  Application memory leaks and/or undefined behaviour may ultimately result if this is not done.  The above object must not be used after it is disallocated.

The `evaluator` object above provides scalar and vector `VARR` arctangent functions as follows:

### Scalar `VARR` Arctangent

```c++
double
   x = 1.0;
// Compute an approximate value of atan(x) using the above VARR evaluator (scalar case):
double const
   atan_x = evaluator.atan(x, evaluator.accelerator);
```

If `|x| > 50` then this implementation clamps `x` into the range `-50 <= x <= 50.`.  This evaluator does not return meaningful values if its inputs `x` are not finite real numbers. 

### Vector `VARR` Arctangent

```c++
{
// Compute the approximate values of atan(x) for an array of values x, using the above 
// VARR evaluator (vector case):
double const
   x[] = { -2., -1., 0., 1., 2. };        // input array
double
   atan_x[5];                             // output array
size_t
   length = 5u;
evaluator.atan_array(
   x, atan_x,
   length,
   evaluator.accelerator
   );                                     // atan_x[i] contains the approximated
                                          // values atan(x[i]), i < 5
}
```

If `|x| > 50` then this implementation clamps `x` into the range `-50 <= x <= 50.`.  This evaluator does not return meaningful values if its inputs `x` are not finite real numbers. 

It is the responsibility of the caller to ensure that the input array (`x`) and the output array (`atan_x`) both have at least the length indicated by `length`, and that both of these arguments are valid non-null pointers.

It is not necessary for the arrays `x` nor `atan_x` to have any specific byte alignments.  If this library is compiled without `AVX` vectorizing extensions enabled, this method will delegate to a scalar `VARR` arctangent function.

## `VARR` Phasor

The header file `varr_phasor.h` declares a `VARR` implementation of the complex phasor function (hereafter `cexp(x)`).  

The values of this function are unit norm complex numbers of the form `exp(i*phi)`, parametrized by angles `phi` in the complex plane.

This implementation uses sampling techniques to compute the complex phasor.  It allocates (and maintains references to) memory for this purpose.

The number of sampling points requested of this implementation is customizable, and the amount of memory allocated by it is approximately proportional to the number of sampling points requested.

The numerical accuracy of this implementation generally increases with the number of sampling points requested.

This library provides several `VARR` phasor implementations.  Chief among these implementations is the vectorizing linear sampling phasor as follows:

### `VARR` Linear Sampling Phasor

#### Allocation

This `VARR` function is allocated and disallocated as per the following example:

```c++
#include "varr_phasor.h"

// Allocate a VARR phasor function with 100,000 sampling points:
size_t const
   number_of_samples = 100000u;
VARRPhasorDEvaluator
   evaluator = linear_interpolating_phasord(number_of_samples);
// Ultimately, disallocate the same:
evaluator.disallocate(&evaluator);
```

The above method `linear_interpolating_phasord` allocates and populates a grid of double complex phasors.  The size of this sampling grid is indicated by its argument, `number_of_samples`, which must be nonzero.

The `disallocate` function above must be called once (once per object as indicated) when this `VARR` function is to be disallocated, never to be used by the application again.  Application memory leaks and/or undefined behaviour may ultimately result if this is not done.  The above object must not be used after it is disallocated.

The `evaluator` object above provides scalar and vector `VARR` phasor functions as follows:

#### Scalar Form

```c++
double
   x = 1.0;
// Compute an approximate value of cexp(i*x) using the above VARR evaluator (scalar case):
double const
   cexp_x = evaluator.phasord(x, evaluator.accelerator);
```

This evaluator does not return meaningful values if its inputs `x` are not finite real numbers. 

It it not necessary for the angle `x` to be normalized or restricted to (eg.) the range `0 <= x <= 2pi`.  Any finite signed values of `x` are acceptable.

#### Vector Form

```c++
{
// Compute the approximate values of cexp(i*x) for an array of values x, using the above 
// VARR evaluator (vector case):
double const
   x[] = { -2., -1., 0., 1., 2. };        // input array
double complex
   cexp_x[5];                             // output array
size_t
   length = 5u;
evaluator.phasord_array(
   x, cexp_x,
   length,
   evaluator.accelerator
   );                                     // cexp_x[i] contains the approximated
                                          // values cexp(i*x[i]), i < 5
}
```

It is the responsibility of the caller to ensure that the input array (`x`) and the output array (`cexp_x`) both have at least the length indicated by `length`, and that both of these arguments are valid non-null pointers.  Note for allocation purposes that `sizeof(double)` is not equal to `sizeof(double complex)` in general.

It is not necessary for the arrays `x` nor `cexp_x` to have any specific byte alignments.  If this library is compiled without `AVX` vectorizing extensions enabled, this method will delegate to a scalar `VARR` phasor function.

## `VARR` Real Exponential

The header file `varr_exp.h` declares a `VARR` implementation of the (real-valued, real-argument) exponential function (hereafter `exp(x)`).  

This implementation uses sampling techniques to compute the real exponential.  It allocates (and maintains references to) memory for this purpose.

The number of sampling points requested of this implementation is customizable, and the amount of memory allocated by it is approximately proportional to the number of sampling points requested.

The numerical accuracy of this implementation generally increases with the number of sampling points requested.

### `VARR` Linear Sampling Real Exponential

#### Allocation

This `VARR` function is allocated and disallocated as per the following example:

```c++
#include "varr_exp.h"

// Allocate a VARR phasor function with 100,000 sampling points:
size_t const
   number_of_samples = 100000u;
VARRExpDEvaluator
   evaluator = shifting_linear_sampling_expd(number_of_samples);
// Ultimately, disallocate the same:
evaluator.disallocate(&evaluator);
```

The above method `shifting_linear_sampling_expd` allocates and populates a grid of real double exponential samples.  The size of this sampling grid is indicated by its argument, `number_of_samples`, which must be nonzero.

The `disallocate` function above must be called once (once per object as indicated) when this `VARR` function is to be disallocated, never to be used by the application again.  Application memory leaks and/or undefined behaviour may ultimately result if this is not done.  The above object must not be used after it is disallocated.

The `evaluator` object above provides scalar and vector `VARR` exponential functions as follows:

#### Scalar Form

```c++
double
   x = 1.0;
// Compute an approximate value of exp(x) using the above VARR evaluator (scalar case):
double const
   exp_x = evaluator.expd(x, evaluator.accelerator);
```

This function does not return meaningful values if the input, `x`, is not a finite real number or if `x` is outside of the range `-1024 <= x <= +709`.

#### Vector Form

```c++
{
// Compute the approximate values of exp(x) for an array of values x, using the above 
// VARR evaluator (vector case):
double const
   x[] = { -2., -1., 0., 1., 2. };        // input array
double complex
   exp_x[5];                              // output array
size_t
   length = 5u;
evaluator.expd_array(
   x, exp_x,
   length,
   evaluator.accelerator
   );                                     // exp_x[i] contains the approximated
                                          // values exp(i*x[i]), i < 5
}
```

It is the responsibility of the caller to ensure that the input array (`x`) and the output array (`exp_x`) both have at least the length indicated by `length`, and that both of these arguments are valid non-null pointers.

It is not necessary for the arrays `x` nor `exp_x` to have any specific byte alignments.  If this library is compiled without `AVX` vectorizing extensions enabled, this method will delegate to a scalar `VARR` phasor function.

This function does not return meaningful values if the inputs, `x`, are not finite real numbers or if any such `x` is outside of the range `-1024 <= x <= +709`.

## `VARR` Real `Sine` Function

The header file `varr_sin.h` declares a `VARR` implementation of the real trigonometric `sine` function (hereafter `sin(x)`).  

This library provides several `VARR` `sine` implementations.  Chief among these implementations is the vectorizing linear sampling `sine` function as follows:

### `VARR` Linear Sampling Real `Sine` Function

This implementation uses sampling techniques to compute the real `sine` function.  It allocates (and maintains references to) memory for this purpose.

The number of sampling points requested of this implementation is customizable, and the amount of memory allocated by it is approximately proportional to the number of sampling points requested.

The numerical accuracy of this implementation generally increases with the number of sampling points requested.

#### Allocation

This `VARR` function is allocated and disallocated as per the following example:

```c++
#include "varr_sin.h"

// Allocate a VARR sine function with 100,000 sampling points:
size_t const
   number_of_samples = 100000u;
VARRSinDEvaluator
   evaluator = sampling_sind(number_of_samples);
// Ultimately, disallocate the same:
evaluator.disallocate(&evaluator);
```

The above method `sampling_sind` allocates and populates a grid of real double `sine` samples.  The size of this sampling grid is indicated by its argument, `number_of_samples`, which must be nonzero.

The `disallocate` function above must be called once (once per object as indicated) when this `VARR` function is to be disallocated, never to be used by the application again.  Application memory leaks and/or undefined behaviour may ultimately result if this is not done.  The above object must not be used after it is disallocated.

The `evaluator` object above provides a scalar `VARR` `sine` functions as follows:

#### Scalar Form

```c++
double
   x = 1.0;
// Compute an approximate value of sin(x) using the above VARR evaluator (scalar case):
double const
   sin_x = evaluator.sind(x, evaluator.accelerator);
```

This function does not return meaningful values if the input, `x`, is not a finite real number.

### `VARR` Cubic Spline `Sine` Function

This implementation uses cubic spline techniques to compute the real `sine` function.  It allocates (and maintains references to) memory for this purpose.

The number of interpolation samples requested of this implementation is customizable, and the amount of memory allocated by it is approximately proportional to the number of samples requested.

The numerical accuracy of this implementation generally increases with the number of interpolation points requested.

#### Allocation

This `VARR` function is allocated and disallocated as per the following example:

```c++
#include "varr_sin.h"

// Allocate a VARR sine function with 100,000 sampling points:
size_t const
   number_of_samples = 100000u;
VARRSinDEvaluator
   evaluator = cubic_spline_sampling_sind(number_of_samples);
// Ultimately, disallocate the same:
evaluator.disallocate(&evaluator);
```

The above method `cubic_spline_sampling_sind` allocates and populates a cubic spline system of sampled real double `sine` samples.  The size of this sampling grid is indicated by its argument, `number_of_samples`, which must be nonzero.

The `disallocate` function above must be called once (once per object as indicated) when this `VARR` function is to be disallocated, never to be used by the application again.  Application memory leaks and/or undefined behaviour may ultimately result if this is not done.  The above object must not be used after it is disallocated.

The `evaluator` object above provides a scalar `VARR` `sine` functions as follows:

#### Scalar Form

```c++
double
   x = 1.0;
// Compute an approximate value of sin(x) using the above VARR evaluator (scalar case):
double const
   sin_x = evaluator.sind(x, evaluator.accelerator);
```

This function does not return meaningful values if the input, `x`, is not a finite real number or if `x` is outside of the range `0 <= x <= 2*pi`.

## `VARR` Real `Cosine` Function

The header file `varr_cos.h` declares a `VARR` implementation of the real trigonometric `cosine` function (hereafter `cos(x)`).  

This library provides several `VARR` `cosine` implementations.  Chief among these implementations is the vectorizing linear sampling `cosine` function as follows:

### `VARR` Linear Sampling Real `Cosine` Function

This implementation uses sampling techniques to compute the real `cosine` function.  It allocates (and maintains references to) memory for this purpose.

The number of sampling points requested of this implementation is customizable, and the amount of memory allocated by it is approximately proportional to the number of sampling points requested.

The numerical accuracy of this implementation generally increases with the number of sampling points requested.

#### Allocation

This `VARR` function is allocated and disallocated as per the following example:

```c++
#include "varr_cos.h"

// Allocate a VARR cosine function with 100,000 sampling points:
size_t const
   number_of_samples = 100000u;
VARRCosDEvaluator
   evaluator = sampling_cosd(number_of_samples);
// Ultimately, disallocate the same:
evaluator.disallocate(&evaluator);
```

The above method `sampling_cosd` allocates and populates a grid of real double `cosine` samples.  The size of this sampling grid is indicated by its argument, `number_of_samples`, which must be nonzero.

The `disallocate` function above must be called once (once per object as indicated) when this `VARR` function is to be disallocated, never to be used by the application again.  Application memory leaks and/or undefined behaviour may ultimately result if this is not done.  The above object must not be used after it is disallocated.

The `evaluator` object above provides a scalar `VARR` `cosine` functions as follows:

#### Scalar Form

```c++
double
   x = 1.0;
// Compute an approximate value of cos(x) using the above VARR evaluator (scalar case):
double const
   cos_x = evaluator.cosd(x, evaluator.accelerator);
```

This function does not return meaningful values if the input, `x`, is not a finite real number.

### `VARR` Cubic Spline `Cosine` Function

This implementation uses cubic spline techniques to compute the real `cosine` function.  It allocates (and maintains references to) memory for this purpose.

The number of interpolation samples requested of this implementation is customizable, and the amount of memory allocated by it is approximately proportional to the number of samples requested.

The numerical accuracy of this implementation generally increases with the number of interpolation points requested.

#### Allocation

This `VARR` function is allocated and disallocated as per the following example:

```c++
#include "varr_sin.h"

// Allocate a VARR cosine function with 100,000 sampling points:
size_t const
   number_of_samples = 100000u;
VARRCosDEvaluator
   evaluator = cubic_spline_sampling_cosd(number_of_samples);
// Ultimately, disallocate the same:
evaluator.disallocate(&evaluator);
```

The above method `cubic_spline_sampling_cosd` allocates and populates a cubic spline system of sampled real double `cosine` samples.  The size of this sampling grid is indicated by its argument, `number_of_samples`, which must be nonzero.

The `disallocate` function above must be called once (once per object as indicated) when this `VARR` function is to be disallocated, never to be used by the application again.  Application memory leaks and/or undefined behaviour may ultimately result if this is not done.  The above object must not be used after it is disallocated.

The `evaluator` object above provides a scalar `VARR` `cosine` functions as follows:

#### Scalar Form

```c++
double
   x = 1.0;
// Compute an approximate value of cos(x) using the above VARR evaluator (scalar case):
double const
   cos_x = evaluator.cosd(x, evaluator.accelerator);
```

This function does not return meaningful values if the input, `x`, is not a finite real number or if `x` is outside of the range `0 <= x <= 2*pi`.

## `VARR` Real Natural Logarithm Function

The header file `varr_log.h` declares a `VARR` implementation of the real natural logarithm function (hereafter `log(x)`).

This library provides several `VARR` `log` implementations.  Chief among these implementations is the vectorizing linear sampling `log` function as follows:

### `VARR` Linear Sampling Real `Log` Function

This implementation uses linear interpolation techniques combined with sampling and a modified De-Bruijn-like method to compute the real natural logarithm.  It allocates and maintains references to memory for this purpose.

The number of sampling points requested of this implementation is customizable, and the amount of memory allocated by it is approximately proportional to the number of sampling points requested.

The numerical accuracy of this implementation generally increases with the number of sampling points requested.

#### Allocation

This `VARR` function is allocated and disallocated as per the following example:

```c++
#include "varr_log.h"

// Allocate a VARR log function with 100,000 sampling points:
size_t const
   number_of_samples = 100000u;
VARRLogDEvaluator
   evaluator = normalizing_linear_sampling_logd(number_of_samples);
// Ultimately, disallocate the same:
evaluator.disallocate(&evaluator);
```

The above method `normalizing_linear_sampling_logd` allocates and populates a grid of real double `log` samples.  The size of this sampling grid is indicated by its argument, `number_of_samples`, which must be nonzero.

The `disallocate` function above must be called once (once per object as indicated) when this `VARR` function is to be disallocated, never to be used by the application again.  Application memory leaks and/or undefined behaviour may ultimately result if this is not done.  The above object must not be used after it is disallocated.

The `evaluator` object above provides a scalar and a vector `VARR` `log` functions as follows:

#### Scalar Form

```c++
double
   x = 1.0;
// Compute an approximate value of logd(x) using the above VARR evaluator (scalar case):
double const
   log_x = evaluator.logd(x, evaluator.accelerator);
```

This function does not return meaningful values if the input, `x`, is not a finite real number or if `x` is outside of the range `+2**-63 <= x <= +2**+63`.

#### Vector Form

```c++
{
// Compute the approximate values of log(x) using the above VARR evaluator (vector case):
double const
   x[] = { 1., 2., 3., 4., 5. };          // input array
double complex
   log_x[5];                              // output array
size_t
   length = 5u;
evaluator.logd_array(
   x, log_x,
   length,
   evaluator.accelerator
   );                                     // log_x[i] contains the approximated
                                          // values log(x[i]), i < 5
}
```

It is the responsibility of the caller to ensure that the input array (`x`) and the output array (`log_x`) both have at least the length indicated by `length`, and that both of these arguments are valid non-null pointers.

It is not necessary for the arrays `x` nor `log_x` to have any specific byte alignments.  If this library is compiled without `AVX` vectorizing extensions enabled, this method will delegate to a scalar `VARR` logarithm function.

This function does not return meaningful values if the inputs, `x`, are not finite real numbers or if any such `x` is outside of the range `+2**-63 <= x <= +2**+63`.

### `VARR` Sublinear Sampling Real `Log` Function

This implementation uses sublinear binning techniques combined with sampling to compute the real natural logarithm.  It allocates and maintains references to memory for this purpose.

The number of sampling points requested of this implementation is customizable, and the amount of memory allocated by it is approximately proportional to the number of sampling points requested.

The numerical accuracy of this implementation generally increases with the number of sampling points requested.  The accuracy of this implementation is typically significantly lower than its linear-interpolating peer for the same number of sampling points requested, but this implementation may be faster in some use cases.

#### Allocation

This `VARR` function is allocated and disallocated as per the following example:

```c++
#include "varr_log.h"

// Allocate a sublinear VARR log function with 100,000 sampling points:
size_t const
   number_of_samples = 100000u;
VARRLogDEvaluator
   evaluator = normalizing_sublinear_sampling_logd(number_of_samples);
// Ultimately, disallocate the same:
evaluator.disallocate(&evaluator);
```

The above method `normalizing_sublinear_sampling_logd` allocates and populates a grid of real double `log` samples.  The size of this sampling grid is indicated by its argument, `number_of_samples`, which must be nonzero.

The `disallocate` function above must be called once (once per object as indicated) when this `VARR` function is to be disallocated, never to be used by the application again.  Application memory leaks and/or undefined behaviour may ultimately result if this is not done.  The above object must not be used after it is disallocated.

The `evaluator` object above provides a scalar and a vector sublinear `VARR` `log` functions as follows:

#### Scalar Form

```c++
double
   x = 1.0;
// Compute an approximate value of logd(x) using the above VARR evaluator (scalar case):
double const
   log_x = evaluator.logd(x, evaluator.accelerator);
```

This function does not return meaningful values if the input, `x`, is not a finite real number or if `x` is outside of the range `+2**-63 <= x <= +2**+63`.

#### Vector Form

```c++
{
// Compute the approximate values of log(x) using the above VARR evaluator (vector case):
double const
   x[] = { 1., 2., 3., 4., 5. };          // input array
double complex
   log_x[5];                              // output array
size_t
   length = 5u;
evaluator.logd_array(
   x, log_x,
   length,
   evaluator.accelerator
   );                                     // log_x[i] contains the approximated
                                          // values log(x[i]), i < 5
}
```

It is the responsibility of the caller to ensure that the input array (`x`) and the output array (`log_x`) both have at least the length indicated by `length`, and that both of these arguments are valid non-null pointers.

It is not necessary for the arrays `x` nor `log_x` to have any specific byte alignments.  If this library is compiled without `AVX` vectorizing extensions enabled, this method will delegate to a scalar `VARR` logarithm function.

This function does not return meaningful values if the inputs, `x`, are not finite real numbers or if any such `x` is outside of the range `+2**-63 <= x <= +2**+63`.

### `VARR` Quad Series Real `Log` Function

This `VARR` function uses a power series to compute the natural logarithm.  The number of terms `N` comprising this power series is customizable.  

This power series evaluates `log(x)` using positive integer powers of `(y - 1)**2/(y + 1)**2` where 
`y = x/10**q` and `q` is selected so that `|y| < 1`.  The rate of convergence of this power series for any given `x` is generally proportional to the value of `N`.

This `VARR` function allocates minimal memory or no memory.  The numerical accuracy of the natural logarithm function that is generated by this method generally increases with the number of iterations requested.

#### Allocation

This `VARR` function is allocated and disallocated as per the following example:

```c++
#include "varr_log.h"

// Allocate a VARR log function using 50 iterations:
size_t const
   N = 50u;
VARRLogDEvaluator
   evaluator = quad_series_logd(N);
// Ultimately, disallocate the same:
evaluator.disallocate(&evaluator);
```

The number of iterations specified (`N`) must be nonzero and must not exceed `500`.

The `disallocate` function above must be called once (once per object as indicated) when this `VARR` function is to be disallocated, never to be used by the application again.  Application memory leaks and/or undefined behaviour may ultimately result if this is not done.  The above object must not be used after it is disallocated.

The `evaluator` object above provides a scalar `VARR` `log` functions as follows:

#### Scalar Form

```c++
double
   x = 1.0;
// Compute an approximate value of logd(x) using the above VARR evaluator (scalar case):
double const
   log_x = evaluator.logd(x, evaluator.accelerator);
```

This function does not return meaningful values if the input, `x`, is not a finite real number of if `x` is not confined to the range `10**-5 < x * < 10**+5`.

## `VARR` Real Sixth Root

The header file `varr_sixth_root.h` declares a `VARR` implementation of the real sixth root function (hereafter `x**1/6` or `pow(x, 1/6)`).

### `VARR` Linear Sampling Real Sixth Root

This implementation uses linear interpolation techniques combined with sampling and a modified De-Bruijn-like method to compute the real sixth root.  It allocates and maintains references to memory for this purpose.

The number of sampling points requested of this implementation is customizable, and the amount of memory allocated by it is approximately proportional to the number of sampling points requested.

The numerical accuracy of this implementation generally increases with the number of sampling points requested.

#### Allocation

This `VARR` function is allocated and disallocated as per the following example:

```c++
#include "varr_sixth_root.h"

// Allocate a VARR pow(x, 1/6) function with 100,000 sampling points:
size_t const
   number_of_samples = 100000u;
VARRSixthRootDEvaluator
   evaluator = linear_sampling_normalizing_sixth_rootd(number_of_samples);
// Ultimately, disallocate the same:
evaluator.disallocate(&evaluator);
```

The above method `linear_sampling_normalizing_sixth_rootd` allocates and populates a grid of real double `x**1/6` samples.  The size of this sampling grid is indicated by its argument, `number_of_samples`, which must be nonzero.

The `disallocate` function above must be called once (once per object as indicated) when this `VARR` function is to be disallocated, never to be used by the application again.  Application memory leaks and/or undefined behaviour may ultimately result if this is not done.  The above object must not be used after it is disallocated.

The `evaluator` object above provides a scalar and a vector `VARR` `x**1/6` functions as follows:

#### Scalar Form

```c++
double
   x = 1.0;
// Compute an approximate value of pow(x, 1/6) using the above VARR evaluator (scalar case):
double const
   root6_x = evaluator.sixthrootd(x, evaluator.accelerator);
```

This function does not return meaningful values if the input, `x`, is not a finite real number of if `x` is not confined to the range `10**-18 < x < 10**+18`.

#### Vector Form

```c++
{
// Compute the approximate values of pow(x, 1./6.) using the above VARR evaluator (vector case):
double const
   x[] = { 1., 2., 3., 4., 5. };          // input array
double complex
   root6_x[5];                            // output array
size_t
   length = 5u;
evaluator.sixthrootd_array(
   x, root6_x,
   length,
   evaluator.accelerator
   );                                     // root6_x[i] contains the approximated
                                          // values pow(x[i], 1./6.), i < 5
}
```

It is the responsibility of the caller to ensure that the input array (`x`) and the output array (`root6_x`) both have at least the length indicated by `length`, and that both of these arguments are valid non-null pointers.

It is not necessary for the arrays `x` nor `root6_x` to have any specific byte alignments.  If this library is compiled without `AVX` vectorizing extensions enabled, this method will delegate to a scalar `VARR` real sixth root function.

This function does not return meaningful values if the inputs, `x`, are not a finite real numbers or if `x` are not confined to the range `10**-18 < x < 10**+18`.
