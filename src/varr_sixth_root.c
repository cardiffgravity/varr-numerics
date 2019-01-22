/*
 *  This file is part of varr-numerics, an experimental variable resolution
 *  primitive numerics library.
 *
 *  Copyright (C) 2019 Cardiff University
 *
 *  varr-numerics is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  varr-numerics is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with varr-numerics.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "varr_sixth_root.h"
#include "varr.h"
#include "varr_internal.h"

#include <math.h>
#include <stdlib.h>
#include <inttypes.h>
#include <string.h>

typedef struct tagSamplingSixthRootDAccelerator
{
   double
      start_x,
      step_x,
      step_x_inverse;
   double
      step_x_inverse_powers[64];
   size_t
      samples;
   double const *
      values;
} SamplingSixthRootDAccelerator;

static
int
__linear_sampling_normalizing_sixth_rootd_disallocate(
   SamplingSixthRootDAccelerator * accelerator
   )
{
   if(accelerator == NULL)
   {
      return 1;
   }
   free((void *) accelerator->values);
   accelerator->start_x = 0x7F800001;
   accelerator->step_x = 0x7F800001;
   accelerator->step_x_inverse = 0x7F800001;
   for(unsigned i = 0u; i< 64u; ++i)
   {
      accelerator->step_x_inverse_powers[i] = 0x7F800001;
   }
   accelerator->samples = (size_t) 0;
   accelerator->values = NULL;
   free(accelerator);
   return 0;
}

static
int
linear_sampling_normalizing_sixth_rootd_disallocate(
   VARRSixthRootDEvaluator * evaluator
   )
{
   if(evaluator == NULL)
   {
      return 1;
   }
   SamplingSixthRootDAccelerator * const
      accelerator = (SamplingSixthRootDAccelerator *)
         evaluator->accelerator;
   return
      __linear_sampling_normalizing_sixth_rootd_disallocate(accelerator);
}

static
SamplingSixthRootDAccelerator *
allocate_linear_sampling_normalizing_sixth_rootd(
   size_t number_of_samples
   )
{
   double const
      step_size = 1.0 / (double) (number_of_samples - (size_t) 1u);
   double * const
      values = (double *) 
         malloc(sizeof(double) * (number_of_samples + 2u)
         );
   for(size_t i = (size_t) 0; i< (number_of_samples + 2u); ++i)
   {
      double const
         x = step_size * (double) i,
         value = pow(x, 1.0/6.0);
      values[i] = value;
      continue;
   }
   SamplingSixthRootDAccelerator * const
      result = (SamplingSixthRootDAccelerator *) malloc(
         sizeof(SamplingSixthRootDAccelerator)
         );
   result->start_x = 0.0;
   result->step_x = step_size;
   result->step_x_inverse = 1.0 / step_size;
   for(uint64_t i = 0u; i< 64u; ++i)
   {
      result->step_x_inverse_powers[i] = 
         result->step_x_inverse / (double) (((uint64_t) 1) << i);
   }
   result->samples = number_of_samples,
   result->values = values;
   return
      result;
}

#ifdef __INTEL_COMPILER
static const int
   __log2_table[128] = {
   0, 1, 59, 0, 0, 2, 60, 0, 48, 0, 54, 0, 0, 3, 61, 0, 0, 40, 49, 0, 28, 0, 
   55, 0, 34, 0, 43, 0, 0, 4, 62, 0, 0, 52, 38, 0, 0, 41, 50, 0, 0, 19, 29, 
   0, 0, 21, 56, 0, 0, 31, 35, 0, 0, 12, 44, 0, 15, 0, 23, 0, 0, 5, 63, 0, 58, 
   0, 0, 47, 53, 0, 0, 39, 0, 27, 0, 33, 42, 0, 0, 51, 37, 0, 0, 18, 0, 20, 0, 
   30, 0, 11, 0, 14, 22, 0, 0, 57, 46, 0, 0, 26, 32, 0, 0, 36, 17, 0, 0, 10, 
   13, 0, 0, 45, 25, 0, 0, 16, 9, 0, 0, 24, 0, 8, 0, 7, 6, 0, 0
   };
#endif

static const double
   // i: 0-100: (1lu << i) ** 1./6.
   __sixth_roots_of_2n[] = {
                            1,
   1.122462048309373017218604,
   1.259921049894873190666544,
   1.414213562373095145474622,
   1.587401051968199361397183,
   1.781797436280678548214951,
                            2,
   2.244924096618746034437208,
   2.519842099789746381333089,
   2.828427124746190290949244,
   3.174802103936399166883575,
   3.563594872561357096429902,
                            4,
   4.489848193237491180695997,
   5.039684199579492762666177,
   5.656854249492380581898487,
   6.349604207872797445588731,
   7.127189745122715081038223,
                            8,
   8.979696386474982361391994,
   10.07936839915898552533235,
   11.31370849898476116379697,
   12.69920841574559489117746,
   14.25437949024543016207645,
                           16,
   17.95939277294997182821135,
   20.15873679831796749795103,
   22.62741699796952232759395,
   25.39841683149119688778228,
   28.50875898049085321872553,
                           32,
   35.91878554589994365642269,
   40.31747359663593499590206,
    45.2548339959390446551879,
   50.79683366298239377556456,
   57.01751796098170643745107,
                           64,
   71.83757109179988731284539,
   80.63494719327186999180412,
    90.5096679918780893103758,
   101.5936673259647875511291,
   114.0350359219634128749021,
                          128,
   143.6751421835997746256908,
   161.2698943865437399836082,
   181.0193359837561786207516,
   203.1873346519295751022582,
   228.0700718439268257498043,
                          256,
    287.350284367199378721125,
   322.5397887730876504974731,
   362.0386719675123572415032,
    406.374669303858922830841,
   456.1401436878539357167028,
                          512,
   574.7005687343987574422499,
   645.0795775461753009949462,
   724.0773439350247144830064,
   812.7493386077178456616821,
   912.2802873757078714334057,
                         1024,
     1149.4011374687975148845,
   1290.159155092350601989892,
   1448.154687870049428966013,
   1625.498677215435691323364
   };

#include "varr_floor_log2.h"

static
double
linear_sampling_normalizing_sixth_rootd_evaluate(
   register double __x,
   register void const * restrict __accelerator
   )
{
   register char const
      do_invert = (__x < 1.0);
   register double const
      x = do_invert ? 1.0/__x : __x;
#ifdef __INTEL_COMPILER
   /*
    * An exponential split pow(x, 1/6) = pow(x / 2**k, 1/6) * 2**(k/6),
    * for the least integer k such that x < 2**k, combined with linear
    * interpolation in the range pow([0-1], 1/6) and a modified De-Bruijn-
    * like method.
    */
   register uint64_t
      __pow2_gt = (uint64_t) x;
   __pow2_gt |= __pow2_gt >> 1;
   __pow2_gt |= __pow2_gt >> 2;
   __pow2_gt |= __pow2_gt >> 4;
   __pow2_gt |= __pow2_gt >> 8;
   __pow2_gt |= __pow2_gt >> 16;
   __pow2_gt |= __pow2_gt >> 32;
   ++__pow2_gt;
   register double const
      normalized_x = x / (double) __pow2_gt;
   register SamplingSixthRootDAccelerator const * const
      accelerator = (SamplingSixthRootDAccelerator const *) 
         __accelerator;
   register double
      step_frac = (normalized_x / accelerator->step_x);
   register uint64_t const
      index = (uint64_t) step_frac;
   step_frac -= index;
   step_frac =
      accelerator->values[index] * (1. - step_frac)
    + accelerator->values[index + (uint64_t) 1] * step_frac;
   register double
      result =
      step_frac * __sixth_roots_of_2n[
         __log2_table[
            (--__pow2_gt * 285673954929480801ul) >> 57
            ]
         ];
#else // assume GCC:
   register unsigned const
      __pow2_exponent = 64u - __builtin_clzl((uint64_t) x);
   SamplingSixthRootDAccelerator const * restrict const
      accelerator = (SamplingSixthRootDAccelerator const *) 
         __accelerator;
   register double
      step_frac = (x * accelerator->step_x_inverse_powers[__pow2_exponent]);
   register uint64_t const
      index = (uint64_t) step_frac;
   step_frac -= index;
   step_frac =
      accelerator->values[index] * (1. - step_frac)
    + accelerator->values[index + (uint64_t) 1] * step_frac;
   register double const
      result = step_frac * __sixth_roots_of_2n[__pow2_exponent];
#endif
   return
      (do_invert ? 1.0/result : result);
}

static
double
sublinear_sampling_normalizing_sixth_rootd_evaluate(
   register double __x,
   register void const * restrict __accelerator
   )
{
   register char const
      do_invert = (__x < 1.0);
   register double const
      x = do_invert ? 1.0/__x : __x;
   register unsigned const
      __pow2_exponent = 64u - __builtin_clzl((uint64_t) x);
   SamplingSixthRootDAccelerator const * restrict const
      accelerator = (SamplingSixthRootDAccelerator const *) 
         __accelerator;
   register double
      step_frac = (x * accelerator->step_x_inverse_powers[__pow2_exponent]);
   register uint64_t const
      index = (uint64_t) step_frac + 1u;
   register double const
      result =
         accelerator->values[index] * __sixth_roots_of_2n[__pow2_exponent];
   return
      (do_invert ? 1.0/result : result);
}

static inline
void
unity_reflection(double * const x)
{
   *x = ((*x >= 1.) ? *x : 1./ *x);
}

static inline
void
conditional_unity_reflection(
   double const * const x,
   double * const y
   )
{
   *y = ((*x >= 1.) ? *y : 1./ *y);
}

static
void
linear_sampling_normalizing_sixth_rootd_batch_evaluate(
   register double const * x,
   register double * out,
   register size_t length,
   register void const * restrict __accelerator
   )
{
   register double const
      * __x = x;
   register double
      * __out = out;
#ifndef __VARR_HAS_AVX__
   for(
      register size_t i = (size_t) 0u;
      i< length;
      ++i
      )
   {
      *__out++ = linear_sampling_normalizing_sixth_rootd_evaluate(
         *__x++,
         __accelerator
         );
   }
#else
   register SamplingSixthRootDAccelerator const * const
      accelerator = ((SamplingSixthRootDAccelerator const *) __accelerator);
   register double const * const
      values = accelerator->values;
   register size_t const
      length_axv_stride = length / __AVX_DOUBLE_STRIDE__;
   double __attribute__((aligned(128)))
      alignment_emulator[__AVX_DOUBLE_STRIDE__];
   // register avxd_array_t
   avxd_array_t
      alpha,
      __values1,
      __values2,
      prefix,
    * const target = (avxd_array_t *) alignment_emulator;
   avxd_array_t const
      __one = _avxd_stride_set_duplicates(1.);
#ifndef __INTEL_COMPILER
   register
#endif
   int64_t
      index_avx[__AVX_DOUBLE_STRIDE__];
   for(register size_t i = (size_t) 0u; i< length_axv_stride; ++i)
   {
      alignment_emulator[0u] = __x[0u];
      alignment_emulator[1u] = __x[1u];
      alignment_emulator[2u] = __x[2u];
      alignment_emulator[3u] = __x[3u];
#ifdef __VARR_USE_AVX512__
      alignment_emulator[4u] = __x[4u];
      alignment_emulator[5u] = __x[5u];
      alignment_emulator[6u] = __x[6u];
      alignment_emulator[7u] = __x[7u];
#endif
      unity_reflection(alignment_emulator + 0u);
      unity_reflection(alignment_emulator + 1u);
      unity_reflection(alignment_emulator + 2u);
      unity_reflection(alignment_emulator + 3u);
#ifdef __VARR_USE_AVX512__
      unity_reflection(alignment_emulator + 4u);
      unity_reflection(alignment_emulator + 5u);
      unity_reflection(alignment_emulator + 6u);
      unity_reflection(alignment_emulator + 7u);
#endif
      index_avx[0u] = 64u - __builtin_clzl((uint64_t) alignment_emulator[0u]);
      index_avx[1u] = 64u - __builtin_clzl((uint64_t) alignment_emulator[1u]);
      index_avx[2u] = 64u - __builtin_clzl((uint64_t) alignment_emulator[2u]);
      index_avx[3u] = 64u - __builtin_clzl((uint64_t) alignment_emulator[3u]);
#ifdef __VARR_USE_AVX512__
      index_avx[4u] = 64u - __builtin_clzl((uint64_t) alignment_emulator[4u]);
      index_avx[5u] = 64u - __builtin_clzl((uint64_t) alignment_emulator[5u]);
      index_avx[6u] = 64u - __builtin_clzl((uint64_t) alignment_emulator[6u]);
      index_avx[7u] = 64u - __builtin_clzl((uint64_t) alignment_emulator[7u]);
#endif
      prefix = _avxd_stride_set(
#ifdef __VARR_USE_AVX512__
         __sixth_roots_of_2n[index_avx[7u]],
         __sixth_roots_of_2n[index_avx[6u]],
         __sixth_roots_of_2n[index_avx[5u]],
         __sixth_roots_of_2n[index_avx[4u]],
         
#endif
         __sixth_roots_of_2n[index_avx[3u]],
         __sixth_roots_of_2n[index_avx[2u]],
         __sixth_roots_of_2n[index_avx[1u]],
         __sixth_roots_of_2n[index_avx[0u]]
         );
      alignment_emulator[0u] *= 
         accelerator->step_x_inverse_powers[index_avx[0u]];
      alignment_emulator[1u] *= 
         accelerator->step_x_inverse_powers[index_avx[1u]];
      alignment_emulator[2u] *= 
         accelerator->step_x_inverse_powers[index_avx[2u]];
      alignment_emulator[3u] *= 
         accelerator->step_x_inverse_powers[index_avx[3u]];
#ifdef __VARR_USE_AVX512__
      alignment_emulator[4u] *= 
         accelerator->step_x_inverse_powers[index_avx[4u]];
      alignment_emulator[5u] *= 
         accelerator->step_x_inverse_powers[index_avx[5u]];
      alignment_emulator[6u] *= 
         accelerator->step_x_inverse_powers[index_avx[6u]];
      alignment_emulator[7u] *= 
         accelerator->step_x_inverse_powers[index_avx[7u]];
#endif
      index_avx[0u] = (int64_t) alignment_emulator[0u];
      index_avx[1u] = (int64_t) alignment_emulator[1u];
      index_avx[2u] = (int64_t) alignment_emulator[2u];
      index_avx[3u] = (int64_t) alignment_emulator[3u];
#ifdef __VARR_USE_AVX512__
      index_avx[4u] = (int64_t) alignment_emulator[4u];
      index_avx[5u] = (int64_t) alignment_emulator[5u];
      index_avx[6u] = (int64_t) alignment_emulator[6u];
      index_avx[7u] = (int64_t) alignment_emulator[7u];
#endif
      alpha = *target -  _avxd_stride_floor(*target);
      __values1 = _avxd_stride_set(
#ifdef __VARR_USE_AVX512__
         values[index_avx[7u]++],
         values[index_avx[6u]++],
         values[index_avx[5u]++],
         values[index_avx[4u]++],
#endif
         values[index_avx[3u]++],
         values[index_avx[2u]++],
         values[index_avx[1u]++],
         values[index_avx[0u]++]
         );
      __values2 = _avxd_stride_set(
#ifdef __VARR_USE_AVX512__
         values[index_avx[7u]],
         values[index_avx[6u]],
         values[index_avx[5u]],
         values[index_avx[4u]],
#endif
         values[index_avx[3u]],
         values[index_avx[2u]],
         values[index_avx[1u]],
         values[index_avx[0u]]
         );
      *target = __values2 * alpha;
      alpha = (__one - alpha);
      *target = (*target + __values1 * alpha);
      *target = prefix * *target;
      conditional_unity_reflection(__x + 0u, alignment_emulator + 0u);
      conditional_unity_reflection(__x + 1u, alignment_emulator + 1u);
      conditional_unity_reflection(__x + 2u, alignment_emulator + 2u);
      conditional_unity_reflection(__x + 3u, alignment_emulator + 3u);
#ifdef __VARR_USE_AVX512__
      conditional_unity_reflection(__x + 4u, alignment_emulator + 4u);
      conditional_unity_reflection(__x + 5u, alignment_emulator + 5u);
      conditional_unity_reflection(__x + 6u, alignment_emulator + 6u);
      conditional_unity_reflection(__x + 7u, alignment_emulator + 7u);
#endif
      __out[0u] = alignment_emulator[0u];
      __out[1u] = alignment_emulator[1u];
      __out[2u] = alignment_emulator[2u];
      __out[3u] = alignment_emulator[3u];
#ifdef __VARR_USE_AVX512__
      __out[4u] = alignment_emulator[4u];
      __out[5u] = alignment_emulator[5u];
      __out[6u] = alignment_emulator[6u];
      __out[7u] = alignment_emulator[7u];
#endif
      __x += __AVX_DOUBLE_STRIDE__;
      __out += __AVX_DOUBLE_STRIDE__;
      continue;
   }
   
   // Remainder loop:
   
   for(
      register size_t i = length_axv_stride * __AVX_DOUBLE_STRIDE__;
      i< length; 
      ++i
      )
   {
      *__out++ = linear_sampling_normalizing_sixth_rootd_evaluate(
         *__x++,
         __accelerator
         );
   }
#endif
   return;
}


static
void
linear_sampling_normalizing_sixth_rootd_batch_burst_evaluate(
   register double const * x,
   register double * out,
   register size_t length,
   register void const * restrict __accelerator
   )
{
   register double const
      * __x = x;
   register double
      * __out = out;
#ifndef __VARR_HAS_AVX__
   for(
      register size_t i = (size_t) 0u;
      i< length;
      ++i
      )
   {
      *__out++ = linear_sampling_normalizing_sixth_rootd_evaluate(
         *__x++,
         __accelerator
         );
   }
#else
   register SamplingSixthRootDAccelerator const * const
      accelerator = ((SamplingSixthRootDAccelerator const *) __accelerator);
   register double const * const
      values = accelerator->values;
   register size_t const
      length_axv_stride = length / __AVX_DOUBLE_STRIDE__;
   static size_t const
      maximum_batch_size = 8192u;
   avxd_array_t const
      __one = _avxd_stride_set_duplicates(1.);
   // register avxd_array_t
   avxd_array_t
      alpha,
      __values1,
      __values2,
    * prefix,
    * target;
   size_t const
      storage_bins = __AVX_DOUBLE_STRIDE__ * maximum_batch_size;
   int64_t
      __attribute__((aligned(128)))
      __int64_workspace[storage_bins];
   double 
      __attribute__((aligned(128)))
      __double_workspace[storage_bins];
   double
      __attribute__((aligned(128)))
      __double_workspace2[storage_bins];
   int64_t *
      int64_workspace = __int64_workspace;
   double
      * double_workspace = __double_workspace,
      * double_workspace2 = __double_workspace2;
   size_t
      samples_remaining = length_axv_stride,
      batch_size =
         (samples_remaining >= maximum_batch_size) ? 
            maximum_batch_size
          : samples_remaining;
   samples_remaining -= batch_size;
   register size_t
      i;
   
   while(1) {
   
   __x = x;
   
   double_workspace = __double_workspace;
   
   memcpy(
      double_workspace,
      __x,
      sizeof(double) * batch_size * __AVX_DOUBLE_STRIDE__
      );
   
   for(
      i = (size_t) 0u;
      i< batch_size;
      ++i,
      (double_workspace += __AVX_DOUBLE_STRIDE__)
      )
   {
      unity_reflection(double_workspace + 0u);
      unity_reflection(double_workspace + 1u);
      unity_reflection(double_workspace + 2u);
      unity_reflection(double_workspace + 3u);
#ifdef __VARR_USE_AVX512__
      unity_reflection(double_workspace + 4u);
      unity_reflection(double_workspace + 5u);
      unity_reflection(double_workspace + 6u);
      unity_reflection(double_workspace + 7u);
#endif
      
      continue;
   }
   
   __x += __AVX_DOUBLE_STRIDE__ * batch_size;
   int64_workspace = __int64_workspace;
   double_workspace = __double_workspace;
   
   for(
      i = (size_t) 0u;
      i< batch_size;
      ++i,
      (int64_workspace += __AVX_DOUBLE_STRIDE__),
      (double_workspace += __AVX_DOUBLE_STRIDE__)
      )
   {
      int64_workspace[0u] =
         64u - __builtin_clzl((uint64_t) double_workspace[0u]);
      int64_workspace[1u] =
         64u - __builtin_clzl((uint64_t) double_workspace[1u]);
      int64_workspace[2u] =
         64u - __builtin_clzl((uint64_t) double_workspace[2u]);
      int64_workspace[3u] =
         64u - __builtin_clzl((uint64_t) double_workspace[3u]);
#ifdef __VARR_USE_AVX512__
      int64_workspace[4u] =
         64u - __builtin_clzl((uint64_t) double_workspace[4u]);
      int64_workspace[5u] =
         64u - __builtin_clzl((uint64_t) double_workspace[5u]);
      int64_workspace[6u] =
         64u - __builtin_clzl((uint64_t) double_workspace[6u]);
      int64_workspace[7u] =
         64u - __builtin_clzl((uint64_t) double_workspace[7u]);
#endif
      
      continue;
   }
   
   int64_workspace = __int64_workspace;
   prefix = (avxd_array_t *) __double_workspace2;
   
   for(
      i = (size_t) 0u;
      i< batch_size;
      ++i,
      (int64_workspace += __AVX_DOUBLE_STRIDE__),
      ++prefix
      )
   {
      *prefix = _avxd_stride_set(
#ifdef __VARR_USE_AVX512__
         __sixth_roots_of_2n[int64_workspace[7u]],
         __sixth_roots_of_2n[int64_workspace[6u]],
         __sixth_roots_of_2n[int64_workspace[5u]],
         __sixth_roots_of_2n[int64_workspace[4u]],
         
#endif
         __sixth_roots_of_2n[int64_workspace[3u]],
         __sixth_roots_of_2n[int64_workspace[2u]],
         __sixth_roots_of_2n[int64_workspace[1u]],
         __sixth_roots_of_2n[int64_workspace[0u]]
         );
      
      continue;
   }
   
   double_workspace = __double_workspace;
   int64_workspace = __int64_workspace;
   
   for(
      i = (size_t) 0u;
      i< batch_size;
      ++i,
      (double_workspace += __AVX_DOUBLE_STRIDE__),
      (int64_workspace += __AVX_DOUBLE_STRIDE__)
      )
   {
      double_workspace[0u] *= 
         accelerator->step_x_inverse_powers[int64_workspace[0u]];
      double_workspace[1u] *= 
         accelerator->step_x_inverse_powers[int64_workspace[1u]];
      double_workspace[2u] *= 
         accelerator->step_x_inverse_powers[int64_workspace[2u]];
      double_workspace[3u] *= 
         accelerator->step_x_inverse_powers[int64_workspace[3u]];
#ifdef __VARR_USE_AVX512__
      double_workspace[4u] *= 
         accelerator->step_x_inverse_powers[int64_workspace[4u]];
      double_workspace[5u] *= 
         accelerator->step_x_inverse_powers[int64_workspace[5u]];
      double_workspace[6u] *= 
         accelerator->step_x_inverse_powers[int64_workspace[6u]];
      double_workspace[7u] *= 
         accelerator->step_x_inverse_powers[int64_workspace[7u]];
#endif
      
      continue;
   }
   
   double_workspace = __double_workspace;
   int64_workspace = __int64_workspace;
   
   for(
      i = (size_t) 0u;
      i< batch_size;
      ++i, 
      (double_workspace += __AVX_DOUBLE_STRIDE__),
      (int64_workspace += __AVX_DOUBLE_STRIDE__)
      )
   {
      int64_workspace[0u] = (int64_t) double_workspace[0u];
      int64_workspace[1u] = (int64_t) double_workspace[1u];
      int64_workspace[2u] = (int64_t) double_workspace[2u];
      int64_workspace[3u] = (int64_t) double_workspace[3u];
#ifdef __VARR_USE_AVX512__
      int64_workspace[4u] = (int64_t) double_workspace[4u];
      int64_workspace[5u] = (int64_t) double_workspace[5u];
      int64_workspace[6u] = (int64_t) double_workspace[6u];
      int64_workspace[7u] = (int64_t) double_workspace[7u];
#endif
      
      continue;
   }
   
   double_workspace = __double_workspace;
   double_workspace2 = __double_workspace2;
   int64_workspace = __int64_workspace;
   
   target = (avxd_array_t *) double_workspace;
   prefix = (avxd_array_t *) double_workspace2;
   
   __x = x;
   
   for(
      i = (size_t) 0u;
      i< batch_size;
      ++i,
      (double_workspace += __AVX_DOUBLE_STRIDE__),
      (int64_workspace += __AVX_DOUBLE_STRIDE__),
      ++prefix,
      ++target
      )
   {
      alpha = *target -  _avxd_stride_floor(*target);
      __values1 = _avxd_stride_set(
#ifdef __VARR_USE_AVX512__
         values[int64_workspace[7u]++],
         values[int64_workspace[6u]++],
         values[int64_workspace[5u]++],
         values[int64_workspace[4u]++],
#endif
         values[int64_workspace[3u]++],
         values[int64_workspace[2u]++],
         values[int64_workspace[1u]++],
         values[int64_workspace[0u]++]
         );
      __values2 = _avxd_stride_set(
#ifdef __VARR_USE_AVX512__
         values[int64_workspace[7u]],
         values[int64_workspace[6u]],
         values[int64_workspace[5u]],
         values[int64_workspace[4u]],
#endif
         values[int64_workspace[3u]],
         values[int64_workspace[2u]],
         values[int64_workspace[1u]],
         values[int64_workspace[0u]]
         );
      *target = __values2 * alpha;
      alpha = (__one - alpha);
      *target = (*target + __values1 * alpha);
      *target = *prefix * *target;
      conditional_unity_reflection(__x + 0u, double_workspace + 0u);
      conditional_unity_reflection(__x + 1u, double_workspace + 1u);
      conditional_unity_reflection(__x + 2u, double_workspace + 2u);
      conditional_unity_reflection(__x + 3u, double_workspace + 3u);
#ifdef __VARR_USE_AVX512__
      conditional_unity_reflection(__x + 4u, double_workspace + 4u);
      conditional_unity_reflection(__x + 5u, double_workspace + 5u);
      conditional_unity_reflection(__x + 6u, double_workspace + 6u);
      conditional_unity_reflection(__x + 7u, double_workspace + 7u);
#endif
      __x += __AVX_DOUBLE_STRIDE__;
      
      continue;
   }
   
   memcpy(
      __out,
      __double_workspace,
      sizeof(double) * batch_size * __AVX_DOUBLE_STRIDE__
      );
   __out += batch_size * __AVX_DOUBLE_STRIDE__;
   
   x += batch_size * __AVX_DOUBLE_STRIDE__;
   
   if(samples_remaining == 0u)
   {
      break;
   }
   
   if(samples_remaining >= batch_size)
      samples_remaining -= batch_size;
   else
   {
      batch_size = samples_remaining;
      samples_remaining = 0u;
   }
   
   continue;
   
   }
   
   __x = x;
   
   // Remainder loop:
   
   for(
      i = length_axv_stride * __AVX_DOUBLE_STRIDE__;
      i< length;
      ++i
      )
   {
      *__out++ = linear_sampling_normalizing_sixth_rootd_evaluate(
         *__x++,
         __accelerator
         );
   }
#endif
   return;
}

static
void
sublinear_sampling_normalizing_sixth_rootd_batch_evaluate(
   register double const * x,
   register double * out,
   register size_t length,
   register void const * restrict __accelerator
   )
{
   register double const
      * __x = x;
   register double
      * __out = out;
#ifndef __VARR_HAS_AVX__
   for(
      register size_t i = (size_t) 0u;
      i< length;
      ++i
      )
      *__out++ = sublinear_sampling_normalizing_sixth_rootd_evaluate(
         *__x++,
         __accelerator
         );
#else
   register SamplingSixthRootDAccelerator const *
      accelerator = ((SamplingSixthRootDAccelerator const *) __accelerator);
   register double const * const
      values = accelerator->values;
   register size_t const
      length_axv_stride = length / __AVX_DOUBLE_STRIDE__;
   double __attribute__((aligned(128)))
      alignment_emulator[__AVX_DOUBLE_STRIDE__];
   // register avxd_array_t
   avxd_array_t
      __value,
      prefix,
    * const target = (avxd_array_t *) alignment_emulator;
#ifndef __INTEL_COMPILER
   register
#endif
   int64_t
      index_avx[__AVX_DOUBLE_STRIDE__];
   for(register size_t i = (size_t) 0u; i< length_axv_stride; ++i)
   {
      alignment_emulator[0u] = __x[0u];
      alignment_emulator[1u] = __x[1u];
      alignment_emulator[2u] = __x[2u];
      alignment_emulator[3u] = __x[3u];
#ifdef __VARR_USE_AVX512__
      alignment_emulator[4u] = __x[4u];
      alignment_emulator[5u] = __x[5u];
      alignment_emulator[6u] = __x[6u];
      alignment_emulator[7u] = __x[7u];
#endif
      unity_reflection(alignment_emulator + 0u);
      unity_reflection(alignment_emulator + 1u);
      unity_reflection(alignment_emulator + 2u);
      unity_reflection(alignment_emulator + 3u);
#ifdef __VARR_USE_AVX512__
      unity_reflection(alignment_emulator + 4u);
      unity_reflection(alignment_emulator + 5u);
      unity_reflection(alignment_emulator + 6u);
      unity_reflection(alignment_emulator + 7u);
#endif
      index_avx[0u] = 64u - __builtin_clzl((uint64_t) alignment_emulator[0u]);
      index_avx[1u] = 64u - __builtin_clzl((uint64_t) alignment_emulator[1u]);
      index_avx[2u] = 64u - __builtin_clzl((uint64_t) alignment_emulator[2u]);
      index_avx[3u] = 64u - __builtin_clzl((uint64_t) alignment_emulator[3u]);
#ifdef __VARR_USE_AVX512__
      index_avx[4u] = 64u - __builtin_clzl((uint64_t) alignment_emulator[4u]);
      index_avx[5u] = 64u - __builtin_clzl((uint64_t) alignment_emulator[5u]);
      index_avx[6u] = 64u - __builtin_clzl((uint64_t) alignment_emulator[6u]);
      index_avx[7u] = 64u - __builtin_clzl((uint64_t) alignment_emulator[7u]);
#endif
      prefix = _avxd_stride_set(
#ifdef __VARR_USE_AVX512__
         __sixth_roots_of_2n[index_avx[7u]],
         __sixth_roots_of_2n[index_avx[6u]],
         __sixth_roots_of_2n[index_avx[5u]],
         __sixth_roots_of_2n[index_avx[4u]],
         
#endif
         __sixth_roots_of_2n[index_avx[3u]],
         __sixth_roots_of_2n[index_avx[2u]],
         __sixth_roots_of_2n[index_avx[1u]],
         __sixth_roots_of_2n[index_avx[0u]]
         );
      alignment_emulator[0u] *= 
         accelerator->step_x_inverse_powers[index_avx[0u]];
      alignment_emulator[1u] *= 
         accelerator->step_x_inverse_powers[index_avx[1u]];
      alignment_emulator[2u] *= 
         accelerator->step_x_inverse_powers[index_avx[2u]];
      alignment_emulator[3u] *= 
         accelerator->step_x_inverse_powers[index_avx[3u]];
#ifdef __VARR_USE_AVX512__
      alignment_emulator[4u] *= 
         accelerator->step_x_inverse_powers[index_avx[4u]];
      alignment_emulator[5u] *= 
         accelerator->step_x_inverse_powers[index_avx[5u]];
      alignment_emulator[6u] *= 
         accelerator->step_x_inverse_powers[index_avx[6u]];
      alignment_emulator[7u] *= 
         accelerator->step_x_inverse_powers[index_avx[7u]];
#endif
      index_avx[0u] = (int64_t) alignment_emulator[0u] + 1u;
      index_avx[1u] = (int64_t) alignment_emulator[1u] + 1u;
      index_avx[2u] = (int64_t) alignment_emulator[2u] + 1u;
      index_avx[3u] = (int64_t) alignment_emulator[3u] + 1u;
#ifdef __VARR_USE_AVX512__
      index_avx[4u] = (int64_t) alignment_emulator[4u] + 1u;
      index_avx[5u] = (int64_t) alignment_emulator[5u] + 1u;
      index_avx[6u] = (int64_t) alignment_emulator[6u] + 1u;
      index_avx[7u] = (int64_t) alignment_emulator[7u] + 1u;
#endif
      __value = _avxd_stride_set(
#ifdef __VARR_USE_AVX512__
         values[index_avx[7u]],
         values[index_avx[6u]],
         values[index_avx[5u]],
         values[index_avx[4u]],
#endif
         values[index_avx[3u]],
         values[index_avx[2u]],
         values[index_avx[1u]],
         values[index_avx[0u]]
         );
      *target = prefix * __value;
      conditional_unity_reflection(__x + 0u, alignment_emulator + 0u);
      conditional_unity_reflection(__x + 1u, alignment_emulator + 1u);
      conditional_unity_reflection(__x + 2u, alignment_emulator + 2u);
      conditional_unity_reflection(__x + 3u, alignment_emulator + 3u);
#ifdef __VARR_USE_AVX512__
      conditional_unity_reflection(__x + 4u, alignment_emulator + 4u);
      conditional_unity_reflection(__x + 5u, alignment_emulator + 5u);
      conditional_unity_reflection(__x + 6u, alignment_emulator + 6u);
      conditional_unity_reflection(__x + 7u, alignment_emulator + 7u);
#endif
      __out[0u] = alignment_emulator[0u];
      __out[1u] = alignment_emulator[1u];
      __out[2u] = alignment_emulator[2u];
      __out[3u] = alignment_emulator[3u];
#ifdef __VARR_USE_AVX512__
      __out[4u] = alignment_emulator[4u];
      __out[5u] = alignment_emulator[5u];
      __out[6u] = alignment_emulator[6u];
      __out[7u] = alignment_emulator[7u];
#endif
      __x += __AVX_DOUBLE_STRIDE__;
      __out += __AVX_DOUBLE_STRIDE__;
      continue;
   }
   
   // Remainder loop:
   
   for(
      register size_t i = length_axv_stride * __AVX_DOUBLE_STRIDE__;
      i< length; 
      ++i
      )
   {
      *__out++ = sublinear_sampling_normalizing_sixth_rootd_evaluate(
         *__x++,
         __accelerator
         );
   }
#endif
   return;
}

VARRSixthRootDEvaluator
linear_sampling_normalizing_sixth_rootd(
   size_t number_of_samples
   )
{
   (void) linear_sampling_normalizing_sixth_rootd_batch_evaluate;
   (void) linear_sampling_normalizing_sixth_rootd_batch_burst_evaluate;
   
   VARRSixthRootDEvaluator
      result;
   SamplingSixthRootDAccelerator * const
      accelerator = allocate_linear_sampling_normalizing_sixth_rootd(
         number_of_samples
         );
   result.accelerator = (void *) accelerator;
   result.sixthrootd = linear_sampling_normalizing_sixth_rootd_evaluate;
   result.sixthrootd_array = 
      linear_sampling_normalizing_sixth_rootd_batch_evaluate;
   result.disallocate = linear_sampling_normalizing_sixth_rootd_disallocate;
   return
      result;
}

VARRSixthRootDEvaluator
sublinear_sampling_normalizing_sixth_rootd(
   size_t number_of_samples
   )
{
   VARRSixthRootDEvaluator
      result;
   SamplingSixthRootDAccelerator * const
      accelerator = allocate_linear_sampling_normalizing_sixth_rootd(
         number_of_samples
         );
   result.accelerator = (void *) accelerator;
   result.sixthrootd = sublinear_sampling_normalizing_sixth_rootd_evaluate;
   result.sixthrootd_array = 
      sublinear_sampling_normalizing_sixth_rootd_batch_evaluate;
   result.disallocate = linear_sampling_normalizing_sixth_rootd_disallocate;
   return
      result;
}
