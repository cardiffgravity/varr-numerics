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

#include "varr_log.h"
#include "varr_internal.h"

#include <math.h>
#include <stdlib.h>
#include <inttypes.h>

typedef struct tagSamplingLogDAccelerator {
   double
      start_x,                                        // Interval [1, 2]
      step_x,
      step_x_inv;
   size_t
      samples;
   double const *
      values;
} SamplingLogDAccelerator;

static
int
__linear_sampling_normalizing_logd_disallocate(
   SamplingLogDAccelerator * accelerator
   )
{
   if(accelerator == NULL)
   {
      return 1;
   }
   free((void *) accelerator->values);
   accelerator->start_x = 0x7F800001;
   accelerator->step_x = 0x7F800001;
   accelerator->step_x_inv = 0x7F800001;
   accelerator->samples = (size_t) 0;
   accelerator->values = NULL;
   free(accelerator);
   return 0;
}

static
int
linear_sampling_normalizing_logd_disallocate(
   VARRLogDEvaluator * evaluator
   )
{
   if(evaluator == NULL)
   {
      return 1;
   }
   SamplingLogDAccelerator * const
      accelerator = (SamplingLogDAccelerator *)
         evaluator->accelerator;
   return
      __linear_sampling_normalizing_logd_disallocate(accelerator);
}

static
SamplingLogDAccelerator *
allocate_linear_sampling_normalizing_logd(
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
         x = 1.0 + step_size * (double) i,
         value = log(x);
      values[i] = value;
      continue;
   }
   SamplingLogDAccelerator * const
      result = (SamplingLogDAccelerator *) malloc(
         sizeof(SamplingLogDAccelerator)
         );
   result->start_x = 1.0;
   result->step_x = step_size;
   result->step_x_inv = 1./step_size;
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
   __log_2k[] = {
                      0.0,
   0.69314718055994528623,
    1.3862943611198905725,
    2.0794415416798357477,
    2.7725887222397811449,
    3.4657359027997265422,
    4.1588830833596714953,
    4.8520302639196168926,
    5.5451774444795622898,
    6.2383246250395076871,
    6.9314718055994530843,
    7.6246189861593984816,
    8.3177661667193429906,
    9.0109133472792883879,
    9.7040605278392337851,
    10.397207708399179182,
     11.09035488895912458,
    11.783502069519069977,
    12.476649250079015374,
    13.169796430638960771,
    13.862943611198906169,
    14.556090791758851566,
    15.249237972318796963,
     15.94238515287874236,
    16.635532333438685981,
    17.328679513998633155,
    18.021826694558576776,
    18.714973875118523949,
     19.40812105567846757,
    20.101268236238414744,
    20.794415416798358365,
    21.487562597358305538,
    22.180709777918249159,
    22.873856958478196333,
    23.567004139038139954,
    24.260151319598087127,
    24.953298500158030748,
    25.646445680717977922,
    26.339592861277921543,
    27.032740041837868716,
    27.725887222397812337,
    28.419034402957755958,
    29.112181583517703132,
    29.805328764077646753,
    30.498475944637593926,
    31.191623125197537547,
    31.884770305757484721,
    32.577917486317431894,
    33.271064666877371963,
    33.964211847437319136,
     34.65735902799726631,
    35.350506208557213483,
    36.043653389117153552,
    36.736800569677100725,
    37.429947750237047899,
    38.123094930796995072,
    38.816242111356935141,
    39.509389291916882314,
    40.202536472476829488,
    40.895683653036776661,
     41.58883083359671673,
    42.281978014156663903,
    42.975125194716611077,
    43.668272375276551145
   };

#include "varr_floor_log2.h"

static inline
void
unity_reflection(double * const x)
{
   *x = ((*x >= 1.) ? *x : 1./ *x);
}

static inline
void
conditional_sign_mirror(
   double const * const x,
   double * const y
   )
{
   *y = ((*x >= 1.) ? *y : -*y);
}

static
double
linear_sampling_normalizing_logd_evaluate(
   register double __x,
   register void const * restrict __accelerator
   )
{
   register char const
      do_invert = (__x < 1.0);
   register double const 
      x = (do_invert ? (1.0/__x) : __x);
#ifdef __INTEL_COMPILER
   /*
    * An exponential split log(x) = log(x/2**k) + log(2**k),
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
   __pow2_gt >>= 1;
   ++__pow2_gt;
   register double const
      normalized_x = (x / (double) __pow2_gt) - 1.0;
   register SamplingLogDAccelerator const * const
      accelerator = (SamplingLogDAccelerator const *) 
         __accelerator;
   register double
      step_frac = (normalized_x / accelerator->step_x);
   register uint64_t const
      index = (uint64_t) step_frac;
   step_frac -= index;
   step_frac =
      accelerator->values[index] * (1. - step_frac)
    + accelerator->values[index + (uint64_t) 1] * step_frac;
   return
      (1 - 2 * do_invert) *
      (
      step_frac
    + __log_2k[
         __log2_table[
            ((__pow2_gt - (uint64_t) 1) * 285673954929480801ul) >> 57
            ]
         ]
      );
#else // assume GCC:
   register unsigned const
      __pow2_exponent = 63u - __builtin_clzl((uint64_t) x);
   register uint64_t
      __pow2_gt = ((uint64_t) 1) << __pow2_exponent;
   register double const
      normalized_x = (x / (double) __pow2_gt) - 1.0;
   register double
      step_frac = (
         normalized_x / 
            ((SamplingLogDAccelerator const *) __accelerator)->step_x
         );
   register uint64_t const
      index = (uint64_t) step_frac;
   step_frac -= index;
   register double const *
      value = ((SamplingLogDAccelerator const *) __accelerator)->values + index;
   step_frac =
      (*value) * (1. - step_frac) + value[1] * step_frac;
   return
      (1 - 2 * do_invert)
    * (step_frac + __log_2k[__pow2_exponent]);
#endif
}

static
double
sublinear_sampling_normalizing_logd_evaluate(
   register double __x,
   register void const * restrict __accelerator
   )
{
   register char const
      do_invert = (__x < 1.0);
   register double const 
      x = (do_invert ? (1.0/__x) : __x);
   register unsigned const
      __pow2_exponent = 63u - __builtin_clzl((uint64_t) x);
   register uint64_t
      __pow2_gt = ((uint64_t) 1) << __pow2_exponent;
   register double const
      normalized_x = (x / (double) __pow2_gt) - 1.0;
   register double
      step_frac = (
         normalized_x / 
            ((SamplingLogDAccelerator const *) __accelerator)->step_x
         );
   register uint64_t const
      index = (uint64_t) step_frac + 1u;
   register double const *
      value = ((SamplingLogDAccelerator const *) __accelerator)->values + index;
   step_frac = *value;
   return
      (1 - 2 * do_invert)
    * (step_frac + __log_2k[__pow2_exponent]);
}

static
void
linear_sampling_normalizing_logd_batch_evaluate(
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
      *__out++ = linear_sampling_normalizing_logd_evaluate(
         *__x++,
         __accelerator
         );
   }
#else
   register SamplingLogDAccelerator const * const
      accelerator = ((SamplingLogDAccelerator const *) __accelerator);
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
      modifier,
    * const target = (avxd_array_t *) alignment_emulator;
   avxd_array_t const
      __one = _avxd_stride_set_duplicates(1.),
      step_x_inv = _avxd_stride_set_duplicates(accelerator->step_x_inv);
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
      index_avx[0u] = 63u - __builtin_clzl((uint64_t) alignment_emulator[0u]);
      index_avx[1u] = 63u - __builtin_clzl((uint64_t) alignment_emulator[1u]);
      index_avx[2u] = 63u - __builtin_clzl((uint64_t) alignment_emulator[2u]);
      index_avx[3u] = 63u - __builtin_clzl((uint64_t) alignment_emulator[3u]);
#ifdef __VARR_USE_AVX512__
      index_avx[4u] = 63u - __builtin_clzl((uint64_t) alignment_emulator[4u]);
      index_avx[5u] = 63u - __builtin_clzl((uint64_t) alignment_emulator[5u]);
      index_avx[6u] = 63u - __builtin_clzl((uint64_t) alignment_emulator[6u]);
      index_avx[7u] = 63u - __builtin_clzl((uint64_t) alignment_emulator[7u]);
#endif
      prefix = _avxd_stride_set(
#ifdef __VARR_USE_AVX512__
         __log_2k[index_avx[7u]],
         __log_2k[index_avx[6u]],
         __log_2k[index_avx[5u]],
         __log_2k[index_avx[4u]],
         
#endif
         __log_2k[index_avx[3u]],
         __log_2k[index_avx[2u]],
         __log_2k[index_avx[1u]],
         __log_2k[index_avx[0u]]
         );
      modifier = _avxd_stride_set(
#ifdef __VARR_USE_AVX512__
         ((uint64_t) 1) << index_avx[7u],
         ((uint64_t) 1) << index_avx[6u],
         ((uint64_t) 1) << index_avx[5u],
         ((uint64_t) 1) << index_avx[4u],
#endif
         ((uint64_t) 1) << index_avx[3u],
         ((uint64_t) 1) << index_avx[2u],
         ((uint64_t) 1) << index_avx[1u],
         ((uint64_t) 1) << index_avx[0u]
         );
      *target = (*target / modifier - __one) * step_x_inv;
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
      *target = prefix + *target;
      conditional_sign_mirror(__x + 0u, alignment_emulator + 0u);
      conditional_sign_mirror(__x + 1u, alignment_emulator + 1u);
      conditional_sign_mirror(__x + 2u, alignment_emulator + 2u);
      conditional_sign_mirror(__x + 3u, alignment_emulator + 3u);
#ifdef __VARR_USE_AVX512__
      conditional_sign_mirror(__x + 4u, alignment_emulator + 4u);
      conditional_sign_mirror(__x + 5u, alignment_emulator + 5u);
      conditional_sign_mirror(__x + 6u, alignment_emulator + 6u);
      conditional_sign_mirror(__x + 7u, alignment_emulator + 7u);
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
      *__out++ = linear_sampling_normalizing_logd_evaluate(
         *__x++,
         __accelerator
         );
   }
#endif
   return;
}

static
void
sublinear_sampling_normalizing_logd_batch_evaluate(
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
      *__out++ = sublinear_sampling_normalizing_logd_evaluate(
         *__x++,
         __accelerator
         );
   }
#else
   register SamplingLogDAccelerator const * const
      accelerator = ((SamplingLogDAccelerator const *) __accelerator);
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
      modifier,
    * const target = (avxd_array_t *) alignment_emulator;
   avxd_array_t const
      __one = _avxd_stride_set_duplicates(1.),
      step_x_inv = _avxd_stride_set_duplicates(accelerator->step_x_inv);
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
      index_avx[0u] = 63u - __builtin_clzl((uint64_t) alignment_emulator[0u]);
      index_avx[1u] = 63u - __builtin_clzl((uint64_t) alignment_emulator[1u]);
      index_avx[2u] = 63u - __builtin_clzl((uint64_t) alignment_emulator[2u]);
      index_avx[3u] = 63u - __builtin_clzl((uint64_t) alignment_emulator[3u]);
#ifdef __VARR_USE_AVX512__
      index_avx[4u] = 63u - __builtin_clzl((uint64_t) alignment_emulator[4u]);
      index_avx[5u] = 63u - __builtin_clzl((uint64_t) alignment_emulator[5u]);
      index_avx[6u] = 63u - __builtin_clzl((uint64_t) alignment_emulator[6u]);
      index_avx[7u] = 63u - __builtin_clzl((uint64_t) alignment_emulator[7u]);
#endif
      prefix = _avxd_stride_set(
#ifdef __VARR_USE_AVX512__
         __log_2k[index_avx[7u]],
         __log_2k[index_avx[6u]],
         __log_2k[index_avx[5u]],
         __log_2k[index_avx[4u]],
         
#endif
         __log_2k[index_avx[3u]],
         __log_2k[index_avx[2u]],
         __log_2k[index_avx[1u]],
         __log_2k[index_avx[0u]]
         );
      modifier = _avxd_stride_set(
#ifdef __VARR_USE_AVX512__
         ((uint64_t) 1) << index_avx[7u],
         ((uint64_t) 1) << index_avx[6u],
         ((uint64_t) 1) << index_avx[5u],
         ((uint64_t) 1) << index_avx[4u],
#endif
         ((uint64_t) 1) << index_avx[3u],
         ((uint64_t) 1) << index_avx[2u],
         ((uint64_t) 1) << index_avx[1u],
         ((uint64_t) 1) << index_avx[0u]
         );
      *target = (*target / modifier - __one) * step_x_inv;
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
      *target = prefix + __value;
      conditional_sign_mirror(__x + 0u, alignment_emulator + 0u);
      conditional_sign_mirror(__x + 1u, alignment_emulator + 1u);
      conditional_sign_mirror(__x + 2u, alignment_emulator + 2u);
      conditional_sign_mirror(__x + 3u, alignment_emulator + 3u);
#ifdef __VARR_USE_AVX512__
      conditional_sign_mirror(__x + 4u, alignment_emulator + 4u);
      conditional_sign_mirror(__x + 5u, alignment_emulator + 5u);
      conditional_sign_mirror(__x + 6u, alignment_emulator + 6u);
      conditional_sign_mirror(__x + 7u, alignment_emulator + 7u);
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
      *__out++ = sublinear_sampling_normalizing_logd_evaluate(
         *__x++,
         __accelerator
         );
   }
#endif
   return;
}

VARRLogDEvaluator
normalizing_linear_sampling_logd(
   size_t number_of_samples
   )
{
   VARRLogDEvaluator
      result;
   SamplingLogDAccelerator * const
      accelerator = allocate_linear_sampling_normalizing_logd(
         number_of_samples
         );
   result.accelerator = (void *) accelerator;
   result.logd = linear_sampling_normalizing_logd_evaluate;
   result.logd_array = linear_sampling_normalizing_logd_batch_evaluate;
   result.disallocate = linear_sampling_normalizing_logd_disallocate;
   return
      result;
}

VARRLogDEvaluator
normalizing_sublinear_sampling_logd(
   size_t number_of_samples
   )
{
   VARRLogDEvaluator
      result;
   SamplingLogDAccelerator * const
      accelerator = allocate_linear_sampling_normalizing_logd(
         number_of_samples
         );
   result.accelerator = (void *) accelerator;
   result.logd = sublinear_sampling_normalizing_logd_evaluate;
   result.logd_array = sublinear_sampling_normalizing_logd_batch_evaluate;
   result.disallocate = linear_sampling_normalizing_logd_disallocate;
   return
      result;
}

static
int
quad_series_logd_disallocate(
   VARRLogDEvaluator * evaluator
   )
{
   return 1;
}

static double const
   __constants[] = {
      2.0/1.0, 
      2.0/3.0, 
      2.0/5.0, 
      2.0/7.0, 
      2.0/9.0, 
      2.0/11.0, 
      2.0/13.0, 
      2.0/15.0, 
      2.0/17.0, 
      2.0/19.0, 
      2.0/21.0, 
      2.0/23.0, 
      2.0/25.0, 
      2.0/27.0, 
      2.0/29.0, 
      2.0/31.0, 
      2.0/33.0, 
      2.0/35.0, 
      2.0/37.0, 
      2.0/39.0, 
      2.0/41.0, 
      2.0/43.0, 
      2.0/45.0, 
      2.0/47.0, 
      2.0/49.0, 
      2.0/51.0, 
      2.0/53.0, 
      2.0/55.0, 
      2.0/57.0, 
      2.0/59.0, 
      2.0/61.0, 
      2.0/63.0, 
      2.0/65.0, 
      2.0/67.0, 
      2.0/69.0, 
      2.0/71.0, 
      2.0/73.0, 
      2.0/75.0, 
      2.0/77.0, 
      2.0/79.0, 
      2.0/81.0, 
      2.0/83.0, 
      2.0/85.0, 
      2.0/87.0, 
      2.0/89.0, 
      2.0/91.0, 
      2.0/93.0, 
      2.0/95.0, 
      2.0/97.0, 
      2.0/99.0, 
      2.0/101.0, 
      2.0/103.0, 
      2.0/105.0, 
      2.0/107.0, 
      2.0/109.0, 
      2.0/111.0, 
      2.0/113.0, 
      2.0/115.0, 
      2.0/117.0, 
      2.0/119.0, 
      2.0/121.0, 
      2.0/123.0, 
      2.0/125.0, 
      2.0/127.0, 
      2.0/129.0, 
      2.0/131.0, 
      2.0/133.0, 
      2.0/135.0, 
      2.0/137.0, 
      2.0/139.0, 
      2.0/141.0, 
      2.0/143.0, 
      2.0/145.0, 
      2.0/147.0, 
      2.0/149.0, 
      2.0/151.0, 
      2.0/153.0, 
      2.0/155.0, 
      2.0/157.0, 
      2.0/159.0, 
      2.0/161.0, 
      2.0/163.0, 
      2.0/165.0, 
      2.0/167.0, 
      2.0/169.0, 
      2.0/171.0, 
      2.0/173.0, 
      2.0/175.0, 
      2.0/177.0, 
      2.0/179.0, 
      2.0/181.0, 
      2.0/183.0, 
      2.0/185.0, 
      2.0/187.0, 
      2.0/189.0, 
      2.0/191.0, 
      2.0/193.0, 
      2.0/195.0, 
      2.0/197.0, 
      2.0/199.0, 
      2.0/201.0, 
      2.0/203.0, 
      2.0/205.0, 
      2.0/207.0, 
      2.0/209.0, 
      2.0/211.0, 
      2.0/213.0, 
      2.0/215.0, 
      2.0/217.0, 
      2.0/219.0, 
      2.0/221.0, 
      2.0/223.0, 
      2.0/225.0, 
      2.0/227.0, 
      2.0/229.0, 
      2.0/231.0, 
      2.0/233.0, 
      2.0/235.0, 
      2.0/237.0, 
      2.0/239.0, 
      2.0/241.0, 
      2.0/243.0, 
      2.0/245.0, 
      2.0/247.0, 
      2.0/249.0, 
      2.0/251.0, 
      2.0/253.0, 
      2.0/255.0, 
      2.0/257.0, 
      2.0/259.0, 
      2.0/261.0, 
      2.0/263.0, 
      2.0/265.0, 
      2.0/267.0, 
      2.0/269.0, 
      2.0/271.0, 
      2.0/273.0, 
      2.0/275.0, 
      2.0/277.0, 
      2.0/279.0, 
      2.0/281.0, 
      2.0/283.0, 
      2.0/285.0, 
      2.0/287.0, 
      2.0/289.0, 
      2.0/291.0, 
      2.0/293.0, 
      2.0/295.0, 
      2.0/297.0, 
      2.0/299.0, 
      2.0/301.0, 
      2.0/303.0, 
      2.0/305.0, 
      2.0/307.0, 
      2.0/309.0, 
      2.0/311.0, 
      2.0/313.0, 
      2.0/315.0, 
      2.0/317.0, 
      2.0/319.0, 
      2.0/321.0, 
      2.0/323.0, 
      2.0/325.0, 
      2.0/327.0, 
      2.0/329.0, 
      2.0/331.0, 
      2.0/333.0, 
      2.0/335.0, 
      2.0/337.0, 
      2.0/339.0, 
      2.0/341.0, 
      2.0/343.0, 
      2.0/345.0, 
      2.0/347.0, 
      2.0/349.0, 
      2.0/351.0, 
      2.0/353.0, 
      2.0/355.0, 
      2.0/357.0, 
      2.0/359.0, 
      2.0/361.0, 
      2.0/363.0, 
      2.0/365.0, 
      2.0/367.0, 
      2.0/369.0, 
      2.0/371.0, 
      2.0/373.0, 
      2.0/375.0, 
      2.0/377.0, 
      2.0/379.0, 
      2.0/381.0, 
      2.0/383.0, 
      2.0/385.0, 
      2.0/387.0, 
      2.0/389.0, 
      2.0/391.0, 
      2.0/393.0, 
      2.0/395.0, 
      2.0/397.0, 
      2.0/399.0, 
      2.0/401.0, 
      2.0/403.0, 
      2.0/405.0, 
      2.0/407.0, 
      2.0/409.0, 
      2.0/411.0, 
      2.0/413.0, 
      2.0/415.0, 
      2.0/417.0, 
      2.0/419.0, 
      2.0/421.0, 
      2.0/423.0, 
      2.0/425.0, 
      2.0/427.0, 
      2.0/429.0, 
      2.0/431.0, 
      2.0/433.0, 
      2.0/435.0, 
      2.0/437.0, 
      2.0/439.0, 
      2.0/441.0, 
      2.0/443.0, 
      2.0/445.0, 
      2.0/447.0, 
      2.0/449.0, 
      2.0/451.0, 
      2.0/453.0, 
      2.0/455.0, 
      2.0/457.0, 
      2.0/459.0, 
      2.0/461.0, 
      2.0/463.0, 
      2.0/465.0, 
      2.0/467.0, 
      2.0/469.0, 
      2.0/471.0, 
      2.0/473.0, 
      2.0/475.0, 
      2.0/477.0, 
      2.0/479.0, 
      2.0/481.0, 
      2.0/483.0, 
      2.0/485.0, 
      2.0/487.0, 
      2.0/489.0, 
      2.0/491.0, 
      2.0/493.0, 
      2.0/495.0, 
      2.0/497.0, 
      2.0/499.0, 
      2.0/501.0, 
      2.0/503.0, 
      2.0/505.0, 
      2.0/507.0, 
      2.0/509.0, 
      2.0/511.0, 
      2.0/513.0, 
      2.0/515.0, 
      2.0/517.0, 
      2.0/519.0, 
      2.0/521.0, 
      2.0/523.0, 
      2.0/525.0, 
      2.0/527.0, 
      2.0/529.0, 
      2.0/531.0, 
      2.0/533.0, 
      2.0/535.0, 
      2.0/537.0, 
      2.0/539.0, 
      2.0/541.0, 
      2.0/543.0, 
      2.0/545.0, 
      2.0/547.0, 
      2.0/549.0, 
      2.0/551.0, 
      2.0/553.0, 
      2.0/555.0, 
      2.0/557.0, 
      2.0/559.0, 
      2.0/561.0, 
      2.0/563.0, 
      2.0/565.0, 
      2.0/567.0, 
      2.0/569.0, 
      2.0/571.0, 
      2.0/573.0, 
      2.0/575.0, 
      2.0/577.0, 
      2.0/579.0, 
      2.0/581.0, 
      2.0/583.0, 
      2.0/585.0, 
      2.0/587.0, 
      2.0/589.0, 
      2.0/591.0, 
      2.0/593.0, 
      2.0/595.0, 
      2.0/597.0, 
      2.0/599.0, 
      2.0/601.0, 
      2.0/603.0, 
      2.0/605.0, 
      2.0/607.0, 
      2.0/609.0, 
      2.0/611.0, 
      2.0/613.0, 
      2.0/615.0, 
      2.0/617.0, 
      2.0/619.0, 
      2.0/621.0, 
      2.0/623.0, 
      2.0/625.0, 
      2.0/627.0, 
      2.0/629.0, 
      2.0/631.0, 
      2.0/633.0, 
      2.0/635.0, 
      2.0/637.0, 
      2.0/639.0, 
      2.0/641.0, 
      2.0/643.0, 
      2.0/645.0, 
      2.0/647.0, 
      2.0/649.0, 
      2.0/651.0, 
      2.0/653.0, 
      2.0/655.0, 
      2.0/657.0, 
      2.0/659.0, 
      2.0/661.0, 
      2.0/663.0, 
      2.0/665.0, 
      2.0/667.0, 
      2.0/669.0, 
      2.0/671.0, 
      2.0/673.0, 
      2.0/675.0, 
      2.0/677.0, 
      2.0/679.0, 
      2.0/681.0, 
      2.0/683.0, 
      2.0/685.0, 
      2.0/687.0, 
      2.0/689.0, 
      2.0/691.0, 
      2.0/693.0, 
      2.0/695.0, 
      2.0/697.0, 
      2.0/699.0, 
      2.0/701.0, 
      2.0/703.0, 
      2.0/705.0, 
      2.0/707.0, 
      2.0/709.0, 
      2.0/711.0, 
      2.0/713.0, 
      2.0/715.0, 
      2.0/717.0, 
      2.0/719.0, 
      2.0/721.0, 
      2.0/723.0, 
      2.0/725.0, 
      2.0/727.0, 
      2.0/729.0, 
      2.0/731.0, 
      2.0/733.0, 
      2.0/735.0, 
      2.0/737.0, 
      2.0/739.0, 
      2.0/741.0, 
      2.0/743.0, 
      2.0/745.0, 
      2.0/747.0, 
      2.0/749.0, 
      2.0/751.0, 
      2.0/753.0, 
      2.0/755.0, 
      2.0/757.0, 
      2.0/759.0, 
      2.0/761.0, 
      2.0/763.0, 
      2.0/765.0, 
      2.0/767.0, 
      2.0/769.0, 
      2.0/771.0, 
      2.0/773.0, 
      2.0/775.0, 
      2.0/777.0, 
      2.0/779.0, 
      2.0/781.0, 
      2.0/783.0, 
      2.0/785.0, 
      2.0/787.0, 
      2.0/789.0, 
      2.0/791.0, 
      2.0/793.0, 
      2.0/795.0, 
      2.0/797.0, 
      2.0/799.0, 
      2.0/801.0, 
      2.0/803.0, 
      2.0/805.0, 
      2.0/807.0, 
      2.0/809.0, 
      2.0/811.0, 
      2.0/813.0, 
      2.0/815.0, 
      2.0/817.0, 
      2.0/819.0, 
      2.0/821.0, 
      2.0/823.0, 
      2.0/825.0, 
      2.0/827.0, 
      2.0/829.0, 
      2.0/831.0, 
      2.0/833.0, 
      2.0/835.0, 
      2.0/837.0, 
      2.0/839.0, 
      2.0/841.0, 
      2.0/843.0, 
      2.0/845.0, 
      2.0/847.0, 
      2.0/849.0, 
      2.0/851.0, 
      2.0/853.0, 
      2.0/855.0, 
      2.0/857.0, 
      2.0/859.0, 
      2.0/861.0, 
      2.0/863.0, 
      2.0/865.0, 
      2.0/867.0, 
      2.0/869.0, 
      2.0/871.0, 
      2.0/873.0, 
      2.0/875.0, 
      2.0/877.0, 
      2.0/879.0, 
      2.0/881.0, 
      2.0/883.0, 
      2.0/885.0, 
      2.0/887.0, 
      2.0/889.0, 
      2.0/891.0, 
      2.0/893.0, 
      2.0/895.0, 
      2.0/897.0, 
      2.0/899.0, 
      2.0/901.0, 
      2.0/903.0, 
      2.0/905.0, 
      2.0/907.0, 
      2.0/909.0, 
      2.0/911.0, 
      2.0/913.0, 
      2.0/915.0, 
      2.0/917.0, 
      2.0/919.0, 
      2.0/921.0, 
      2.0/923.0, 
      2.0/925.0, 
      2.0/927.0, 
      2.0/929.0, 
      2.0/931.0, 
      2.0/933.0, 
      2.0/935.0, 
      2.0/937.0, 
      2.0/939.0, 
      2.0/941.0, 
      2.0/943.0, 
      2.0/945.0, 
      2.0/947.0, 
      2.0/949.0, 
      2.0/951.0, 
      2.0/953.0, 
      2.0/955.0, 
      2.0/957.0, 
      2.0/959.0, 
      2.0/961.0, 
      2.0/963.0, 
      2.0/965.0, 
      2.0/967.0, 
      2.0/969.0, 
      2.0/971.0, 
      2.0/973.0, 
      2.0/975.0, 
      2.0/977.0, 
      2.0/979.0, 
      2.0/981.0, 
      2.0/983.0, 
      2.0/985.0, 
      2.0/987.0, 
      2.0/989.0, 
      2.0/991.0, 
      2.0/993.0, 
      2.0/995.0, 
      2.0/997.0, 
      2.0/999.0,
      2.0/1001.0,
      2.0/1003.0
   };

static
double
quad_series_logd_evaluate(
   register double x,
   void const * restrict __accelerator
   )
{
   static double const
      __log10 = 2.302585092994045901093614,             // log(10.0)
      __log_base10_of_2 = 0.301029995663981198017467;   // log10(2.0)
   static double const
      __pow10s[] = {
         1./10000., 1./1000., 1./100., 1./10., 
         1., 10., 100., 1000., 10000., 100000.
         };
   register int const
      number_of_iterations = (int) (ptrdiff_t) __accelerator;
   register double const
      __pow2_exponent = log2(x);
   register int const
      __pow10_exponent = (int) floor(__pow2_exponent * __log_base10_of_2);
   register double
      convergent = x / __pow10s[__pow10_exponent + 4];
   register double
      yp = (convergent - 1.0) / (convergent + 1.0);
   register double const
      y2 = (yp * yp);
   convergent = yp * __constants[0u];
   for(int i = 1u; i< number_of_iterations; ++i)
   {
      convergent += (yp *= y2) * __constants[i];
   }
   return
      __pow10_exponent * __log10 + convergent;
}

VARRLogDEvaluator
quad_series_logd(size_t number_of_iterations)
{
   VARRLogDEvaluator
      result;
   result.accelerator =
      (void *) (ptrdiff_t) (int) (number_of_iterations + 1u);
   result.logd = quad_series_logd_evaluate;
   result.logd_array = NULL;
   result.disallocate = quad_series_logd_disallocate;
   return
      result;
}
