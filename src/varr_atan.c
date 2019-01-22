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

#include "varr_atan.h"
#include "varr_internal.h"

#include <math.h>
#include <stdlib.h>

typedef struct tagSamplingAtanDAccelerator
{
   double
      start_x,
      step_x,
      step_x_inv;
   size_t
      samples;
   double const *
      values;
} SamplingAtanDAccelerator;

static
int
__sampling_atand_disallocate(
   SamplingAtanDAccelerator * accelerator
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
sampling_atand_disallocate(
   VARRAtanDEvaluator * evaluator
   )
{
   if(evaluator == NULL)
   {
      return 0;
   }
   SamplingAtanDAccelerator * const
      accelerator = (SamplingAtanDAccelerator *)
         evaluator->accelerator;
   return
      __sampling_atand_disallocate(accelerator);
}

static double const
   __sampling_atand_lower_limit = -50.,
   __sampling_atand_upper_limit = +50.;

static
SamplingAtanDAccelerator *
allocate_sampling_atand(
   size_t number_of_samples
   )
{
   double const
      range = (__sampling_atand_upper_limit - __sampling_atand_lower_limit),
      step_size = 
         range / (double) (number_of_samples - (size_t) 1u);
   ++number_of_samples;
   double * const
      values = (double *) malloc(sizeof(double) * number_of_samples);
   for(size_t i = (size_t) 0; i< number_of_samples; ++i)
   {
      double const
         x = __sampling_atand_lower_limit + step_size * (double) i,
         value = atan(x);
      values[i] = value;
      continue;
   }
   SamplingAtanDAccelerator * const
      result = (SamplingAtanDAccelerator *) malloc(
         sizeof(SamplingAtanDAccelerator)
         );
   result->start_x = __sampling_atand_lower_limit;
   result->step_x = step_size;
   result->step_x_inv = 1.0 / step_size;
   result->samples = number_of_samples,
   result->values = values;
   return
      result;
}

static inline
double
__sampling_atand_evaluate(
   register double x,
   SamplingAtanDAccelerator const * restrict accelerator
   )
{
   if(x < __sampling_atand_lower_limit)
   {
      x = __sampling_atand_lower_limit;
   }
   if(x > __sampling_atand_upper_limit)
   {
      x = __sampling_atand_upper_limit;
   }
   x -= __sampling_atand_lower_limit;
   register double
      step_frac = x * accelerator->step_x_inv;
   register uint64_t const
      index = (uint64_t) step_frac;
   step_frac -= floor(step_frac);
   double const * restrict
      __values = (accelerator->values + index);
   return
      *__values * (1. - step_frac) + __values[1] * step_frac;
}

static
double
sampling_atand_evaluate(
   register double x,
   register void const * restrict __accelerator
   )
{
   return
      __sampling_atand_evaluate(
         x,
         (SamplingAtanDAccelerator const *) __accelerator
         );
}

static
void
sampling_atand_batch_evaluate(
   register double const * __x,
   register double * out,
   register size_t length,
   register void const * restrict __accelerator
   )
{
#ifndef __VARR_HAS_AVX__
   for(register size_t i = (size_t) 0u; i< length; ++i)
   {
      *out++ = __sampling_atand_evaluate(
         *__x++,
         ((SamplingAtanDAccelerator const *) __accelerator)
         );
   }
#else
   register double const
      __step_x_inv =
         ((SamplingAtanDAccelerator const *) __accelerator)->step_x_inv;
   register double const * const
      values = ((SamplingAtanDAccelerator const *) __accelerator)->values;
   register size_t const
      length_axv_stride = length / __AVX_DOUBLE_STRIDE__;
#ifdef __INTEL_COMPILER
   for(size_t i = (size_t) 0u; i< length; ++i)
   {
      double
         x = __x[i];
      if(x < __sampling_atand_lower_limit)
         x = __sampling_atand_lower_limit;
      if(x > __sampling_atand_upper_limit)
         x = __sampling_atand_upper_limit;
      x -= __sampling_atand_lower_limit;
      register double
         step_frac = x * __step_x_inv;
      register unsigned const
         index = (unsigned) step_frac;
      step_frac -= floor(step_frac);
      double const * restrict
         __values = (values + index);
      out[i] =
         *__values * (1. - step_frac) + __values[1] * step_frac;
      /*
         note: icc -march=native -xHost:
         remark #15388: vectorization support: reference __x[i] has aligned 
                        access   [ ./src/varr_atan.c(148,14) ]
         remark #15388: vectorization support: reference out[i] has aligned 
                        access   [ ./src/varr_atan.c(161,7) ]
         remark #15411: vectorization support: conversion from float to int 
                        will be emulated   [ ./src/varr_atan.c(157,29) ]
         remark #15328: vectorization support: indirect load was emulated for 
                        the variable <*__values>, 64-bit indexed, part of 
                        address is read from memory   
                        [ ./src/varr_atan.c(162,11) ]
         remark #15328: vectorization support: indirect load was emulated for 
                        the variable <*(__values+8)>, 64-bit indexed, part of 
                        address is read from memory   
                        [ ./src/varr_atan.c(162,41) ]
         remark #15305: vectorization support: vector length 4
         remark #15309: vectorization support: normalized vectorization 
                        overhead 0.084
         remark #15300: LOOP WAS VECTORIZED
         remark #15448: unmasked aligned unit stride loads: 1 
         remark #15449: unmasked aligned unit stride stores: 1 
         remark #15462: unmasked indexed (or gather) loads: 2 
         remark #15475: --- begin vector cost summary ---
         remark #15476: scalar cost: 147 
         remark #15477: vector cost: 38.750 
         remark #15478: estimated potential speedup: 3.780 
         remark #15482: vectorized math library calls: 1 
         remark #15487: type converts: 1 
         remark #15488: --- end vector cost summary ---
      */
   }
#else       // assume GCC:
   double __attribute__((aligned(128)))
      alignment_emulator[__AVX_DOUBLE_STRIDE__];
   register avxd_array_t const
      __mavxd_sampling_atand_lower_limit = _avxd_stride_set(
#ifdef __VARR_USE_AVX512__
         __sampling_atand_lower_limit,
         __sampling_atand_lower_limit,
         __sampling_atand_lower_limit,
         __sampling_atand_lower_limit,
#endif
         __sampling_atand_lower_limit,
         __sampling_atand_lower_limit,
         __sampling_atand_lower_limit,
         __sampling_atand_lower_limit
         ),
      __mavxd_sampling_atand_upper_limit = _avxd_stride_set(
#ifdef __VARR_USE_AVX512__
         __sampling_atand_upper_limit,
         __sampling_atand_upper_limit,
         __sampling_atand_upper_limit,
         __sampling_atand_upper_limit,
#endif
         __sampling_atand_upper_limit,
         __sampling_atand_upper_limit,
         __sampling_atand_upper_limit,
         __sampling_atand_upper_limit
         ),
       __one = _avxd_stride_set(
#ifdef __VARR_USE_AVX512__
         1.0,
         1.0,
         1.0,
         1.0,
#endif
         1.0,
         1.0,
         1.0,
         1.0
         ),
       __step_x_inv_vector = _avxd_stride_set(
#ifdef __VARR_USE_AVX512__
         __step_x_inv,
         __step_x_inv,
         __step_x_inv,
         __step_x_inv,
#endif
         __step_x_inv,
         __step_x_inv,
         __step_x_inv,
         __step_x_inv
         );
   register avxd_array_t
      alpha,
      __values1,
      __values2,
    * const target = (avxd_array_t *) alignment_emulator;
   register uint64_t
      indexing[__AVX_DOUBLE_STRIDE__];
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
      *target = _avxd_stride_max(*target, __mavxd_sampling_atand_lower_limit);
      *target = _avxd_stride_min(*target, __mavxd_sampling_atand_upper_limit);
      *target = *target - __mavxd_sampling_atand_lower_limit;
      *target = *target * __step_x_inv_vector;
      indexing[0u] = (uint64_t) alignment_emulator[0u];
      indexing[1u] = (uint64_t) alignment_emulator[1u];
      indexing[2u] = (uint64_t) alignment_emulator[2u];
      indexing[3u] = (uint64_t) alignment_emulator[3u];
#ifdef __VARR_USE_AVX512__
      indexing[4u] = (uint64_t) alignment_emulator[4u];
      indexing[5u] = (uint64_t) alignment_emulator[5u];
      indexing[6u] = (uint64_t) alignment_emulator[6u];
      indexing[7u] = (uint64_t) alignment_emulator[7u];
#endif
      alpha = *target - _avxd_stride_floor(*target);
      __values1 = _avxd_stride_set(
#ifdef __VARR_USE_AVX512__
         values[indexing[7u]++],
         values[indexing[6u]++],
         values[indexing[5u]++],
         values[indexing[4u]++],
#endif
         values[indexing[3u]++],
         values[indexing[2u]++],
         values[indexing[1u]++],
         values[indexing[0u]++]
         );
      __values2 = _avxd_stride_set(
#ifdef __VARR_USE_AVX512__
         values[indexing[7u]],
         values[indexing[6u]],
         values[indexing[5u]],
         values[indexing[4u]],
#endif
         values[indexing[3u]],
         values[indexing[2u]],
         values[indexing[1u]],
         values[indexing[0u]]
         );
      *target = __values2 * alpha;
      alpha = __one - alpha;
      *target = *target + (__values1 * alpha);
      out[0u] = alignment_emulator[0u];
      out[1u] = alignment_emulator[1u];
      out[2u] = alignment_emulator[2u];
      out[3u] = alignment_emulator[3u];
#ifdef __VARR_USE_AVX512__
      out[4u] = alignment_emulator[4u];
      out[5u] = alignment_emulator[5u];
      out[6u] = alignment_emulator[6u];
      out[7u] = alignment_emulator[7u];
#endif
      __x += __AVX_DOUBLE_STRIDE__;
      out += __AVX_DOUBLE_STRIDE__;
      continue;
   }
   
   // Remainder loop:

   for(
      register size_t i = length_axv_stride * __AVX_DOUBLE_STRIDE__; 
      i< length; 
      ++i
      )
   {
      *out++ = __sampling_atand_evaluate(
         *__x++,
         ((SamplingAtanDAccelerator const *) __accelerator)
         );
   }
#endif
#endif
   return;
}

VARRAtanDEvaluator
clamping_linear_interpolating_atand(size_t number_of_samples)
{
   VARRAtanDEvaluator
      result;
   SamplingAtanDAccelerator * const
      accelerator = allocate_sampling_atand(
         number_of_samples
         );
   result.accelerator = (void *) accelerator;
   result.atan = sampling_atand_evaluate;
   result.atan_array = sampling_atand_batch_evaluate;
   result.disallocate = sampling_atand_disallocate;
   return
      result;
}
