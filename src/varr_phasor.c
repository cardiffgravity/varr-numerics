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

#include "varr_phasor.h"
#include "varr_internal.h"
#include "varr_sin.h"
#include "varr_cos.h"
#include "math.h"

#include <stdlib.h>

typedef struct tagCubicSplineInterpolatingPhasorDAccelerator
{
   VARRSinDEvaluator
      sind_evaluator;
   VARRCosDEvaluator
      cosd_evaluator;
} CubicSplineInterpolatingPhasorDAccelerator;

static
int
__cubic_spline_sampling_phasord_disallocate(
   CubicSplineInterpolatingPhasorDAccelerator * accelerator
   )
{
   if(accelerator == NULL)
   {
      return 1;
   }
   int
      result = accelerator->sind_evaluator.disallocate(
         &accelerator->sind_evaluator
         );
   result &= accelerator->cosd_evaluator.disallocate(
      &accelerator->cosd_evaluator
      );
   free(accelerator);
   return
      result;
}

static
int
cubic_spline_sampling_phasord_disallocate(
   VARRPhasorDEvaluator * evaluator
   )
{
   if(evaluator == NULL)
   {
      return 1;
   }
   CubicSplineInterpolatingPhasorDAccelerator * const
      accelerator = (CubicSplineInterpolatingPhasorDAccelerator *)
         evaluator->accelerator;
   return
      __cubic_spline_sampling_phasord_disallocate(accelerator);
}

static
CubicSplineInterpolatingPhasorDAccelerator *
allocate_cubic_spline_sampling_phasord(
   size_t number_of_samples
   )
{
   CubicSplineInterpolatingPhasorDAccelerator * const
      result = (CubicSplineInterpolatingPhasorDAccelerator *) malloc(
         sizeof(CubicSplineInterpolatingPhasorDAccelerator)
         );
   result->sind_evaluator = 
      cubic_spline_sampling_sind(number_of_samples);
   result->cosd_evaluator = 
      cubic_spline_sampling_cosd(number_of_samples);
   return
      result;
}

static
double complex
cubic_spline_sampling_phasord_evaluate(
   double phi,
   void const * __accelerator
   )
{
   CubicSplineInterpolatingPhasorDAccelerator const * const
      accelerator = (CubicSplineInterpolatingPhasorDAccelerator const *)
         __accelerator;
   double complex const
      result = 
         accelerator->cosd_evaluator.cosd(
            phi,
            accelerator->cosd_evaluator.accelerator
            )
       + accelerator->sind_evaluator.sind(
            phi,
            accelerator->sind_evaluator.accelerator
            ) * I;
   return
      result;
}

VARRPhasorDEvaluator
cubic_spline_sampling_phasord(size_t number_of_samples)
{
   VARRPhasorDEvaluator
      result;
   CubicSplineInterpolatingPhasorDAccelerator * const
      accelerator = allocate_cubic_spline_sampling_phasord(
         number_of_samples
         );
   result.accelerator = (void *) accelerator;
   result.phasord = cubic_spline_sampling_phasord_evaluate;
   result.phasord_array = NULL;
   result.disallocate = cubic_spline_sampling_phasord_disallocate;
   return
      result;
}

#include "varr_internal.h"

typedef struct 
   tagLinearInterpolatingPhasorDAccelerator
{
   double
      start_phi,
      step_phi,
      step_phi_inv;
   size_t
      samples;
   double complex const * __attribute__((aligned(128)))
      values;
} __attribute__((aligned(128))) LinearInterpolatingPhasorDAccelerator;

static
int
__linear_interpolating_phasord_disallocate(
   LinearInterpolatingPhasorDAccelerator * accelerator
   )
{
   if(accelerator == NULL)
   {
      return 1;
   }
   free((void *) (accelerator->values - 2u));
   accelerator->start_phi = 0x7F800001;
   accelerator->step_phi = 0x7F800001;
   accelerator->step_phi_inv = 0x7F800001;
   accelerator->samples = (size_t) 0;
   accelerator->values = NULL;
   free(accelerator);
   return 0;
}

static
int
linear_interpolating_phasord_disallocate(
   VARRPhasorDEvaluator * evaluator
   )
{
   if(evaluator == NULL)
   {
      return 1;
   }
   LinearInterpolatingPhasorDAccelerator * const
      accelerator = (LinearInterpolatingPhasorDAccelerator *)
         evaluator->accelerator;
   return
      __linear_interpolating_phasord_disallocate(accelerator);
}

static double const
   __2pi = 2.0 * M_PI;

static inline
double
fmod_phase_2pi(double phase)
{
   return 
      (phase >= 0.0) ? fmod(phase, __2pi) : (__2pi + fmod(phase, __2pi));
}

static inline
double complex
__linear_interpolating_phasord_evaluate(
   register double phi,
   register LinearInterpolatingPhasorDAccelerator const * accelerator
   )
{
   phi = fmod_phase_2pi(phi);
   register double
      step_frac = phi * accelerator->step_phi_inv;
   register size_t const
      index = (size_t) step_frac;
   step_frac -= floor(step_frac);
   return
      (
      accelerator->values[index] * (1. - step_frac) + 
    + accelerator->values[index + 1u] * step_frac
      );
}

static
double complex
linear_interpolating_phasord_evaluate(
   register double phi,
   register void const * __accelerator
   )
{
   return
      __linear_interpolating_phasord_evaluate(
         phi,
         (LinearInterpolatingPhasorDAccelerator const *) 
            __accelerator
         );
}

static
LinearInterpolatingPhasorDAccelerator *
allocate_linear_interpolating_PhasorD(
   size_t number_of_samples
   )
{
   double const
      step_size = __2pi / (double) number_of_samples,
      step_size_inv = 1.0 / step_size;
   double complex *
      values = (double complex *) malloc(
         sizeof(double complex) * (number_of_samples + 4u)
         );
   values[0] = cexp(I * (-2.*step_size));
   values[1] = cexp(I * (-step_size));
   values += 2u;
   for(size_t i = (size_t) 0; i<= (number_of_samples + 1u); ++i)
   {
      double const
         phi = step_size * (double) i;
      double complex const
         value = cexp(I * phi);
      values[i] = value;
      continue;
   }
   void *
      __result;
   int const
      alloc_result = 
         posix_memalign(
            &__result, 64,
            sizeof(LinearInterpolatingPhasorDAccelerator)
            );
   if(alloc_result)
      return NULL;
   LinearInterpolatingPhasorDAccelerator * const
      result = (LinearInterpolatingPhasorDAccelerator *) __result;
   result->start_phi = 0.0;
   result->step_phi = step_size;
   result->step_phi_inv = step_size_inv;
   result->samples = number_of_samples,
   result->values = values;
   return
      result;
}

#include <stdio.h>

static
void
linear_interpolating_phasord_batch_evaluate(
   register double const * restrict x,
   register double complex * restrict out,
   register size_t length,
   register void const * restrict __accelerator
   )
{
   register LinearInterpolatingPhasorDAccelerator const * const
      accelerator = 
          ((LinearInterpolatingPhasorDAccelerator const *) __accelerator);
#ifndef __VARR_HAS_AVX__
   for(register size_t i = (size_t) 0u; i< length; ++i, ++x)
   {
      out[i] = __linear_interpolating_phasord_evaluate(
         *x,
         accelerator
         );
   }
#else
   static double const
      __1_over_2pi = 1.0 / (2.0 * M_PI);
   register double complex const * const
      values = accelerator->values;
   register size_t const
      number_of_samples = accelerator->samples;
   register size_t const
      length_avx_stride = length / __AVX_DOUBLE_STRIDE__;
   register avxd_array_t
      alpha0,
      alpha1;
   register avxd_array_t const
      __one = _avxd_stride_set_duplicates(1.);
   register avxd_array_t const
      step_phi_inv = _avxd_stride_set_duplicates(accelerator->step_phi_inv);
   double __attribute__((aligned(128)))
      workspace[__AVX_DOUBLE_STRIDE__ * (size_t) 2u];
   double complex __attribute__((aligned(128)))
      alignment_emulator[__AVX_DOUBLE_STRIDE__];
   register avxd_array_t
      * const lower = (avxd_array_t *) workspace,
      * const upper = (lower + 1u);
   register avxd_array_t
      * __out = (avxd_array_t *) alignment_emulator;
   register double complex
      * export_target = out;
#ifndef __INTEL_COMPILER
   register
#endif
   size_t
      indexing[__AVX_DOUBLE_STRIDE__];
   register double
      x0, x1;
#ifdef __VARR_USE_AVX512__
   register double
      x2, x3;
#endif
   for(
      register size_t i = 0u;
      (i++) < length_avx_stride;
      x += __AVX_DOUBLE_STRIDE__
      )
   {
#ifdef __VARR_USE_AVX512__
      x3 = x[3] - (x[3] < 0.) * ((int64_t) (x[3] * __1_over_2pi) - 1) * __2pi;
      x2 = x[2] - (x[2] < 0.) * ((int64_t) (x[2] * __1_over_2pi) - 1) * __2pi;
#endif
      x1 = x[1] - (x[1] < 0.) * ((int64_t) (x[1] * __1_over_2pi) - 1) * __2pi;
      x0 = x[0] - (x[0] < 0.) * ((int64_t) (x[0] * __1_over_2pi) - 1) * __2pi;
      *lower = _avxd_stride_set(
#ifdef __VARR_USE_AVX512__
         x3, x3,
         x2, x2,
#endif
         x1, x1,
         x0, x0
         ); 
      *lower *= step_phi_inv;
#ifdef __VARR_USE_AVX512__
      x3 = x[7] - (x[7] < 0.) * ((int64_t) (x[7] * __1_over_2pi) - 1) * __2pi;
      x2 = x[6] - (x[6] < 0.) * ((int64_t) (x[6] * __1_over_2pi) - 1) * __2pi;
      x1 = x[5] - (x[5] < 0.) * ((int64_t) (x[5] * __1_over_2pi) - 1) * __2pi;
      x0 = x[4] - (x[4] < 0.) * ((int64_t) (x[4] * __1_over_2pi) - 1) * __2pi;
#else
      x1 = x[3] - (x[3] < 0.) * ((int64_t) (x[3] * __1_over_2pi) - 1) * __2pi;
      x0 = x[2] - (x[2] < 0.) * ((int64_t) (x[2] * __1_over_2pi) - 1) * __2pi;
#endif
      *upper = _avxd_stride_set(
#ifdef __VARR_USE_AVX512__
         x3, x3,
         x2, x2,
#endif
         x1, x1,
         x0, x0
         ); 
      *upper *= step_phi_inv;
      indexing[0u] = ((size_t) workspace[0u]) % number_of_samples;
      indexing[1u] = ((size_t) workspace[2u]) % number_of_samples;
      indexing[2u] = ((size_t) workspace[4u]) % number_of_samples;
      indexing[3u] = ((size_t) workspace[6u]) % number_of_samples;
#ifdef __VARR_USE_AVX512__
      indexing[4u] = ((size_t) workspace[8u]) % number_of_samples;
      indexing[5u] = ((size_t) workspace[10u]) % number_of_samples;
      indexing[6u] = ((size_t) workspace[12u]) % number_of_samples;
      indexing[7u] = ((size_t) workspace[14u]) % number_of_samples;
#endif
      alpha0 = *lower - _avxd_stride_floor(*lower);
      alpha1 = *upper - _avxd_stride_floor(*upper);
#ifdef __VARR_USE_AVX512__
      *(double complex *) workspace        = values[indexing[0u]++];
      *(double complex *) (workspace + 2u) = values[indexing[1u]++];
      *(double complex *) (workspace + 4u) = values[indexing[2u]++];
      *(double complex *) (workspace + 6u) = values[indexing[3u]++];
      *(double complex *) (workspace + 8u) = values[indexing[0u]];
      *(double complex *) (workspace + 10u) = values[indexing[1u]];
      *(double complex *) (workspace + 12u) = values[indexing[2u]];
      *(double complex *) (workspace + 14u) = values[indexing[3u]];
      *__out = *upper * alpha0;
      *__out = (*__out + (*lower * (__one - alpha0)));
      *(double complex *) workspace        = values[indexing[4u]++];
      *(double complex *) (workspace + 2u) = values[indexing[5u]++];
      *(double complex *) (workspace + 4u) = values[indexing[6u]++];
      *(double complex *) (workspace + 6u) = values[indexing[7u]++];
      *(double complex *) (workspace + 8u) = values[indexing[4u]];
      *(double complex *) (workspace + 10u) = values[indexing[5u]];
      *(double complex *) (workspace + 12u) = values[indexing[6u]];
      *(double complex *) (workspace + 14u) = values[indexing[7u]];
      __out[1] = *upper * alpha1;
      __out[1] = (__out[1] + (*lower * (__one - alpha1)));
#else
      *(double complex *) workspace        = values[indexing[0u]++];
      *(double complex *) (workspace + 2u) = values[indexing[1u]++];
      *(double complex *) (workspace + 4u) = values[indexing[0u]];
      *(double complex *) (workspace + 6u) = values[indexing[1u]];
      *__out = *upper * alpha0;
      *__out = *__out + (*lower * (__one - alpha0));
      *(double complex *) workspace        = values[indexing[2u]++];
      *(double complex *) (workspace + 2u) = values[indexing[3u]++];
      *(double complex *) (workspace + 4u) = values[indexing[2u]];
      *(double complex *) (workspace + 6u) = values[indexing[3u]];
      __out[1] = *upper * alpha1;
      __out[1] = __out[1] + (*lower * (__one - alpha1));
#endif
      export_target[0] = alignment_emulator[0];
      export_target[1] = alignment_emulator[1];
      export_target[2] = alignment_emulator[2];
      export_target[3] = alignment_emulator[3];
#ifdef __VARR_USE_AVX512__
      export_target[4] = alignment_emulator[4];
      export_target[5] = alignment_emulator[5];
      export_target[6] = alignment_emulator[6];
      export_target[7] = alignment_emulator[7];
#endif
      export_target += __AVX_DOUBLE_STRIDE__;
      continue;
   }
   
   // Remainder loop:
   
   for(
      register size_t i = length_avx_stride * __AVX_DOUBLE_STRIDE__;
      i< length;
      ++i, ++x
      )
   {
      out[i] = __linear_interpolating_phasord_evaluate(
         *x,
         accelerator
         );
   }
#endif
   return;
}

VARRPhasorDEvaluator
linear_interpolating_phasord(size_t number_of_samples)
{
   VARRPhasorDEvaluator
      result;
   LinearInterpolatingPhasorDAccelerator * const
      accelerator = allocate_linear_interpolating_PhasorD(
         number_of_samples
         );
   result.accelerator = (void *) accelerator;
   result.phasord = linear_interpolating_phasord_evaluate;
   result.phasord_array = linear_interpolating_phasord_batch_evaluate;
   result.disallocate = linear_interpolating_phasord_disallocate;
   return
      result;
}
