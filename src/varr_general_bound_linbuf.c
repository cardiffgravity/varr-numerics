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

#include "varr_general_bound_linbuf.h"
#include "varr_internal.h"

#include <stdlib.h>
#include <math.h>
#include <inttypes.h>

typedef struct tagVARRBoundGLBEvaluator
{
   double
      start_x,
      end_x,
      step_x,
      step_x_inverse,
      interval_length,
      interval_length_inv,
      f_min_x,
      f_max_x,
      normalization;
   size_t
      samples;
   double const *
      linbuf_values;
   double (* delegate) (double, void *);
   void *
      delegate_argument;
} VARRBoundGLBEvaluator;

static
int
disallocate_evaluator(
   VARRBoundGLBEvaluator * accelerator
   )
{
   if(accelerator == NULL)
   {
      return 1;
   }
   
   static double const
      __none = 0x7F800001;
   
   if(accelerator->linbuf_values)
   {
      free((void *) accelerator->linbuf_values);
   }
   
   accelerator->start_x = __none;
   accelerator->end_x = __none;
   
   accelerator->step_x = __none;
   accelerator->step_x_inverse = __none;
   
   accelerator->interval_length = __none;
   accelerator->interval_length_inv = __none;
   
   accelerator->f_min_x = __none;
   accelerator->f_max_x = __none;
   
   accelerator->normalization = __none;
   
   accelerator->samples = (size_t) 0;
   accelerator->linbuf_values = NULL;
   accelerator->delegate = NULL;
   accelerator->delegate_argument = NULL;
   
   free(accelerator);
   
   return 0;
}

static
int
disallocate(
   VARRBoundGLBAccelerator * object
   )
{
   if(object == NULL)
   {
      return 1;
   }
   
   VARRBoundGLBEvaluator * const
      evaluator =
         (VARRBoundGLBEvaluator *) object->accelerator;
   
   return
      disallocate_evaluator(evaluator);
}

static
VARRBoundGLBEvaluator *
allocate(
   size_t number_of_samples,
   double min_x,
   double max_x,
   double (* delegate) (double x, void *),
   void * delegate_argument
   )
{
   double const
      interval_size = (max_x - min_x),
      step_size =
         interval_size / (double) (number_of_samples - (size_t) 1u);
   double * const
      linbuf_values =
         (double *)
         malloc(
            sizeof(double) * (number_of_samples + 2u)
            );
   for(
      size_t i = (size_t) 0;
      i< (number_of_samples + 2u);
      ++i
      )
   {
      double const
         x = min_x + step_size * (double) i,
         value = delegate(x, delegate_argument);
      linbuf_values[i] = value;
      
      continue;
   }
   
   VARRBoundGLBEvaluator * const
      result = (VARRBoundGLBEvaluator *) malloc(
         sizeof(VARRBoundGLBEvaluator)
         );
   
   result->start_x = min_x;
   result->end_x = max_x;
   
   result->step_x = step_size;
   result->step_x_inverse = 1.0 / step_size;
   
   result->interval_length = (max_x - min_x);
   result->interval_length_inv = 1./result->interval_length;
   
   result->f_min_x = delegate(min_x, delegate_argument);
   result->f_max_x = delegate(max_x, delegate_argument);
   
   result->normalization = result->step_x_inverse;
   
   result->samples = number_of_samples,
   result->linbuf_values = linbuf_values;
   result->delegate = delegate;
   result->delegate_argument = delegate_argument;
   
   return
      result;
}

static
double
evaluate_delegate(
   register double x,
   VARRBoundGLBEvaluator const * restrict accelerator
   )
{
   if(x <= accelerator->start_x)
   {
      return
         accelerator->f_min_x;
   }
   if(x >= accelerator->end_x)
   {
      return
         accelerator->f_max_x;
   }
   register double
      normalized_x =
         (x - accelerator->start_x) * accelerator->normalization;
   register int64_t
      index = (int64_t) floor(normalized_x);
   register double const * restrict
      value = (accelerator->linbuf_values + index);
   register double
      alpha = (normalized_x - index);
   return
      (*value * (1.0 - alpha) + value[1] * alpha);
}

static
double
evaluate_scalar(
   register double x,
   register void const * restrict accelerator
   )
{
   return
      evaluate_delegate(
         x,
         (VARRBoundGLBEvaluator const * const) accelerator
         );
}

static
void
batch_evaluate(
   register double const * __x,
   register double * out,
   register size_t length,
   register void const * restrict __accelerator
   )
{
   VARRBoundGLBEvaluator const * const
      accelerator =
         ((VARRBoundGLBEvaluator const *) __accelerator);
#ifndef __VARR_HAS_AVX__
   for(
      register size_t i = (size_t) 0u;
      i< length;
      ++i
      )
   {
      *out++ = evaluate_delegate(
         *__x++,
         accelerator
         );
      
      continue;
   }
   
   return;
#else
   register double const * const
      linbuf_values =
         ((VARRBoundGLBEvaluator const *) __accelerator)->linbuf_values;
   register size_t const
      length_axv_stride = length / __AVX_DOUBLE_STRIDE__;
   double __attribute__((aligned(128)))
      alignment_emulator[__AVX_DOUBLE_STRIDE__];
   avxd_array_t
      alpha,
      __linbuf_values1,
      __linbuf_values2,
    * const target = (avxd_array_t *) alignment_emulator;
   avxd_array_t
      __one = _avxd_stride_set_duplicates(1.0),
      avxd_lower_limit = _avxd_stride_set_duplicates(accelerator->start_x),
      avxd_upper_limit = _avxd_stride_set_duplicates(accelerator->end_x),
      avxd_normalization =
         _avxd_stride_set_duplicates(accelerator->normalization);
#ifndef __INTEL_COMPILER
   register
#endif
   int64_t
      index_avx[__AVX_DOUBLE_STRIDE__];
   for(
      register size_t i = (size_t) 0u;
      i< length_axv_stride;
      ++i
      )
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
      *target = _avxd_stride_max(*target, avxd_lower_limit);
      *target = _avxd_stride_min(*target, avxd_upper_limit);
      *target -= avxd_lower_limit;
      *target *= avxd_normalization;
      alpha = (*target - _avxd_stride_floor(*target));
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
      __linbuf_values1 = _avxd_stride_set(
#ifdef __VARR_USE_AVX512__
         linbuf_values[index_avx[7u]++],
         linbuf_values[index_avx[6u]++],
         linbuf_values[index_avx[5u]++],
         linbuf_values[index_avx[4u]++],
#endif
         linbuf_values[index_avx[3u]++],
         linbuf_values[index_avx[2u]++],
         linbuf_values[index_avx[1u]++],
         linbuf_values[index_avx[0u]++]
         );
      __linbuf_values2 = _avxd_stride_set(
#ifdef __VARR_USE_AVX512__
         linbuf_values[index_avx[7u]],
         linbuf_values[index_avx[6u]],
         linbuf_values[index_avx[5u]],
         linbuf_values[index_avx[4u]],
#endif
         linbuf_values[index_avx[3u]],
         linbuf_values[index_avx[2u]],
         linbuf_values[index_avx[1u]],
         linbuf_values[index_avx[0u]]
         );
      *target = __linbuf_values2 * alpha;
      alpha = (__one - alpha);
      *target = (*target + __linbuf_values1 * alpha);
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
   
   //
   // Remainder loop:
   //
   
   for(
      register size_t i =
         (length_axv_stride * __AVX_DOUBLE_STRIDE__);
      i< length;
      ++i
      )
   {
      *out++ = evaluate_delegate(
         *__x++,
         accelerator
         );
      
      continue;
   }
   
   return;
#endif
}

VARRBoundGLBAccelerator
bound_general_linbuf(
   size_t number_of_samples,
   double min_x,
   double max_x,
   double (* delegate) (double, void *),
   void * delegate_argument
   )
{
   VARRBoundGLBAccelerator
      result;
   VARRBoundGLBEvaluator * const
      evaluator =
         allocate(
            number_of_samples,
            min_x,
            max_x,
            delegate,
            delegate_argument
            );
   result.accelerator = (void *) evaluator;
   result.scalar = evaluate_scalar;
   result.batch = batch_evaluate;
   result.disallocate = disallocate;
   return
      result;
}
