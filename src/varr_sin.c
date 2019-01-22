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

#include "varr_sin.h"
#include "varr_internal.h"

#include <math.h>
#include <stdlib.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

typedef struct tagSamplingSinDAccelerator
{
   double
      start_x,
      step_x;
   size_t
      samples;
   double const *
      values;
} SamplingSinDAccelerator;

static
int
__sampling_sind_disallocate(
   SamplingSinDAccelerator * accelerator
   )
{
   if(accelerator == NULL)
   {
      return 1;
   }
   free((void *) accelerator->values);
   accelerator->start_x = 0x7F800001;
   accelerator->step_x = 0x7F800001;
   accelerator->samples = (size_t) 0;
   accelerator->values = NULL;
   free(accelerator);
   return 0;
}

static
int
sampling_sind_disallocate(
   VARRSinDEvaluator * evaluator
   )
{
   if(evaluator == NULL)
   {
      return 1;
   }
   SamplingSinDAccelerator * const
      accelerator = (SamplingSinDAccelerator *)
         evaluator->accelerator;
   return
      __sampling_sind_disallocate(accelerator);
}

static double const
   __2pi = 2.0 * M_PI;

static
SamplingSinDAccelerator *
allocate_sampling_sind(
   size_t number_of_samples
   )
{
   double const
      step_size = __2pi / (double) (number_of_samples - (size_t) 1u);
   double * const
      values = (double *) malloc(sizeof(double) * number_of_samples);
   for(size_t i = (size_t) 0; i< number_of_samples; ++i)
   {
      double const
         x = step_size * (double) i,
         value = sin(x);
      values[i] = value;
      continue;
   }
   SamplingSinDAccelerator * const
      result = (SamplingSinDAccelerator *) malloc(
         sizeof(SamplingSinDAccelerator)
         );
   result->start_x = 0.0;
   result->step_x = step_size;
   result->samples = number_of_samples,
   result->values = values;
   return
      result;
}

static
double
sampling_sind_evaluate(
   register double x,
   register void const * __accelerator
   )
{
   if((x <= 0.0) || (x >= __2pi))
   {
      return 0.0;
   }
   register SamplingSinDAccelerator const * const
      accelerator = (SamplingSinDAccelerator const *) __accelerator;
   register size_t const
      index = (size_t) (x / accelerator->step_x);
   register double
      step_frac = x / accelerator->step_x;
   step_frac -= floor(step_frac);
   return
      (
      accelerator->values[index] * (1. - step_frac) + 
    + accelerator->values[index + 1] * step_frac
      );
}

VARRSinDEvaluator
sampling_sind(size_t number_of_samples)
{
   VARRSinDEvaluator
      result;
   SamplingSinDAccelerator * const
      accelerator = allocate_sampling_sind(
         number_of_samples
         );
   result.accelerator = (void *) accelerator;
   result.sind = sampling_sind_evaluate;
   result.disallocate = sampling_sind_disallocate;
   return
      result;
}

typedef struct tagCubicSplineSamplingSinDAccelerator
{
   SamplingSinDAccelerator const *
      base;
   gsl_interp_accel *
      gsl_accelerator;
   gsl_spline *
      gsl_spline;
   double const *
      x_nodes;
} CubicSplineSamplingSinDAccelerator;

static
int
cubic_spline_sampling_sind_disallocate(
   VARRSinDEvaluator * evaluator
   )
{
   if(evaluator == NULL)
   {
      return 1;
   }
   CubicSplineSamplingSinDAccelerator * const
      accelerator = (CubicSplineSamplingSinDAccelerator *)
         evaluator->accelerator;
   if(accelerator == NULL)
   {
      return 1;
   }
   int
      base_disallocator_result =
         __sampling_sind_disallocate(
            (SamplingSinDAccelerator *)
            accelerator->base
            );
   if(base_disallocator_result != 0)
   {
      return base_disallocator_result;
   }
   
   gsl_spline_free(accelerator->gsl_spline);
   accelerator->gsl_spline = NULL;
   
   gsl_interp_accel_free(accelerator->gsl_accelerator);
   accelerator->gsl_accelerator = NULL;
   
   free((double *) accelerator->x_nodes);
   accelerator->x_nodes = NULL;
   
   free(accelerator);
   evaluator->accelerator = NULL;
   
   return 0;
}

static
double
cubic_spline_sampling_sind_evaluate(
   double x,
   void const * __accelerator
   )
{
   if((x <= 0.0) || (x >= __2pi))
   {
      return 0.0;
   }
   CubicSplineSamplingSinDAccelerator const * const
      accelerator = (CubicSplineSamplingSinDAccelerator const *) 
         __accelerator;
   return
      gsl_spline_eval(
         accelerator->gsl_spline,
         x,
         accelerator->gsl_accelerator
         );
}

VARRSinDEvaluator
cubic_spline_sampling_sind(size_t number_of_samples)
{
   VARRSinDEvaluator
      result;
   
   CubicSplineSamplingSinDAccelerator * const
      cubic_accelerator = 
         (CubicSplineSamplingSinDAccelerator *)
            malloc(sizeof(CubicSplineSamplingSinDAccelerator));
   SamplingSinDAccelerator * const
      base_accelerator = allocate_sampling_sind(
         number_of_samples
         );
   
   cubic_accelerator->base = base_accelerator;
   
   cubic_accelerator->gsl_accelerator = gsl_interp_accel_alloc();
   cubic_accelerator->gsl_spline = 
      gsl_spline_alloc(gsl_interp_cspline, number_of_samples);
   
   double * const
      x_nodes = (double *) malloc(sizeof(double) * number_of_samples);
   
   for(size_t i = (size_t) 0u; i< number_of_samples; ++i)
   {
      x_nodes[i] = base_accelerator->step_x * i;
   }
   
   cubic_accelerator->x_nodes = x_nodes;
   
   gsl_spline_init(
      cubic_accelerator->gsl_spline,
      x_nodes,
      cubic_accelerator->base->values,
      number_of_samples
      );
   
   result.accelerator = (void *) cubic_accelerator;
   result.sind = cubic_spline_sampling_sind_evaluate;
   result.disallocate = cubic_spline_sampling_sind_disallocate;
   return
      result;
}

