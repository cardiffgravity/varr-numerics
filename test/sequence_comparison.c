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

#include "varr.h"
#include "varr_utils.h"

#include "timings.h"

#include <time.h>
#include <float.h>
#include <stdlib.h>
#include <math.h>

// #define VERBOSE

#include <stdio.h>

static
int
__sort(void const * first, void const * second)
{
   double
      __first = *(double const *) first,
      __second = *(double const *) second;
   if(__first > __second)
   {
      return 1;
   }
   else if(__first < __second)
   {
      return -1;
   }
   else
   {
      return 0;
   }
}

double
evaluate_numerical_differerence_exponential_weights(
   double const * reference_sequence,
   double const * approximate_sequence,
   size_t sequence_length,
   double alpha_scale
   )
{
   double
      * const fabs_reference_sequence = 
         (double *) malloc(sizeof(double) * sequence_length),
      * const numerical_differences =
         (double *) malloc(sizeof(double) * sequence_length);
   for(size_t i = (size_t) 0u; i< sequence_length; ++i)
   {
      fabs_reference_sequence[i] = fabs(reference_sequence[i]);
   }
   qsort(fabs_reference_sequence, sequence_length, sizeof(double), __sort);
   double const
      median_abs = 
         0.5 * fabs_reference_sequence[(unsigned) floor(sequence_length/2)]
       + 0.5 * fabs_reference_sequence[(unsigned) ceil(sequence_length/2)];
   // Relative numerical differences:
   double const
      alpha = alpha_scale / median_abs,
      alpha_sq = alpha * alpha;
#ifdef VERBOSE
   printf(
      "Note: median scale (machine evaluations): %g, normalization: %g\n",
      median_abs,
      alpha
      );
#endif
   double 
      worst_numerical_difference = 0.;
#ifdef VERBOSE
   if(sequence_length > 0)
   {
      printf(
         "                    Ref. Value"
         "                  Approx. Value"
         "                     Num. Diff."
         "                       (Weight)"
         "                          Worst\n"
         );
   }
#endif
   for(size_t i = (size_t) 0u; i< sequence_length; ++i)
   {
      double
         relative_difference = 
            (  
               (reference_sequence[i] == 0.0) && 
               (approximate_sequence[i] == 0.0)
            ) ? 
            0.0 :
            (
            fabs(reference_sequence[i]) > fabs(approximate_sequence[i])) ?
               fabs(1.0 - approximate_sequence[i] / reference_sequence[i]) :
               fabs(1.0 - reference_sequence[i] / approximate_sequence[i]
            );
      double const
         weight = (
            1.0 - exp(
               -alpha_sq * reference_sequence[i] * reference_sequence[i]
               -alpha_sq * approximate_sequence[i] * approximate_sequence[i]
               )
            );
      relative_difference *= weight;
      numerical_differences[i] = relative_difference;
      if(worst_numerical_difference < relative_difference)
      {
         worst_numerical_difference = relative_difference;
      }
#ifdef VERBOSE
      printf(
         "%30.20g %30.20g %30.20g %30.20g %30.20g %s\n",
         reference_sequence[i],
         approximate_sequence[i],
         relative_difference,
         weight,
         worst_numerical_difference,
         ((worst_numerical_difference == relative_difference) ? 
            "*" : " ")
         );
#endif
      continue;
   }
   free(fabs_reference_sequence);
   free(numerical_differences);
   return
      worst_numerical_difference;
}

static
double
pow2(double x)
{
   return (x * x);
}

double
evaluate_numerical_differerence_exponential_weights_complex(
   double complex const * reference_sequence,
   double complex const * approximate_sequence,
   size_t sequence_length,
   double alpha_scale
   )
{
   double
      * const abs_reference_sequence = 
         (double *) malloc(sizeof(double) * sequence_length),
      * const numerical_differences =
         (double *) malloc(sizeof(double) * sequence_length);
   for(size_t i = (size_t) 0u; i< sequence_length; ++i)
   {
      abs_reference_sequence[i] = cabs(reference_sequence[i]);
   }
   qsort(abs_reference_sequence, sequence_length, sizeof(double), __sort);
   double const
      median_abs = 
         0.5 * abs_reference_sequence[(unsigned) floor(sequence_length/2)]
       + 0.5 * abs_reference_sequence[(unsigned) ceil(sequence_length/2)];
   // Relative numerical differences:
   double const
      alpha = alpha_scale / median_abs,
      alpha_sq = alpha * alpha;
#ifdef VERBOSE
   printf(
      "Note: median scale (machine evaluations): %g, normalization: %g\n",
      median_abs,
      alpha
      );
#endif
   double 
      worst_numerical_difference = 0.;
#ifdef VERBOSE
   if(sequence_length > 0)
      printf(
         "             Ref. Value (Re.)"
         "             Ref. Value (Im.)"
         "           Approx. Value (Re.)"
         "          Approx. Value (Im.)"
         "                   Num. Diff."
         "                     (Weight)"
         "                        Worst\n"
         );
#endif
   for(size_t i = (size_t) 0u; i< sequence_length; ++i)
   {
      double
         relative_difference = 
            (  
               (reference_sequence[i] == 0.0) && 
               (approximate_sequence[i] == 0.0)
            ) ? 
            0.0 :
            (
            cabs(reference_sequence[i]) > cabs(approximate_sequence[i])) ?
               cabs(1.0 - approximate_sequence[i] / reference_sequence[i]) :
               cabs(1.0 - reference_sequence[i] / approximate_sequence[i]
            );
      double const
         weight = (
            1.0 - exp(
               -alpha_sq * pow2(cabs(reference_sequence[i]))
               -alpha_sq * pow2(cabs(approximate_sequence[i]))
               )
            );
      relative_difference *= weight;
      numerical_differences[i] = relative_difference;
      if(worst_numerical_difference < relative_difference)
      {
         worst_numerical_difference = relative_difference;
      }
#ifdef VERBOSE
      printf(
         "%28.20g+%28.20gi %28.20g+%28.20gi %28.20g %28.20g %28.20g %s\n",
         creal(reference_sequence[i]),
         cimag(reference_sequence[i]),
         creal(approximate_sequence[i]),
         cimag(approximate_sequence[i]),
         relative_difference,
         weight,
         worst_numerical_difference,
         ((worst_numerical_difference == relative_difference) ? 
            "*" : " ")
         );
#endif
      continue;
   }
   free(abs_reference_sequence);
   free(numerical_differences);
   return
      worst_numerical_difference;
}
