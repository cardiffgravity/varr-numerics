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
#include "sequence_generation.h"
#include "sequence_comparison.h"

#include <time.h>
#include <stdlib.h>
#include <stdio.h>

double
evaluate_performancedc(
   double range_x_lower,
   double range_x_upper,
   size_t number_of_samples,
   char use_logarithmic_sampling,
   double complex (* machine_stock_implementation) (double),
   double complex (* test_implementation) (double)
   )
{
   double complex
      * const standard_evaluations = 
         (double complex *) malloc(
            sizeof(double complex) * number_of_samples
            ),
      * const approximations =
         (double complex *) malloc(
            sizeof(double complex) * number_of_samples
            );
   double const * const
      x_samples = 
         generate_samples(
            range_x_lower, range_x_upper, 
            number_of_samples,
            use_logarithmic_sampling
            );
   {
   clock_t const
      begin = clock();
   for(size_t i = (size_t) 0u; i< number_of_samples; ++i)
   {
      double complex const
         v = machine_stock_implementation(x_samples[i]);
      standard_evaluations[i] = v;
      continue;
   }
   clock_t const
      end = clock();
   printf(
      "Timing: machine stock evaluation: %g\n",
      (end - begin) * (1000. / CLOCKS_PER_SEC)
      );
   }
   {
   clock_t const
      begin = clock();
   for(size_t i = (size_t) 0u; i< number_of_samples; ++i)
   {
      double complex const
         v = test_implementation(x_samples[i]);
      approximations[i] = v;
      continue;
   }
   clock_t const
      end = clock();
   printf(
      "Timing: custom evaluation: %g\n",
      (end - begin) * (1000. / CLOCKS_PER_SEC)
      );
   }
   double const
      worst_numerical_difference = 
         evaluate_numerical_differerence_exponential_weights_complex(
            standard_evaluations,
            approximations,
            number_of_samples,
            1.e3
            );
   printf(
      "Worst (relative) numerical difference: %g (%lu tests in the range "
      "[%g, %g], %s sampling)\n",
      worst_numerical_difference,
      number_of_samples,
      range_x_lower,
      range_x_upper,
      (use_logarithmic_sampling ? "logarithmic" : "linear")
      );
   free(standard_evaluations);
   free(approximations);
   free((void *) x_samples);
   return
      worst_numerical_difference;
}

double
evaluate_batch_performancedc(
   double range_x_lower,
   double range_x_upper,
   size_t number_of_samples,
   char use_logarithmic_sampling,
   void (* machine_stock_implementation) (
      double const *, double complex *, size_t
      ),
   void (* test_implementation) (
      double const *, double complex *, size_t
      )
   )
{
   void
      * __standard_evaluations,
      * __approximations;
   // Create poor alignments:
   posix_memalign(
      &__standard_evaluations, 64,
      sizeof(double complex) * (number_of_samples + 1)
      );
   posix_memalign(
      &__approximations, 64, sizeof(double complex) * 
      (number_of_samples + 2)
      );
   double complex
      * const standard_evaluations = 
         (double complex * const) (__standard_evaluations + 1),
      * const approximations = 
         (double complex * const) (__approximations + 2);
   double const * const
      x_samples = 
         generate_samples(
            range_x_lower, range_x_upper,
            number_of_samples,
            use_logarithmic_sampling
            );
   {
   clock_t const
      begin = clock();
   machine_stock_implementation(
      x_samples, standard_evaluations, number_of_samples
      );
   clock_t const
      end = clock();
   printf(
      "Timing: machine stock evaluation: %g\n",
      (end - begin) * (1000. / CLOCKS_PER_SEC)
      );
   }
   {
   clock_t const
      begin = clock();
   test_implementation(
      x_samples, approximations, number_of_samples
      );
   clock_t const
      end = clock();
   printf(
      "Timing: custom evaluation: %g\n",
      (end - begin) * (1000. / CLOCKS_PER_SEC)
      );
   }
   double const
      worst_numerical_difference = 
         evaluate_numerical_differerence_exponential_weights_complex(
            standard_evaluations,
            approximations,
            number_of_samples,
            1.e3
            );
   printf(
      "Worst (relative) numerical difference: %g (%lu tests in the range "
      "[%g, %g], %s sampling)\n",
      worst_numerical_difference,
      number_of_samples,
      range_x_lower,
      range_x_upper,
      (use_logarithmic_sampling ? "logarithmic" : "linear")
      );
   free(__standard_evaluations);
   free(__approximations);
   free((void *) x_samples);
   return
      worst_numerical_difference;
}
