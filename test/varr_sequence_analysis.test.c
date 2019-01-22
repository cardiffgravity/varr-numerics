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
#include "varr_test.h"
#include "sequence_generation.h"
#include "sequence_comparison.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct tagNumericalResult {
   double
      observed_numerical_difference,
      maximum_acceptable_numerical_difference;
} NumericalResult;

static
NumericalResult
test_linear_sequence_generation(void)
{
   NumericalResult
      result;
   result.maximum_acceptable_numerical_difference = 1.e-14;
   result.observed_numerical_difference = -1.0;
   double const * const
      sequence = 
         generate_samples(
            -10.,
            +10.,
            (size_t) 21u,
            0
            );
   double const
      expected_sequence[21] = {
         -10., -9., -8., -7., -6., -5., -4., -3., -2., -1., 0.,
         1., 2., 3., 4., 5., 6., 7., 8., 9., 10.
         };
   for(size_t i = (size_t) 0u; i< (size_t) 21u; ++i)
   {
      double const
         observed = sequence[i],
         expected = expected_sequence[i],
         absolute_difference = fabs(observed - expected);
      if(absolute_difference > result.observed_numerical_difference)
         result.observed_numerical_difference = absolute_difference;
      continue;
   }
   free((void *) sequence);
   printf(
      "Numerical discrepancy: observed: %g, maximum allowed: %g\n",
      result.observed_numerical_difference,
      result.maximum_acceptable_numerical_difference
      );
   return
      result;
}

static
NumericalResult
test_log_sequence_generation(void)
{
   NumericalResult
      result;
   result.maximum_acceptable_numerical_difference = 1.e-10;
   result.observed_numerical_difference = -1.0;
   double const * const
      sequence = 
         generate_samples(
            1.e-2,
            1.e3,
            (size_t) 5u,
            1
            );
   double const
      expected_sequence[5] = {
         1.e-2,
         0.177827941,
         3.16227766,
         56.23413252,
         1000.
         };
   for(size_t i = (size_t) 0u; i< (size_t) 5u; ++i)
   {
      double const
         observed = sequence[i],
         expected = expected_sequence[i],
         relative_difference = 
            fabs(1. - fmin(observed, expected) / fmax(observed, expected));
      if(relative_difference > result.observed_numerical_difference)
         result.observed_numerical_difference = relative_difference;
      continue;
   }
   free((void *) sequence);
   printf(
      "Numerical discrepancy: observed: %g, maximum allowed: %g\n",
      result.observed_numerical_difference,
      result.maximum_acceptable_numerical_difference
      );
   return
      result;
}

static
NumericalResult
test_exponential_sequence_comparison(void)
{
   NumericalResult
      result;
   result.maximum_acceptable_numerical_difference = 1.e-10;
   result.observed_numerical_difference = -1.0;
   double const
      first_sequence[5] = {
         -1., 2., -3., 9., -5.e-5
         },
      second_sequence[5] = {
         6., -7., -8., 4., 10.e-5
         };
   double const
      observed = 
         evaluate_numerical_differerence_exponential_weights(
            first_sequence,
            second_sequence,
            (size_t) 5u,
            1.e3
            );
   /*
         Expected calculation:
         1.   |1.0 - (-1 / 6.)|         ~ 1.1666666..
         2.   |1.0 - (2. / -7.)|        ~ 1.285714286..
         3.   |1.0 - (-3. / -8.)|       ~ 0.625
         4.   |1.0 - (4. / 9.)|         ~ 0.5555555..
         5.   |1.0 - (-5e-5 / .10e-5)|  ~ 1.5
         absmedian(first sequence) = 2.
         Exponential weight (w): (alpha / 2.)**2 = 250000..
         Significances:
         1.   1.-exp(-w*36)      ~ 1.0
         2.   1.-exp(-w*53)      ~ 1.0
         3.   1.-exp(-w*73)      ~ 1.0
         4.   1.-exp(-w*97)      ~ 1.0
         5.   1.-exp(-w*1.25e-8) ~ 3.0120122268e-3
         Weighed differences:
         1.   1.1666666.               *
         2.   1.285714286.             *
         3.   1.375
         4.   0.5555555
         5.   1.5 * 3.0120122268e-3 = 4.6801..e-3
              -------------------
              1.285714286 = 9/7
    */
   static double const
      expected = 9./7.;
   result.observed_numerical_difference = fabs(observed - expected);
   result.maximum_acceptable_numerical_difference = 5.e-16;
   printf(
      "Numerical discrepancy: observed: %g, maximum allowed: %g\n",
      result.observed_numerical_difference,
      result.maximum_acceptable_numerical_difference
      );
   return
      result;
}

UnitTestResult
test_varr_sequence_analysis(void)
{
   UnitTestResult
      result = create_test_results();
   
   {
   UnitTestResult
      test_outcome = create_test_results();
   declare_start_of_unit_test();
   printf("sequence analysis numerical tests (linear generation):\n");
   NumericalResult const
      report = test_linear_sequence_generation();
   test_outcome.test_message =
      create_message_specific_to_numerical_error_test_case(
         "varr-sequence-analysis",
         "sequence generation (linear)",
         report.observed_numerical_difference,
         report.maximum_acceptable_numerical_difference
         );
   update_test_results_for_numerical_error_test_case(
      &test_outcome,
      report.observed_numerical_difference,
      report.maximum_acceptable_numerical_difference
      );
   declare_end_of_unit_test();
   combine_test_results(test_outcome, &result);
   }
   
   {
   UnitTestResult
      test_outcome = create_test_results();
   declare_start_of_unit_test();
   printf("Sequence analysis numerical tests (logarithmic generation):\n");
   NumericalResult const
      report = test_log_sequence_generation();
   test_outcome.test_message =
      create_message_specific_to_numerical_error_test_case(
         "varr-sequence-analysis",
         "sequence generation (logscale)",
         report.observed_numerical_difference,
         report.maximum_acceptable_numerical_difference
         );
   update_test_results_for_numerical_error_test_case(
      &test_outcome,
      report.observed_numerical_difference,
      report.maximum_acceptable_numerical_difference
      );
   declare_end_of_unit_test();
   combine_test_results(test_outcome, &result);
   }
   
   {
   UnitTestResult
      test_outcome = create_test_results();
   declare_start_of_unit_test();
   printf("Sequence comparison numerical tests (exponential weight):\n");
   NumericalResult const
      report = test_exponential_sequence_comparison();
   test_outcome.test_message =
      create_message_specific_to_numerical_error_test_case(
         "varr-sequence-analysis",
         "sequence comparison (exponential)",
         report.observed_numerical_difference,
         report.maximum_acceptable_numerical_difference
         );
   update_test_results_for_numerical_error_test_case(
      &test_outcome,
      report.observed_numerical_difference,
      report.maximum_acceptable_numerical_difference
      );
   declare_end_of_unit_test();
   combine_test_results(test_outcome, &result);
   }
   
   return
      result;
}
