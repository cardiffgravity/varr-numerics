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
#include "timings.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

static
VARRLogDEvaluator
   logd_evaluator;

static
double
logd_evaluate(double x)
{
   return
      logd_evaluator.logd(x, logd_evaluator.accelerator);
}

static
void
logd_machine_batch_evaluate(
   double const * restrict in,
   double * restrict out,
   size_t length
   )
{
   for(size_t i = 0u; i< length; ++i)
      out[i] = log(in[i]);
   return;
}

static
void
logd_batch_evaluate(
   double const * restrict in,
   double * restrict out,
   size_t length
   )
{
   logd_evaluator.logd_array(
      in,
      out,
      length,
      logd_evaluator.accelerator
      );
   return;
}

static
double
logd_test(void)
{
   double
      numerical_error = -1.;
   
   printf("Case 1: normalizing sampling log:\n");
   logd_evaluator = 
      normalizing_linear_sampling_logd(100000u);
   numerical_error = fmax(
      numerical_error,
      evaluate_performanced(
         1.e-10,
         1.e+10,
         10000000u,
         1,
         log,
         logd_evaluate
         )
      );
   logd_evaluator.disallocate(&logd_evaluator);
   
   printf("Case 2: quad series convergent:\n");
   logd_evaluator = 
      quad_series_logd(100u);
   numerical_error = fmax(
      numerical_error,
      evaluate_performanced(
         1.e-3,
         1.e+3,
         10000000u,
         1,
         log,
         logd_evaluate
         )
      );
   logd_evaluator.disallocate(&logd_evaluator);
   
   printf("Case 3: batch normalizing sampling log:\n");
   logd_evaluator = 
      normalizing_linear_sampling_logd(100000u);
   numerical_error = fmax(
      numerical_error,
      evaluate_batch_performanced(
         1.e-10,
         1.e+10,
         50000000u,
         1,
         logd_machine_batch_evaluate,
         logd_batch_evaluate,
         1
         )
      ); 
   logd_evaluator.disallocate(&logd_evaluator);
   
   return
      numerical_error;
}

static
double
sublinear_logd_test(void)
{
   double
      numerical_error = -1.;

   printf("Case 4: batch sublinear sampling log:\n");
   logd_evaluator = 
      normalizing_sublinear_sampling_logd(5000000u);
   numerical_error = fmax(
      numerical_error,
      evaluate_performanced(
         1.e-10,
         1.e+10,
         10000000u,
         1,
         log,
         logd_evaluate
         )
      );
   logd_evaluator.disallocate(&logd_evaluator);
   
   printf("Case 5: batch sublinear sampling log:\n");
   logd_evaluator = 
      normalizing_sublinear_sampling_logd(5000000u);
   numerical_error = fmax(
      numerical_error,
      evaluate_batch_performanced(
         1.e-10,
         1.e+10,
         50000000u,
         1,
         logd_machine_batch_evaluate,
         logd_batch_evaluate,
         1
         )
      ); 
   logd_evaluator.disallocate(&logd_evaluator);
   
   return
      numerical_error;
}

UnitTestResult
test_varr_log(void)
{
   UnitTestResult
      result = create_test_results();
   
   {
   declare_start_of_unit_test();
   printf("log(x) numerical tests:\n");
   static double const
      worst_allowed_numerical_error = 9.63e-10;
   double const
      numerical_error = logd_test();
   result.test_message = create_message_specific_to_numerical_error_test_case(
      "varr-log/natural-logarithm",
      "Sampling evaluation",
      numerical_error,
      worst_allowed_numerical_error
      );
   update_test_results_for_numerical_error_test_case(
      &result, numerical_error, worst_allowed_numerical_error
      );
   declare_end_of_unit_test();
   }
   
   {
   declare_start_of_unit_test();
   printf("log(x) (sublinear) numerical tests:\n");
   static double const
      worst_allowed_numerical_error = 1.56e-05;
   double const
      numerical_error = sublinear_logd_test();
   result.test_message = create_message_specific_to_numerical_error_test_case(
      "varr-log/natural-logarithm",
      "Sublinear sampling evaluation",
      numerical_error,
      worst_allowed_numerical_error
      );
   update_test_results_for_numerical_error_test_case(
      &result, numerical_error, worst_allowed_numerical_error
      );
   declare_end_of_unit_test();
   }
   
   return
      result;
}
