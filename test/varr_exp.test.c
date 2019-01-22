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
VARRExpDEvaluator
   expd_evaluator;

static
double
expd_evaluate(double x)
{
   return
      expd_evaluator.expd(x, expd_evaluator.accelerator);
}

static
void
expd_machine_batch_evaluate(
   double const * in,
   double * out,
   size_t length
   )
{
   for(size_t i = 0u; i< length; ++i)
      out[i] = exp(in[i]);
   return;
}

static
void
expd_batch_evaluate(
   double const * in,
   double * out,
   size_t length
   )
{
   expd_evaluator.expd_array(
      in,
      out,
      length,
      expd_evaluator.accelerator
      );
   return;
}

static
double
expd_test_no_remainder_loop(void)
{
   double
      numerical_error = -1.0;
   
   expd_evaluator = 
      shifting_linear_sampling_expd(3000000u);
   
   printf("Scalar evaluation:\n");
   numerical_error = 
      evaluate_performanced(
         -10.0,
         +10.0,
         50000000u,
         0,
         exp,
         expd_evaluate
         );
   
   // batch evaluation:
   printf("Vector evaluation:\n");
   numerical_error = fmax(
      numerical_error,
      evaluate_batch_performanced(
         -10.0,
         +10.0,
         50000000u,
         0,
         expd_machine_batch_evaluate,
         expd_batch_evaluate,
         0
         )
      );
   
   expd_evaluator.disallocate(&expd_evaluator);
   
   return
      numerical_error;
}

static
double
expd_test_with_remainder_loop(void)
{
   double
      numerical_error = -1.0;
   
   expd_evaluator = 
      shifting_linear_sampling_expd(3000000u);
   
   printf("Scalar evaluation:\n");
   numerical_error = 
      evaluate_performanced(
         -10.0,
         +10.0,
         50000003u,
         0,
         exp,
         expd_evaluate
         );
   
   // batch evaluation:
   printf("Vector evaluation:\n");
   numerical_error = fmax(
      numerical_error,
      evaluate_batch_performanced(
         -10.0,
         +10.0,
         50000003u,
         0,
         expd_machine_batch_evaluate,
         expd_batch_evaluate,
         0
         )
      );
   
   expd_evaluator.disallocate(&expd_evaluator);
   
   return
      numerical_error;
}

static
double
expd_test_with_remainder_loop_in_place_transform(void)
{
   double
      numerical_error = -1.0;
   
   expd_evaluator = 
      shifting_linear_sampling_expd(3000000u);
   
   // batch evaluation:
   printf("Vector evaluation:\n");
   numerical_error = fmax(
      numerical_error,
      evaluate_batch_performanced(
         -10.0,
         +10.0,
         50000003u,
         0,
         expd_machine_batch_evaluate,
         expd_batch_evaluate,
         1
         )
      );
   
   expd_evaluator.disallocate(&expd_evaluator);
   
   return
      numerical_error;
}

UnitTestResult
test_varr_exp(void)
{
   UnitTestResult
      result = create_test_results();
   
   printf("exp(x) numerical tests:\n");
   
   static double const
      worst_allowed_numerical_error = 1.47e-14;
   
   {
   UnitTestResult
      test_outcome = create_test_results();
   declare_start_of_unit_test();
   double const
      numerical_error = expd_test_no_remainder_loop();
   test_outcome.test_message = 
      create_message_specific_to_numerical_error_test_case(
         "varr-exp/exp(x)",
         "Sampling evaluation (no remainder loop)",
         numerical_error,
         worst_allowed_numerical_error
         );
   update_test_results_for_numerical_error_test_case(
      &test_outcome, numerical_error, worst_allowed_numerical_error
      );
   declare_end_of_unit_test();
   combine_test_results(test_outcome, &result);
   }
   
   {
   UnitTestResult
      test_outcome = create_test_results();
   declare_start_of_unit_test();
   double const
      numerical_error = expd_test_with_remainder_loop();
   test_outcome.test_message = 
      create_message_specific_to_numerical_error_test_case(
         "varr-exp/exp(x)",
         "Sampling evaluation (with remainder loop)",
         numerical_error,
         worst_allowed_numerical_error
         );
   update_test_results_for_numerical_error_test_case(
      &test_outcome, numerical_error, worst_allowed_numerical_error
      );
   declare_end_of_unit_test();
   combine_test_results(test_outcome, &result);
   }
   
   {
   UnitTestResult
      test_outcome = create_test_results();
   declare_start_of_unit_test();
   double const
      numerical_error = expd_test_with_remainder_loop_in_place_transform();
   test_outcome.test_message = 
      create_message_specific_to_numerical_error_test_case(
         "varr-exp/exp(x)",
         "Sampling evaluation (with remainder loop, in-place)",
         numerical_error,
         worst_allowed_numerical_error
         );
   update_test_results_for_numerical_error_test_case(
      &test_outcome, numerical_error, worst_allowed_numerical_error
      );
   declare_end_of_unit_test();
   combine_test_results(test_outcome, &result);
   }
   
   return
      result;
}
