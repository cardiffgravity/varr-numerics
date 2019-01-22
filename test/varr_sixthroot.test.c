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

static
VARRSixthRootDEvaluator
   sixthrootd_evaluator;

static
double
sixthrootd_evaluate(double x)
{
   return
      sixthrootd_evaluator.sixthrootd(x, sixthrootd_evaluator.accelerator);
}

static
double
sixthrootd_machine_evaluate(double x)
{
   return
      pow(x, 1.0/6.0);
}

static
void
sixthrootd_machine_batch_evaluate(
   double const * restrict in,
   double * restrict out,
   size_t length
   )
{
   for(size_t i = 0u; i< length; ++i)
      out[i] = pow(in[i], 1./6.);
   return;
}

static
void
sixthrootd_batch_evaluate(
   double const * restrict in,
   double * restrict out,
   size_t length
   )
{
   sixthrootd_evaluator.sixthrootd_array(
      in,
      out,
      length,
      sixthrootd_evaluator.accelerator
      );
   return;
}

static
double
sixthrootd_test(void)
{
   sixthrootd_evaluator =
      linear_sampling_normalizing_sixth_rootd(3000000u);
   double const
      worst_numerical_difference = 
         evaluate_performanced(
            1.e-18,
            1.e18,         // note 2**64 ~ 1.8e19
            50000000u,
            1,
            sixthrootd_machine_evaluate,
            sixthrootd_evaluate
            );
   sixthrootd_evaluator.disallocate(&sixthrootd_evaluator);
   return
      worst_numerical_difference;
}

static
double
sixthrootd_test_batch_evaluation(void)
{
   double
      numerical_error = -1.0;
   
   sixthrootd_evaluator =
      linear_sampling_normalizing_sixth_rootd(3000000u);
   
   // batch evaluation:
   printf("Vector evaluation:\n");
   numerical_error = fmax(
      numerical_error,
      evaluate_batch_performanced(
         1.e-18,
         1.e18,         // note 2**64 ~ 1.8e19
         50000000u,
         1,
         sixthrootd_machine_batch_evaluate,
         sixthrootd_batch_evaluate,
         0
         )
      );
   
   sixthrootd_evaluator.disallocate(&sixthrootd_evaluator);
   
   return
      numerical_error;
}

static
double
sixthrootd_test_batch_evaluation_with_remainder_loop(void)
{
   double
      numerical_error = -1.0;
   
   sixthrootd_evaluator =
      linear_sampling_normalizing_sixth_rootd(3000000u);
   
   // batch evaluation:
   printf("Vector evaluation:\n");
   numerical_error = fmax(
      numerical_error,
      evaluate_batch_performanced(
         1.e-18,
         1.e18,         // note 2**64 ~ 1.8e19
         50000003u,
         1,
         sixthrootd_machine_batch_evaluate,
         sixthrootd_batch_evaluate,
         0
         )
      );
   
   sixthrootd_evaluator.disallocate(&sixthrootd_evaluator);
   
   return
      numerical_error;
}

static
double
sixthrootd_test_batch_evaluation_with_remainder_loop_in_place(void)
{
   double
      numerical_error = -1.0;
   
   sixthrootd_evaluator =
      linear_sampling_normalizing_sixth_rootd(3000000u);
   
   // batch evaluation:
   printf("Vector evaluation:\n");
   numerical_error = fmax(
      numerical_error,
      evaluate_batch_performanced(
         1.e-18,
         1.e18,         // note 2**64 ~ 1.8e19
         50000003u,
         1,
         sixthrootd_machine_batch_evaluate,
         sixthrootd_batch_evaluate,
         1
         )
      );
   
   sixthrootd_evaluator.disallocate(&sixthrootd_evaluator);
   
   return
      numerical_error;
}

static
double
sixthrootd_sublinear_test(void)
{
   sixthrootd_evaluator =
      sublinear_sampling_normalizing_sixth_rootd(3000000u);
   double const
      worst_numerical_difference = 
         evaluate_performanced(
            1.e-18,
            1.e18,         // note 2**64 ~ 1.8e19
            50000000u,
            1,
            sixthrootd_machine_evaluate,
            sixthrootd_evaluate
            );
   sixthrootd_evaluator.disallocate(&sixthrootd_evaluator);
   return
      worst_numerical_difference;
}

static
double
sixthrootd_test_sublinear_batch_evaluation(void)
{
   double
      numerical_error = -1.0;
   
   sixthrootd_evaluator =
      sublinear_sampling_normalizing_sixth_rootd(3000000u);
   
   // batch evaluation:
   printf("Vector evaluation:\n");
   numerical_error = fmax(
      numerical_error,
      evaluate_batch_performanced(
         1.e-18,
         1.e18,         // note 2**64 ~ 1.8e19
         50000003u,
         1,
         sixthrootd_machine_batch_evaluate,
         sixthrootd_batch_evaluate,
         1
         )
      );
   
   sixthrootd_evaluator.disallocate(&sixthrootd_evaluator);
   
   return
      numerical_error;
}

UnitTestResult
test_varr_sixthroot(void)
{
   UnitTestResult
      result = create_test_results();
   
   static double const
      worst_allowed_numerical_error = 8.66e-15;
   
   printf("pow(x, 1/6) numerical tests:\n");
   
   {
   UnitTestResult
      test_outcome = create_test_results();
   declare_start_of_unit_test();
   double const
      numerical_error = sixthrootd_test();
   test_outcome.test_message = 
      create_message_specific_to_numerical_error_test_case(
         "varr-sixthroot/pow(x, 1/6)",
         "Sampling evaluation",
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
      numerical_error = sixthrootd_test_batch_evaluation();
   test_outcome.test_message = 
      create_message_specific_to_numerical_error_test_case(
         "varr-sixthroot/pow(x, 1/6) (batch evaluation)",
         "Sampling evaluation",
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
      numerical_error = sixthrootd_test_batch_evaluation_with_remainder_loop();
   test_outcome.test_message = 
      create_message_specific_to_numerical_error_test_case(
         "varr-sixthroot/pow(x, 1/6) (batch evaluation, with remainder loop)",
         "Sampling evaluation",
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
      numerical_error = 
         sixthrootd_test_batch_evaluation_with_remainder_loop_in_place();
   test_outcome.test_message = 
      create_message_specific_to_numerical_error_test_case(
         "varr-sixthroot/pow(x, 1/6) "
         "(batch evaluation, with remainder loop, in-place)",
         "Sampling evaluation",
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
   static double const
      worst_allowed_numerical_error = 1.12e-07;
   UnitTestResult
      test_outcome = create_test_results();
   declare_start_of_unit_test();
   double const
      numerical_error = 
         sixthrootd_sublinear_test();
   test_outcome.test_message = 
      create_message_specific_to_numerical_error_test_case(
         "varr-sixthroot/pow(x, 1/6) "
         "(scalar evaluation, with remainder loop, in-place, sublinear)",
         "Sampling evaluation",
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
   static double const
      worst_allowed_numerical_error = 1.12e-07;
   UnitTestResult
      test_outcome = create_test_results();
   declare_start_of_unit_test();
   double const
      numerical_error = 
         sixthrootd_test_sublinear_batch_evaluation();
   test_outcome.test_message = 
      create_message_specific_to_numerical_error_test_case(
         "varr-sixthroot/pow(x, 1/6) "
         "(batch evaluation, with remainder loop, in-place, sublinear)",
         "Sampling evaluation",
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
