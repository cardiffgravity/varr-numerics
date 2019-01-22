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

#define M_PI (3.14159265358979323846)

static
double
cos_ax_delegate(
   double x,
   void * p
   )
{
   return
      cos(x * *((double const *) p));
}

static
void *
cos_ax_delegate_fixed_p_value;

static
double
cos_ax_delegate_fixed_p(double x)
{
   return
      cos_ax_delegate(x, cos_ax_delegate_fixed_p_value);
}

static
void
cos_ax_delegate_fixed_p_machine_batch(
   double const * in,
   double * out,
   size_t length
   )
{
   for(size_t i = 0u; i< length; ++i)
   {
      out[i] = cos_ax_delegate_fixed_p(in[i]);
   }
   
   return;
}

static
VARRBoundGLBAccelerator
   delegate_evaluator;

static
double
delegate_evaluate(double x)
{
   return
      delegate_evaluator.scalar(
         x,
         delegate_evaluator.accelerator
         );
}

static
void
delegate_batch_evaluate(
   double const * in,
   double * out,
   size_t length
   )
{
   delegate_evaluator.batch(
      in,
      out,
      length,
      delegate_evaluator.accelerator
      );
   
   return;
}

static
double
delegate_test_no_remainder_loop(void)
{
   double
      numerical_error = -1.0;
   
   double const
      min_x = 0.,
      max_x = 2. * M_PI;
   
   double
      argument = 2.2;
   cos_ax_delegate_fixed_p_value = &argument;
   
   delegate_evaluator =
      bound_general_linbuf(
         3000000u,
         min_x,
         max_x,
         cos_ax_delegate,
         cos_ax_delegate_fixed_p_value
         );
   
   printf("Scalar evaluation:\n");
   
   numerical_error = fmax(
      numerical_error,
      evaluate_performanced(
         min_x,
         max_x,
         50000000u,
         0,
         cos_ax_delegate_fixed_p,
         delegate_evaluate
         )
      );
   
   //
   // batch evaluation:
   //
   
   printf("Vector evaluation:\n");
   
   numerical_error = fmax(
      numerical_error,
      evaluate_batch_performanced(
         min_x,
         max_x,
         50000000u,
         0,
         cos_ax_delegate_fixed_p_machine_batch,
         delegate_batch_evaluate,
         0
         )
      );
   
   delegate_evaluator.disallocate(&delegate_evaluator);
   
   return
      numerical_error;
}

static
double
delegate_test_with_remainder_loop(void)
{
   double
      numerical_error = -1.0;
   
   double const
      min_x = 0.,
      max_x = 2. * M_PI;
   
   double
      argument = 2.2;
   cos_ax_delegate_fixed_p_value = &argument;
   
   delegate_evaluator =
      bound_general_linbuf(
         3000000u,
         min_x,
         max_x,
         cos_ax_delegate,
         cos_ax_delegate_fixed_p_value
         );
   
   printf("Scalar evaluation:\n");
   
   numerical_error = fmax(
      numerical_error,
      evaluate_performanced(
         min_x,
         max_x,
         50000003u,
         0,
         cos_ax_delegate_fixed_p,
         delegate_evaluate
         )
      );
   
   //
   // batch evaluation:
   //
   
   printf("Vector evaluation:\n");
   
   numerical_error = fmax(
      numerical_error,
      evaluate_batch_performanced(
         min_x,
         max_x,
         50000003u,
         0,
         cos_ax_delegate_fixed_p_machine_batch,
         delegate_batch_evaluate,
         0
         )
      );
   
   delegate_evaluator.disallocate(&delegate_evaluator);
   
   return
      numerical_error;
}

static
double
delegate_test_with_remainder_loop_in_place_transform(void)
{
   double
      numerical_error = -1.0;
   
   double const
      min_x = 0.,
      max_x = 2. * M_PI;
   
   double
      argument = 2.2;
   cos_ax_delegate_fixed_p_value = &argument;
   
   delegate_evaluator =
      bound_general_linbuf(
         3000000u,
         min_x,
         max_x,
         cos_ax_delegate,
         cos_ax_delegate_fixed_p_value
         );
   
   printf("Scalar evaluation:\n");
   
   numerical_error = fmax(
      numerical_error,
      evaluate_performanced(
         min_x,
         max_x,
         50000003u,
         0,
         cos_ax_delegate_fixed_p,
         delegate_evaluate
         )
      );
   
   //
   // batch evaluation:
   //
   
   printf("Vector evaluation:\n");
   
   numerical_error = fmax(
      numerical_error,
      evaluate_batch_performanced(
         min_x,
         max_x,
         50000003u,
         0,
         cos_ax_delegate_fixed_p_machine_batch,
         delegate_batch_evaluate,
         1
         )
      );
   
   delegate_evaluator.disallocate(&delegate_evaluator);
   
   return
      numerical_error;
}

UnitTestResult
test_varr_general_bound_linbuf(void)
{
   UnitTestResult
      result = create_test_results();
   
   printf("general bound linbuf numerical tests:\n");
   
   {
   UnitTestResult
      test_outcome = create_test_results();
   declare_start_of_unit_test();
   static double const
      worst_allowed_numerical_error = 5.66e-12;
   double const
      numerical_error = delegate_test_no_remainder_loop();
   test_outcome.test_message = 
      create_message_specific_to_numerical_error_test_case(
         "varr general bound linbuf (cos specialization)",
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
   static double const
      worst_allowed_numerical_error = 5.66e-12;
   double const
      numerical_error = delegate_test_with_remainder_loop();
   test_outcome.test_message = 
      create_message_specific_to_numerical_error_test_case(
         "varr general bound linbuf (cos specialization)",
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
   static double const
      worst_allowed_numerical_error = 5.66e-12;
   double const
      numerical_error = delegate_test_with_remainder_loop_in_place_transform();
   test_outcome.test_message = 
      create_message_specific_to_numerical_error_test_case(
         "varr general bound linbuf (cos specialization)",
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
