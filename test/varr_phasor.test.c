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

#include "varr_test.h"
#include "varr_internal.h"
#include "timings.h"

#include "varr_phasor.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

static
VARRPhasorDEvaluator
   phasord_evaluator;

static
complex double
phasord_machine_evaluate(double x)
{
   return cexp(+ (1.0*I) * x);
}

static
complex double
phasord_evaluate(double x)
{
   return
      phasord_evaluator.phasord(x, phasord_evaluator.accelerator);
}

static
void
phasord_machine_batch_evaluate(
   double const * in,
   double complex * out,
   size_t length
   )
{
   for(size_t i = 0u; i< length; ++i)
      out[i] = cexp(+ (1.0*I) * in[i]);
   return;
}

static
void
phasord_batch_evaluate(
   double const * in,
   double complex * out,
   size_t length
   )
{
   phasord_evaluator.phasord_array(
      in,
      out,
      length,
      phasord_evaluator.accelerator
      );
   return;
}

static
double
phasord_test_no_remainder_loop(void)
{
   phasord_evaluator =
      linear_interpolating_phasord(3000000u);
      
   printf("Scalar evaluation:\n");
   double
      numerical_error = 
         evaluate_performancedc(
            -6.0 * M_PI,
            6.0 * M_PI,
            10000000u,
            0,
            phasord_machine_evaluate,
            phasord_evaluate
            );
   
   // batch evaluation:
   printf("Vector evaluation:\n");
   numerical_error = fmax(
      numerical_error,
         evaluate_batch_performancedc(
            -6.0 * M_PI,
            6.0 * M_PI,
            10000000u,
            0,
            phasord_machine_batch_evaluate,
            phasord_batch_evaluate
            )
         );
   
   phasord_evaluator.disallocate(&phasord_evaluator);
   
   return
      numerical_error;
}

static
double
phasord_test_with_remainder_loop(void)
{
   phasord_evaluator =
      linear_interpolating_phasord(3000000u);
      
   printf("Scalar evaluation:\n");
   double
      numerical_error = 
         evaluate_performancedc(
            -6.0 * M_PI,
            6.0 * M_PI,
            10000003u,
            0,
            phasord_machine_evaluate,
            phasord_evaluate
            );
   
   // batch evaluation:
   printf("Vector evaluation:\n");
   numerical_error = fmax(
      numerical_error,
         evaluate_batch_performancedc(
            -6.0 * M_PI,
            6.0 * M_PI,
            10000003u,
            0,
            phasord_machine_batch_evaluate,
            phasord_batch_evaluate
            )
         );
   
   phasord_evaluator.disallocate(&phasord_evaluator);
   
   return
      numerical_error;
}

UnitTestResult
test_varr_phasor(void)
{
   UnitTestResult
      result = create_test_results();
   
   printf("phasor(phi) = exp(i*phi) numerical tests:\n");
   
   {
   UnitTestResult
      test_outcome = create_test_results();
   declare_start_of_unit_test();
   static double const
      worst_allowed_numerical_error = 5.5e-13;
   double const
      numerical_error = phasord_test_no_remainder_loop();
   test_outcome.test_message = 
      create_message_specific_to_numerical_error_test_case(
         "varr-phasor/phasor(phi) = exp(i*phi)",
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
      worst_allowed_numerical_error = 5.5e-13;
   double const
      numerical_error = phasord_test_with_remainder_loop();
   test_outcome.test_message = 
      create_message_specific_to_numerical_error_test_case(
         "varr-phasor/phasor(phi) = exp(i*phi)",
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
   
   return
      result;
}
