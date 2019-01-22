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

#include "varr_sin.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

static
VARRSinDEvaluator
   sind_evaluator;

static
double
sind_evaluate(double x)
{
   return
      sind_evaluator.sind(x, sind_evaluator.accelerator);
}

static
double
sind_test(void)
{
   sind_evaluator = sampling_sind(500000u);
   double const
      numerical_error =
         evaluate_performanced(
            0.0,
            2.0 * M_PI,
            10000000u,
            0,
            sin,
            sind_evaluate
            );
   sind_evaluator.disallocate(&sind_evaluator);
   return
      numerical_error;
}

UnitTestResult
test_varr_sin(void)
{
   UnitTestResult
      result = create_test_results();
   
   {
   declare_start_of_unit_test();
   printf("sin(x) numerical tests:\n");
   static double const
      worst_allowed_numerical_error = 2.1e-11;
   double const
      numerical_error = sind_test();
   result.test_message = create_message_specific_to_numerical_error_test_case(
      "varr-sin/sin(x)",
      "Sampling evaluation",
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
