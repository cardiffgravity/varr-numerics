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
VARR3Over4DEvaluator
   threequartersd_evaluator;

static
double
threequartersd_evaluate(double x)
{
   return
      threequartersd_evaluator.threequartersd(
         x, threequartersd_evaluator.accelerator
         );
}

static
double
threequartersd_machine_evaluate(double x)
{
   return
      pow(x, 3./4.);
}

static
double
threequartersd_test(void)
{
   threequartersd_evaluator =
      linear_sampling_normalizing_3over4d(3000000u);
   double const
      worst_numerical_difference = 
         evaluate_performanced(
            1.e-18,
            1.e18,         // note 2**64 ~ 1.8e19
            100000000u,
            1,
            threequartersd_machine_evaluate,
            threequartersd_evaluate
            );
   threequartersd_evaluator.disallocate(&threequartersd_evaluator);
   return
      worst_numerical_difference;
}

UnitTestResult
test_varr_3_over_4(void)
{
   UnitTestResult
      result = create_test_results();
   
   {
   declare_start_of_unit_test();
   printf("pow(x, 3/4) numerical tests:\n");
   static double const
      worst_allowed_numerical_error = 1.12e-14;
   double const
      numerical_error = threequartersd_test();
   result.test_message = create_message_specific_to_numerical_error_test_case(
      "varr-3_over_4/pow(x, 3/4)",
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
