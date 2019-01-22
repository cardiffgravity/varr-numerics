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

#ifndef __VARR_TEST_H__
#define __VARR_TEST_H__

#include <stdlib.h>

typedef struct tagUnitTestResult {
   size_t
      number_of_tests_evaluated,
      number_of_tests_failed;
   char *
      test_message;
} UnitTestResult;

UnitTestResult
create_test_results(void);

void
declare_start_of_unit_test(void);

void
declare_end_of_unit_test(void);

void
destroy_test_results(UnitTestResult *);

/*
 * Note that the first argument will be released, and assumed invalid,
 * following a call to this function:
 */
void
combine_test_results(UnitTestResult new_results, UnitTestResult * restrict);

void
print_test_results(UnitTestResult const *);

char *
create_unit_test_message(
   char const * test_group_name,
   char const * unit_test_name,
   char const * message
   );

char *
create_message_specific_to_numerical_error_test_case(
   char const * test_group_name,
   char const * unit_test_name,
   double observed_numerical_error,
   double worst_allowed_numerical_error
   );

void
update_test_results_for_numerical_error_test_case(
   UnitTestResult *,
   double observed_numerical_error,
   double worst_allowed_numerical_error
   );

UnitTestResult
test_varr_sixthroot(void);

UnitTestResult
test_varr_3_over_4(void);

UnitTestResult
test_varr_phasor(void);

UnitTestResult
test_varr_atan(void);

UnitTestResult
test_varr_log(void);

UnitTestResult
test_varr_sin(void);

UnitTestResult
test_varr_exp(void);

UnitTestResult
test_varr_general_bound_linbuf(void);

UnitTestResult
test_varr_sequence_analysis(void);

#endif /* __VARR_TEST_H__ */
