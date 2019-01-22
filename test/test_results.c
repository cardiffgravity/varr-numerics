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

#include <string.h>

UnitTestResult
create_test_results(void)
{
   UnitTestResult
      result;
   
   result.number_of_tests_evaluated = (size_t) 0u;
   result.number_of_tests_failed = (size_t) 0u;
   
   result.test_message = NULL;
   
   return
      result;
}

static size_t
   __test_index = (size_t) 1u;

#include <stdio.h>

void
declare_start_of_unit_test(void)
{
   fflush(stdout);
   printf("-------- Test %lu Begins --------\n", __test_index);
   fflush(stdout);
   return;
}

void
declare_end_of_unit_test(void)
{
   fflush(stdout);
   printf("-------- Test %lu Ends --------\n", __test_index++);
   fflush(stdout);
   return;
}

void
destroy_test_results(UnitTestResult * test_result)
{
   if(test_result->test_message)
   {
      free(test_result->test_message);
      test_result->test_message = NULL;
   }
   return;
}

void
combine_test_results(
   UnitTestResult new_results,
   UnitTestResult * restrict target
   )
{
   target->number_of_tests_evaluated += 
      new_results.number_of_tests_evaluated;
   target->number_of_tests_failed +=
      new_results.number_of_tests_failed;
   
   if(target->test_message == NULL)
   {
      target->test_message = new_results.test_message;
   }
   else if(new_results.test_message != NULL)
   {
      char const * const
         previous_message = target->test_message;
      char const * const
         new_message = new_results.test_message;
      size_t const
         previous_message_len = strlen(previous_message),
         new_message_len = strlen(new_message);
      target->test_message = 
         (char *) malloc(
            previous_message_len
          + new_message_len
          + 2u
            );
      strcpy(target->test_message, previous_message);
      target->test_message[previous_message_len] = '\n';
      char *
         cursor = (target->test_message + previous_message_len + 1u);
      strcpy(cursor, new_message);
      cursor += new_message_len;
      *cursor++ = '\0';
      free((void *) previous_message);
      free((void *) new_message);
      new_results.test_message = NULL;
   }
   
   return;
}

#include <stdio.h>

void
print_test_results(UnitTestResult const * result)
{
   printf(
      "-------- Result --------\n"
      "Number of tests evaluated: %lu\n"
      "Number of failed tests: %lu\n"
      "Notes and messages:\n%s\n"
      ,
      result->number_of_tests_evaluated,
      result->number_of_tests_failed,
      (result->test_message == NULL) ? "none" : result->test_message
      );
   return;
}

char *
create_unit_test_message(
   char const * test_group_name,
   char const * unit_test_name,
   char const * message
   )
{
   size_t const
      len_test_group_name = 
         strlen(test_group_name),
      len_unit_test_name = 
         strlen(unit_test_name),
      len_message =
         strlen(message);
   char *
      result = 
         (char *) malloc(
            len_test_group_name
          + len_unit_test_name
          + len_message
          + 5u
            );
   strcpy(result, test_group_name);
   char *
      cursor = result + len_test_group_name;
   *cursor++ = ':';
   *cursor++ = ' ';
   strcpy(cursor, unit_test_name);
   cursor += len_unit_test_name;
   *cursor++ = ':';
   *cursor++ = ' ';
   strcpy(cursor, message);
   return
      result;
}

char *
create_message_specific_to_numerical_error_test_case(
   char const * test_group_name,
   char const * unit_test_name,
   double observed_numerical_error,
   double worst_allowed_numerical_error
   )
{
   char
      message[1000u];
   if(observed_numerical_error >= worst_allowed_numerical_error)
   {
      sprintf(
         message,
         "*** Test failed: observed numerical error: %g; maximum allowed: %g. "
         "Threshold exceeded by %g%%.",
         observed_numerical_error,
         worst_allowed_numerical_error,
         100.0 * (observed_numerical_error / worst_allowed_numerical_error - 1.)
         );
   }
   else
   {
      sprintf(
         message,
         "Pass: observed numerical error: %g; maximum allowed: %g. "
         "Threshold missed by %g%%.",
         observed_numerical_error,
         worst_allowed_numerical_error,
         100.0 * (1. - observed_numerical_error / worst_allowed_numerical_error)
         );
   }
   return 
      create_unit_test_message(
         test_group_name,
         unit_test_name,
         message
         );
}

void
update_test_results_for_numerical_error_test_case(
   UnitTestResult * test_results,
   double observed_numerical_error,
   double worst_allowed_numerical_error
   )
{
   test_results->number_of_tests_evaluated++;
   if(observed_numerical_error >= worst_allowed_numerical_error)
      test_results->number_of_tests_failed++;
   return;
}
