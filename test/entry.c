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

#include <stdlib.h>

int
main(int argc, char ** argv) {
   UnitTestResult
      result = create_test_results();
   
   combine_test_results(test_varr_general_bound_linbuf(), &result);
   
   combine_test_results(test_varr_3_over_4(), &result);
   
   combine_test_results(test_varr_sixthroot(), &result);
   
   combine_test_results(test_varr_sequence_analysis(), &result);
   
   combine_test_results(test_varr_log(), &result);
   
   combine_test_results(test_varr_atan(), &result);
   
   combine_test_results(test_varr_phasor(), &result);
   
   combine_test_results(test_varr_sin(), &result);
   
   combine_test_results(test_varr_exp(), &result);
   
   print_test_results(&result);
   
   destroy_test_results(&result);
   
   return
      EXIT_SUCCESS;
}
