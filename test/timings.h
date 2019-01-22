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
#ifndef __VARR_TIMINGS_H__
#define __VARR_TIMINGS_H__

#include <complex.h>

double
evaluate_performanced(
   double range_x_lower,
   double range_x_upper,
   size_t number_of_samples,
   char use_logarithmic_sampling,
   double (* machine_stock_implementation) (double),
   double (* test_implementation) (double)
   );

double
evaluate_batch_performanced(
   double range_x_lower,
   double range_x_upper,
   size_t number_of_samples,
   char use_logarithmic_sampling,
   void (* machine_stock_implementation) (
      double const *, double *, size_t
      ),
   void (* test_implementation) (
      double const *, double *, size_t
      ),
   unsigned char compute_in_place
   );

double
evaluate_performancedc(
   double range_x_lower,
   double range_x_upper,
   size_t number_of_samples,
   char use_logarithmic_sampling,
   complex double (* machine_stock_implementation) (double),
   complex double (* test_implementation) (double)
   );

double
evaluate_batch_performancedc(
   double range_x_lower,
   double range_x_upper,
   size_t number_of_samples,
   char use_logarithmic_sampling,
   void (* machine_stock_implementation) (
      double const *, double complex *, size_t
      ),
   void (* test_implementation) (
      double const *, double complex *, size_t
      )
   );

#endif /* __VARR_TIMINGS_H__ */
