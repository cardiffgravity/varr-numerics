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

#include "sequence_generation.h"

#include <stdlib.h>
#include <math.h>

double *
generate_samples(
   double range_x_lower,
   double range_x_upper,
   size_t number_of_samples,
   char use_logarithmic_sampling
   )
{
   void *
      __x_samples;
   posix_memalign(
      &__x_samples, 64, sizeof(double) * number_of_samples
      );
   double * const
      x_samples = (double *) __x_samples;
   if(use_logarithmic_sampling == 0)
   {
      double const
         x_step = 
            (range_x_upper - range_x_lower) / 
            (double) (number_of_samples - (size_t) 1);
      for(size_t i = (size_t) 0u; i< number_of_samples; ++i)
      {
         double const
            x = range_x_lower + x_step * i;
         x_samples[i] = x;
         continue;
      }
   }
   else
   {
      double const
         log_x_lower = log(range_x_lower),
         log_x_upper = log(range_x_upper);
      double const
         x_step = 
            (log_x_upper - log_x_lower) / 
            (double) (number_of_samples - (size_t) 1);
      for(size_t i = (size_t) 0u; i< number_of_samples; ++i)
      {
         double const
            x = log_x_lower + x_step * i;
         x_samples[i] = exp(x);
         continue;
      }
   }
   return
      x_samples;
}
