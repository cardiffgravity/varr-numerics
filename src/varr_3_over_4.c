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

#include "varr_3_over_4.h"
#include "varr.h"

#include <math.h>
#include <stdlib.h>
#include <inttypes.h>

typedef struct tagSampling3Over4DAccelerator
{
   double
      start_x,
      step_x,
      step_x_inverse;
   double
      step_x_inverse_powers[64];
   size_t
      samples;
   double const *
      values;
} Sampling3Over4DAccelerator;

static
int
__linear_sampling_normalizing_threequartersd_disallocate(
   Sampling3Over4DAccelerator * accelerator
   )
{
   if(accelerator == NULL)
   {
      return 1;
   }
   free((void *) accelerator->values);
   accelerator->start_x = 0x7F800001;
   accelerator->step_x = 0x7F800001;
   accelerator->step_x_inverse = 0x7F800001;
   for(unsigned i = 0u; i< 64u; ++i)
   {
      accelerator->step_x_inverse_powers[i] = 0x7F800001;
   }
   accelerator->samples = (size_t) 0;
   accelerator->values = NULL;
   free(accelerator);
   return 0;
}

static
int
linear_sampling_normalizing_threequartersd_disallocate(
   VARR3Over4DEvaluator * evaluator
   )
{
   if(evaluator == NULL)
   {
      return 1;
   }
   Sampling3Over4DAccelerator * const
      accelerator = (Sampling3Over4DAccelerator *)
         evaluator->accelerator;
   return
      __linear_sampling_normalizing_threequartersd_disallocate(accelerator);
}

static
Sampling3Over4DAccelerator *
allocate_linear_sampling_normalizing_threequartersd(
   size_t number_of_samples
   )
{
   double const
      step_size = 1.0 / (double) (number_of_samples - (size_t) 1u);
   double * const
      values = (double *) 
         malloc(sizeof(double) * (number_of_samples + 2u)
         );
   for(size_t i = (size_t) 0; i< (number_of_samples + 2u); ++i)
   {
      double const
         x = step_size * (double) i,
         value = pow(x, 3./4.0);
      values[i] = value;
      continue;
   }
   Sampling3Over4DAccelerator * const
      result = (Sampling3Over4DAccelerator *) malloc(
         sizeof(Sampling3Over4DAccelerator)
         );
   result->start_x = 0.0;
   result->step_x = step_size;
   result->step_x_inverse = 1.0 / step_size;
   for(uint64_t i = 0u; i< 64u; ++i)
   {
      result->step_x_inverse_powers[i] = 
         result->step_x_inverse / (double) (((uint64_t) 1) << i);
   }
   result->samples = number_of_samples,
   result->values = values;
   return
      result;
}

#ifdef __INTEL_COMPILER
static const int
   __log2_table[128] = {
   0, 1, 59, 0, 0, 2, 60, 0, 48, 0, 54, 0, 0, 3, 61, 0, 0, 40, 49, 0, 28, 0, 
   55, 0, 34, 0, 43, 0, 0, 4, 62, 0, 0, 52, 38, 0, 0, 41, 50, 0, 0, 19, 29, 
   0, 0, 21, 56, 0, 0, 31, 35, 0, 0, 12, 44, 0, 15, 0, 23, 0, 0, 5, 63, 0, 58, 
   0, 0, 47, 53, 0, 0, 39, 0, 27, 0, 33, 42, 0, 0, 51, 37, 0, 0, 18, 0, 20, 0, 
   30, 0, 11, 0, 14, 22, 0, 0, 57, 46, 0, 0, 26, 32, 0, 0, 36, 17, 0, 0, 10, 
   13, 0, 0, 45, 25, 0, 0, 16, 9, 0, 0, 24, 0, 8, 0, 7, 6, 0, 0
   };
#endif

static const double
   // i: 0-100: (1lu << i) ** 3./4.
   __3o4_roots[] = {
                            1.,
         1.6817928305074290041,
         2.8284271247461902909,
         4.7568284600108841076,
                            8.,
         13.454342644059432033,
         22.627416997969522328,
         38.054627680087072861,
                           64.,
         107.63474115247545626,
         181.01933598375617862,
         304.43702144069658289,
                          512.,
         861.07792921980365008,
          1448.154687870049429,
         2435.4961715255726631,
                         4096.,
         6888.6234337584292007,
         11585.237502960395432,
         19483.969372204581305,
                        32768.,
         55108.987470067433605,
         92681.900023683163454,
         155871.75497763665044,
                       262144.,
         440871.89976053946884,
         741455.20018946530763,
         1246974.0398210932035,
                      2097152.,
         3526975.1980843157507,
          5931641.601515722461,
          9975792.318568745628,
                     16777216.,
         28215801.584674526006,
         47453132.812125779688,
         79806338.548549965024,
                    134217728.,
         225726412.67739620805,
         379625062.49700623751,
         638450708.38839972019,
                   1073741824.,
         1805811301.4191696644,
         3037000499.9760499001,
         5107605667.1071977615,
                   8589934592.,
         14446490411.353357315,
           24296003999.8083992,
         40860845336.857582092,
                  68719476736.,
         115571923290.82685852,
          194368031998.4671936,
         326886762694.86065674,
                 549755813888.,
         924575386326.61486816,
         1554944255987.7375488,
         2615094101558.8852539,
                4398046511104.,
         7396603090612.9189453,
         12439554047901.900391,
         20920752812471.082031,
               35184372088832.,
         59172824724903.351562,
         99516432383215.203125,
         167366022499768.65625
   };

#include "varr_floor_log2.h"

static
double
linear_sampling_normalizing_threequartersd_evaluate(
   register double __x,
   register void const * restrict __accelerator
   )
{
   register char const
      do_invert = (__x < 1.0);
   register double const
      x = do_invert ? 1.0/__x : __x;
#ifdef __INTEL_COMPILER
   //
   // An exponential split pow(x, 3/4.) = pow(x / 2**k, 3/4.) * 2**(k*3./4.),
   // for the least integer k such that x < 2**k, combined with linear
   // interpolation in the range pow([0-1], 3/4.) and a modified De-Bruijn-
   // like method.
   //
   register uint64_t
      __pow2_gt = (uint64_t) x;
   __pow2_gt |= __pow2_gt >> 1;
   __pow2_gt |= __pow2_gt >> 2;
   __pow2_gt |= __pow2_gt >> 4;
   __pow2_gt |= __pow2_gt >> 8;
   __pow2_gt |= __pow2_gt >> 16;
   __pow2_gt |= __pow2_gt >> 32;
   ++__pow2_gt;
   register double const
      normalized_x = x / (double) __pow2_gt;
   register Sampling3Over4DAccelerator const * const
      accelerator = (Sampling3Over4DAccelerator const *) 
         __accelerator;
   register double
      step_frac = (normalized_x / accelerator->step_x);
   register uint64_t const
      index = (uint64_t) step_frac;
   step_frac -= index;
   step_frac =
      accelerator->values[index] * (1. - step_frac)
    + accelerator->values[index + (uint64_t) 1] * step_frac;
   register double
      result =
      step_frac * __3o4_roots[
         __log2_table[
            (--__pow2_gt * 285673954929480801ul) >> 57
            ]
         ];
#else // assume GCC:
   register unsigned const
      __pow2_exponent = 64u - __builtin_clzl((uint64_t) x);
   Sampling3Over4DAccelerator const * restrict const
      accelerator = (Sampling3Over4DAccelerator const *) 
         __accelerator;
   register double
      step_frac = (x * accelerator->step_x_inverse_powers[__pow2_exponent]);
   register uint64_t const
      index = (uint64_t) step_frac;
   step_frac -= index;
   step_frac =
      accelerator->values[index] * (1. - step_frac)
    + accelerator->values[index + (uint64_t) 1] * step_frac;
   register double const
      result = step_frac * __3o4_roots[__pow2_exponent];
#endif
   return
      (do_invert ? 1.0/result : result);
}

VARR3Over4DEvaluator
linear_sampling_normalizing_3over4d(
   size_t number_of_samples
   )
{
   VARR3Over4DEvaluator
      result;
   Sampling3Over4DAccelerator * const
      accelerator = allocate_linear_sampling_normalizing_threequartersd(
         number_of_samples
         );
   result.accelerator = (void *) accelerator;
   result.threequartersd = linear_sampling_normalizing_threequartersd_evaluate;
   result.disallocate = linear_sampling_normalizing_threequartersd_disallocate;
   return
      result;
}
