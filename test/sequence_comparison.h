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

#ifndef __VARR_SEQUENCE_COMPARISON_H__
#define __VARR_SEQUENCE_COMPARISON_H__

#include <stddef.h>

double
evaluate_numerical_differerence_exponential_weights(
   double const * reference_sequence,
   double const * approximate_sequence,
   size_t sequence_length,
   double alpha_scale
   );

double
evaluate_numerical_differerence_exponential_weights_complex(
   double complex const * reference_sequence,
   double complex const * approximate_sequence,
   size_t sequence_length,
   double alpha_scale
   );

#endif /* __VARR_SEQUENCE_COMPARISON_H__ */
