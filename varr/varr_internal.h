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

#ifndef __VARR_INTERNAL_H__
#define __VARR_INTERNAL_H__

#include "varr.h"

#ifdef __VARR_USE_AVX512__
#ifndef __VARR_HAS_AVX__
#define __VARR_HAS_AVX__
#endif
#endif

#ifdef __VARR_HAS_AVX__

#include <x86intrin.h>

#ifndef __VARR_USE_AVX512__
typedef
   __m256d
   avxd_array_t;
#define __AVX_DOUBLE_STRIDE__ ((size_t) 4u)
#define _avxd_stride_max _mm256_max_pd
#define _avxd_stride_min _mm256_min_pd
#define _avxd_stride_floor _mm256_floor_pd
#define _avxd_stride_set _mm256_set_pd
#define _avxd_stride_set_duplicates _mm256_set1_pd
#else
typedef
   __m512d
   avxd_array_t;
#define __AVX_DOUBLE_STRIDE__ ((size_t) 8u)
#define _avxd_stride_max _mm512_max_pd
#define _avxd_stride_min _mm512_min_pd
#define _avxd_stride_floor _mm512_floor_pd
#define _avxd_stride_set _mm512_set_pd
#define _avxd_stride_set_duplicates _mm512_set1_pd
#endif
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#endif /* __VARR_INTERNAL_H__ */
