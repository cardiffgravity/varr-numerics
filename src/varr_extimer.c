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

#include "varr_extimer.h"

VARRExecutionTimer
varr_start_timer(void)
{
   VARRExecutionTimer
      result;
   clock_gettime(CLOCK_REALTIME, &result.begin_);
   return
      result;
}

double
varr_end_timer(VARRExecutionTimer from)
{
   clock_gettime(CLOCK_REALTIME, &from.end_);
   double const
      time_taken =
            ( from.end_.tv_sec - from.begin_.tv_sec )
          + ( from.end_.tv_nsec - from.begin_.tv_nsec ) * 1.e-9;
   return
      time_taken;
}
