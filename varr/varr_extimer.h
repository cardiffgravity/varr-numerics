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

#ifndef __VARR_EXTIMER_H__
#define __VARR_EXTIMER_H__

#include <time.h>

typedef struct tagVARRExecutionTimer {
   struct timespec
      begin_,
      end_;
} VARRExecutionTimer;

/*
 * Returns an object that includes timer data for a VARR timer operation.
 * When the returned object is provided as an argument to varr_end_timer,
 * that function returns the time in seconds elapsed since the creation of
 * its argument.
 */
VARRExecutionTimer
varr_start_timer(void);

double
varr_end_timer(VARRExecutionTimer);

#endif /* __VARR_EXTIMER_H__ */
