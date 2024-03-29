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

#ifndef __VARR_PHASOR_H__
#define __VARR_PHASOR_H__

#include <stddef.h>
#include <complex.h>

/*
 * An object that provides a (double precision) VARR implementation of the 
 * complex phasor function.
 * 
 * This implementation provides scalar and vector complex phasor functions as 
 * follows:
 *    i.    'phasord' - a function accepting one real number, x, and returning 
 *          a VARR approximation of cexp(i*x);
 *    ii.   'phasord_array' - a vectorized version of (i.).
 */
typedef struct tagVARRPhasorDEvaluator {
   void const * accelerator;
   
   /*
    * A VARR complex phasor function.  The first argument to this method is the
    * number, x, for which cexp(i*x) is to be evaluated and returned.
    *
    * This function does not return meaningful values if the input, x, is not 
    * a finite real number.
    */
   double complex (* phasord) (double phi, void const * accelerator);
   
   /*
    * A vectorized VARR complex phasor function.  The meaning of arguments to 
    * this method are as follows:
    * 
    *    i..   'in' - a const array of doubles, being the values of x to 
    *          process;
    *    ii.   'out' - an array of double complex of (at least) the same length 
    *          as (i.), being the values of cexp(i*x) to assign;
    *    iii.  'length' - the length of the arrays (i.) and (ii.) above.
    *,   iv.   'accelerator' - the value of the above enclosed (void *) 
    *          accelerator.
    *
    * This function does not return meaningful values if the inputs (i.) are
    * not finite real numbers.  It is the responsibility of the caller to ensure
    * that the arrays (i.) and (ii.) both have size (iii.).  Note for 
    * allocation purposes that sizeof(double) is not equal to sizeof(double 
    * complex) in general.
    *
    * It is not necessary for the arrays (i.) or (ii.) to have any specific
    * alignments.  If this library is compiled without AVX extensions enabled,
    * this method will delegate to the above scalar VARR implementation of 
    * complex phasor.
    */
   void (* phasord_array) (
      double const * in, double complex * out,
      size_t length, 
      void const * accelerator
      );
   
   int (* disallocate) (struct tagVARRPhasorDEvaluator *);
   
} VARRPhasorDEvaluator;

/*
 * Returns an object that provides a VARR implementation of the complex
 * phasor function.   See documentation for the type VARRPhasorDEvaluator
 * for further information.
 * 
 * This function delegates to a pair of real-valued cubic spline evaluators, one
 * for the sine function and another for the cosine.  The value of the complex
 * phasor is assembled of the results of these delegates.  These delegates both 
 * allocate sampling grids of real (double) values.  The sizes
 * of these sampling grids are as indicated by the argument 'number_of_samples',
 * which must be nonzero.  The amount of memory allocated by this function is
 * approximately proportional to the number of sampling points requested.  The
 * numerical accuracy of the complex phasor function that is generated by this
 * method generally increases with the number of sampling points requested.
 *
 * See documentation of the VARR cubic spline sine evaluator, and of the VARR
 * cubic spline cosine evaluator, for further information.
 *
 * The evaluator returned by this method does not provide a vector
 * phasord_array function (VARRPhasorDEvaluator::phasord_array), and as such no 
 * calls to phasord_array should be made against the evaluator that is returned
 * by this method.
 */
VARRPhasorDEvaluator
cubic_spline_sampling_phasord(size_t number_of_samples);

/*
 * Returns an object that provides a VARR implementation of the complex
 * phasor function.   See documentation for the type VARRPhasorDEvaluator
 * for further information.
 * 
 * This function allocates a sampling grid of complex phasor values.  The size 
 * of this sampling grid is indicated by the argument 'number_of_samples', which
 * must be nonzero.  The amount of memory allocated by this function is
 * approximately proportional to the number of sampling points requested.  The
 * numerical accuracy of the complex phasor function that is generated by this
 * method generally increases with the number of sampling points requested.
 */
VARRPhasorDEvaluator
linear_interpolating_phasord(size_t number_of_samples);

#endif /* __VARR_PHASOR_H__ */
