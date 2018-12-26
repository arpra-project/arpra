/*
 * arpra-impl.h -- Private header file, used in compilation.
 *
 * Copyright 2016-2018 James Paul Turner.
 *
 * This file is part of the Arpra library.
 *
 * The Arpra library is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * The Arpra library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with the Arpra library. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef ARPRA_IMPL_H
#define ARPRA_IMPL_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif // HAVE_CONFIG_H

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <math.h>

#include <arpra.h>
#include <arpra_ode.h>

// Default precisions.
#define ARPRA_DEFAULT_PRECISION 53
#define ARPRA_DEFAULT_INTERNAL_PRECISION 256

// Temp buffers.
#define ARPRA_BUFFER_RESIZE_FACTOR 256

// Use Min-Range linear approximation
#define ARPRA_MIN_RANGE 1

/*
 * If ARPRA_TIGHT_MUL is defined, then the linear approximation of the quadratic term of
 * arpra_mul (in mul.c) is defined the same as in (26) of:
 *
 * S. M. Rump and M. Kashiwagi, Implementation and improvements of affine arithmetic,
 * Nonlinear Theory an Its Applications, IEICE, vol. 6, no. 3, pp. 341-359, 2015.
 *
 * Otherwise it is trivially defined as the product of the radii of x and y:
 *
 * \sum^{n}_{i=1} x_{i} \sum^{n}_{i=1} y_{i}
 */

#define ARPRA_TIGHT_MUL 1

// Internal helper functions.
void arpra_helper_error_ulp (arpra_mpfr *error, const arpra_mpfr *x);
void arpra_helper_error_half_ulp (arpra_mpfr *error, const arpra_mpfr *x);
int arpra_helper_term (arpra_mpfr *z, const arpra_mpfr *x, const arpra_mpfr *y,
                       const arpra_mpfr *alpha, const arpra_mpfr *beta,
                       const arpra_mpfr *gamma);
int arpra_helper_mpfr_sum (arpra_mpfr *z, arpra_mpfr *x,
                           const arpra_uint n, const mpfr_rnd_t rnd);
int arpra_helper_mpfr_sumabs (arpra_mpfr *z, arpra_mpfr *x,
                              const arpra_uint n, const mpfr_rnd_t rnd);
arpra_mpfr **arpra_helper_buffer_mpfr_ptr (arpra_uint n);
arpra_mpfr *arpra_helper_buffer_mpfr (arpra_uint n);

#endif // ARPRA_IMPL_H
