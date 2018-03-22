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
#ifdef WITH_MPFI
#include <arpra2mpfi.h>
#endif // WITH_MPFI

// Internal precision default value.
#define ARPRA_DEFAULT_INTERNAL_PRECISION 128

// Use tighter arpra_mul approximation.
#define ARPRA_TIGHT_MUL

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

// Internal helper functions.
arpra_int arpra_term (mpfr_ptr z, mpfr_srcptr x, mpfr_srcptr y,
                        mpfr_srcptr alpha, mpfr_srcptr beta, mpfr_srcptr gamma);
void arpra_error (mpfr_ptr error, mpfr_srcptr x);

#endif // ARPRA_IMPL_H
