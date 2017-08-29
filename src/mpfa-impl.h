/*
 * mpfa.h -- MPFA private header file, used in compilation.
 *
 * Copyright 2016-2017 James Paul Turner.
 *
 * This file is part of the MPFA library.
 *
 * The MPFA library is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * The MPFA library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with the MPFA library. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef MPFA_IMPL_H
#define MPFA_IMPL_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif // HAVE_CONFIG_H

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>

#include <mpfa.h>
#ifdef WITH_MPFI
#include <mpfa2mpfi.h>
#endif // WITH_MPFI

#define MPFA_DEFAULT_INTERNAL_PREC 128

/*
 * If MPFA_TIGHT_MUL is defined, then the linear approximation of the quadratic term of
 * mpfa_mul (in mul.c) is defined the same as in (26) of:
 *
 * S. M. Rump and M. Kashiwagi, Implementation and improvements of affine arithmetic,
 * Nonlinear Theory an Its Applications, IEICE, vol. 6, no. 3, pp. 341-359, 2015.
 *
 * Otherwise it is trivially defined as the product of the radii of x and y:
 *
 * \sum^{n}_{i=1} x_{i} \sum^{n}_{i=1} y_{i}
 */

#define MPFA_TIGHT_MUL

#endif // MPFA_IMPL_H
