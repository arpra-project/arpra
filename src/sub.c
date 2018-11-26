/*
 * sub.c -- Subtract one Arpra range from another.
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

#include "arpra-impl.h"

void arpra_sub (arpra_range *z, const arpra_range *x, const arpra_range *y)
{
    arpra_mpfr alpha, beta, gamma, delta;
    arpra_mpfi z_range;

    // Initialise vars.
    mpfr_init2(&alpha, 2);
    mpfr_init2(&beta, 2);
    mpfr_init2(&gamma, 2);
    mpfr_init2(&delta, 2);
    mpfi_init2(&z_range, z->precision);
    mpfr_set_si(&alpha, 1, MPFR_RNDN);
    mpfr_set_si(&beta, -1, MPFR_RNDN);
    mpfr_set_si(&gamma, 0, MPFR_RNDN);
    mpfr_set_si(&delta, 0, MPFR_RNDN);

    // MPFI subtraction
    mpfi_sub(&z_range, &(x->true_range), &(y->true_range));

    // z = x - y
    arpra_affine_2(z, x, y, &alpha, &beta, &gamma, &delta);

    // Compute true range.
    mpfi_intersect(&(z->true_range), &(z->true_range), &z_range);

    // Clear vars.
    mpfr_clear(&alpha);
    mpfr_clear(&beta);
    mpfr_clear(&gamma);
    mpfr_clear(&delta);
    mpfi_clear(&z_range);
}
