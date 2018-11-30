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
    arpra_mpfr temp1, temp2, *error;
    arpra_mpfi ia_range;
    arpra_prec prec_internal;

    // Initialise vars.
    prec_internal = arpra_get_internal_precision();
    mpfr_init2(&alpha, 2);
    mpfr_init2(&beta, 2);
    mpfr_init2(&gamma, 2);
    mpfr_init2(&delta, 2);
    mpfr_init2(&temp1, prec_internal);
    mpfr_init2(&temp2, prec_internal);
    mpfi_init2(&ia_range, z->precision);
    mpfr_set_si(&alpha, 1, MPFR_RNDN);
    mpfr_set_si(&beta, -1, MPFR_RNDN);
    mpfr_set_si(&gamma, 0, MPFR_RNDN);
    mpfr_set_si(&delta, 0, MPFR_RNDN);

    // MPFI subtraction
    mpfi_sub(&ia_range, &(x->true_range), &(y->true_range));

    // z = x - y
    arpra_affine_2(z, x, y, &alpha, &beta, &gamma, &delta);

    // Trim error term if Arpra range fully contains IA range.
    error = &(z->deviations[z->nTerms - 1]);
    if (mpfr_less_p(&(z->true_range.left), &(ia_range.left))
        && mpfr_greater_p(&(z->true_range.right), &(ia_range.right))) {
        mpfr_sub(&temp1, &(ia_range.left), &(z->true_range.left), MPFR_RNDD);
        mpfr_sub(&temp2, &(z->true_range.right), &(ia_range.right), MPFR_RNDD);
        mpfr_min(&temp1, &temp1, &temp2, MPFR_RNDD);
        mpfr_sub(error, error, &temp1, MPFR_RNDU);
        if (mpfr_cmp_ui(error, 0) < 0) {
            mpfr_set_ui(error, 0, MPFR_RNDN);
        }
    }
    mpfi_intersect(&(z->true_range), &(z->true_range), &ia_range);

    // Clear vars.
    mpfr_clear(&alpha);
    mpfr_clear(&beta);
    mpfr_clear(&gamma);
    mpfr_clear(&delta);
    mpfr_clear(&temp1);
    mpfr_clear(&temp2);
    mpfi_clear(&ia_range);
}
