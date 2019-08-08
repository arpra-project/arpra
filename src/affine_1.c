/*
 * affine_1.c -- Compute a univariate affine function of an Arpra range.
 *
 * Copyright 2016-2019 James Paul Turner.
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

void arpra_affine_1 (arpra_range *y, const arpra_range *x1,
                     mpfr_srcptr alpha, mpfr_srcptr gamma, mpfr_srcptr delta)
{
    arpra_uint i;
    arpra_range yy;
    mpfr_t temp, error;
    arpra_prec prec_internal;

    // Initialise vars.
    prec_internal = arpra_get_internal_precision();
    arpra_init2(&yy, y->precision);
    mpfr_init2(temp, prec_internal);
    mpfr_init2(error, prec_internal);
    mpfr_set_zero(error, 1);
    mpfr_set_zero(&(yy.radius), 1);

    // y[0] = (alpha * x1[0]) + (gamma)
    arpra_helper_mpfr_f3(error, &mpfr_fma, &(yy.centre), alpha, &(x1->centre), gamma, MPFR_RNDN);

    // Allocate memory for deviation terms.
    yy.symbols = malloc((x1->nTerms + 1) * sizeof(arpra_uint));
    yy.deviations = malloc((x1->nTerms + 1) * sizeof(arpra_mpfr));

    for (i = 0; i < x1->nTerms; i++) {
        mpfr_init2(&(yy.deviations[i]), prec_internal);

        // y[i] = (alpha * x1[i])
        yy.symbols[i] = x1->symbols[i];
        arpra_helper_mpfr_f2(error, &mpfr_mul, &(yy.deviations[i]), alpha, &(x1->deviations[i]), MPFR_RNDN);

        // Add term to radius.
        mpfr_abs(temp, &(yy.deviations[i]), MPFR_RNDU);
        mpfr_add(&(yy.radius), &(yy.radius), temp, MPFR_RNDU);
    }

    // Add delta to error.
    mpfr_add(error, error, delta, MPFR_RNDU);

    // Store numerical error term.
    yy.symbols[i] = arpra_helper_next_symbol();
    yy.deviations[i] = *error;
    mpfr_add(&(yy.radius), &(yy.radius), &(yy.deviations[i]), MPFR_RNDU);
    yy.nTerms = i + 1;

    // Clear vars, and set y.
    mpfr_clear(temp);
    arpra_clear(y);
    *y = yy;
}
