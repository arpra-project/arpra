/*
 * affine_2.c -- Compute a bivariate affine function of two Arpra ranges.
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

void arpra_affine_2 (arpra_range *y, const arpra_range *x1, const arpra_range *x2,
                     mpfr_srcptr a, mpfr_srcptr b, mpfr_srcptr c, mpfr_srcptr d)
{
    arpra_uint i, ix1, ix2;
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

    // y[0] = (a * x1[0]) + (b * x2[0]) + (c)
    arpra_helper_mpfr_fmmaa(&(yy.centre), &(x1->centre), &(x2->centre), a, b, c, error);

    // Allocate memory for deviation terms.
    yy.symbols = malloc((x1->nTerms + x2->nTerms + 1) * sizeof(arpra_uint));
    yy.deviations = malloc((x1->nTerms + x2->nTerms + 1) * sizeof(arpra_mpfr));

    for (i = 0, ix1 = 0, ix2 = 0; (ix1 < x1->nTerms) || (ix2 < x2->nTerms); i++) {
        mpfr_init2(&(yy.deviations[i]), prec_internal);

        if ((ix2 == x2->nTerms) || ((ix1 < x1->nTerms) && (x1->symbols[ix1] < x2->symbols[ix2]))) {
            // y[i] = (a * x1[i])
            yy.symbols[i] = x1->symbols[ix1];
            arpra_helper_mpfr_f2(&mpfr_mul, &(yy.deviations[i]), a, &(x1->deviations[ix1]), MPFR_RNDN, error);
            ix1++;
        }
        else if ((ix1 == x1->nTerms) || ((ix2 < x2->nTerms) && (x2->symbols[ix2] < x1->symbols[ix1]))) {
            // y[i] = (b * x2[i])
            yy.symbols[i] = x2->symbols[ix2];
            arpra_helper_mpfr_f2(&mpfr_mul, &(yy.deviations[i]), b, &(x2->deviations[ix2]), MPFR_RNDN, error);
            ix2++;
        }
        else {
            // y[i] = (a * x1[i]) + (b * x2[i])
            yy.symbols[i] = x1->symbols[ix1];
            arpra_helper_mpfr_fmma(&(yy.deviations[i]), &(x1->deviations[ix1]), &(x2->deviations[ix2]), a, b, error);
            ix1++;
            ix2++;
        }

        // Add term to radius.
        mpfr_abs(temp, &(yy.deviations[i]), MPFR_RNDU);
        mpfr_add(&(yy.radius), &(yy.radius), temp, MPFR_RNDU);
    }

    // Add d to error.
    mpfr_add(error, error, d, MPFR_RNDU);

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
