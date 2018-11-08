/*
 * set_mpfi.c -- Set an Arpra range with an MPFI interval.
 *
 * Copyright 2017-2018 James Paul Turner.
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

void arpra_set_mpfi (arpra_range *z, mpfi_srcptr x)
{
    arpra_prec prec, prec_internal;
    arpra_mpfr temp;

    // Handle domain violations.
    if (mpfi_nan_p(x)) {
        arpra_set_nan(z);
        return;
    }
    if (mpfi_inf_p(x)) {
        arpra_set_inf(z);
        return;
    }

    // Initialise vars.
    prec = arpra_get_precision(z);
    prec_internal = arpra_get_internal_precision();
    mpfr_init2(&temp, prec_internal);
    mpfr_set_prec(&(z->radius), prec_internal);

    // z_0 = (x_lo + x_hi) / 2
    mpfr_add(&(z->centre), &(x->left), &(x->right), MPFR_RNDN);
    mpfr_div_2ui(&(z->centre), &(z->centre), 1, MPFR_RNDN);

    // rad(z) = max{(z_0 - x_lo), (x_hi - z_0)}
    mpfr_sub(&(z->radius), &(z->centre), &(x->left), MPFR_RNDU);
    mpfr_sub(&temp, &(x->right), &(z->centre), MPFR_RNDU);
    mpfr_max(&(z->radius), &(z->radius), &temp, MPFR_RNDU);

    // Clear existing deviation terms.
    arpra_clear_terms(z);

    // If radius is nonzero:
    if (!mpfr_zero_p(&(z->radius))) {
        // Allocate one deviation term.
        z->nTerms = 1;
        z->symbols = malloc(sizeof(arpra_uint));
        z->deviations = malloc(sizeof(arpra_range));

        // Set deviation term.
        z->symbols[0] = arpra_next_symbol();
        mpfr_init2(&(z->deviations[0]), prec);
        mpfr_set(&(z->deviations[0]), &(z->radius), MPFR_RNDU);
    }

    // Handle domain violations.
    if (mpfr_nan_p(&(z->centre)) || mpfr_nan_p(&(z->deviations[0]))) {
        arpra_set_nan(z);
    }
    else if (mpfr_inf_p(&(z->centre)) || mpfr_inf_p(&(z->deviations[0]))) {
        arpra_set_inf(z);
    }

    // Clear vars.
    mpfr_clear(&temp);
}
