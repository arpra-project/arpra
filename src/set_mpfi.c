/*
 * set_mpfi.c -- Set an affine form with an MPFI interval.
 *
 * Copyright 2017 James Paul Turner.
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

#include "mpfa-impl.h"

void mpfa_set_mpfi (mpfa_ptr z, mpfi_srcptr x)
{
    mpfa_prec_t prec, prec_internal;
    mpfr_t temp;

    // Init temp var.
    prec = mpfa_get_prec(z);
    prec_internal = mpfa_get_internal_prec();
    mpfr_init2(temp, prec_internal);

    // z_0 = (x_lo + x_hi) / 2
    mpfr_add(&(z->centre), &(x->left), &(x->right), MPFR_RNDN);
    mpfr_div_2ui(&(z->centre), &(z->centre), 1, MPFR_RNDN);

    // rad(z) = max{(z_0 - x_lo), (x_hi - z_0)}
    mpfr_sub(&(z->radius), &(z->centre), &(x->left), MPFR_RNDU);
    mpfr_sub(temp, &(x->right), &(z->centre), MPFR_RNDU);
    mpfr_max(&(z->radius), &(z->radius), temp, MPFR_RNDU);

    // Clear existing noise terms.
    mpfa_clear_terms(z);

    // If radius is nonzero:
    if (!mpfr_zero_p(&(z->radius))) {
        // Allocate one noise term.
        z->nTerms = 1;
        z->symbols = malloc(sizeof(mpfa_uint_t));
        z->deviations = malloc(sizeof(mpfa_t));

        // Set noise term.
        z->symbols[0] = mpfa_next_sym();
        mpfr_init2(&(z->deviations[0]), prec);
        mpfr_set(&(z->deviations[0]), &(z->radius), MPFR_RNDU);
    }
}
