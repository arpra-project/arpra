/*
 * affine_1.c -- Compute a univariate affine function of an affine form.
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

#include "mpfa-impl.h"

void mpfa_affine_1 (mpfa_ptr z, mpfa_srcptr x, mpfr_srcptr alpha, mpfr_srcptr gamma, mpfr_srcptr delta)
{
    mpfa_uint_t xTerm, zTerm;
    mpfr_t temp, error;
    mpfa_prec_t prec, prec_internal;
    mpfa_t zNew;

    // Init temp vars, and set internal precision.
    prec = mpfa_get_prec(z);
    prec_internal = mpfa_get_internal_prec();
    mpfr_init2(temp, prec_internal);
    mpfr_init2(error, prec_internal);
    mpfr_init2(&(zNew->centre), prec);
    mpfr_init2(&(zNew->radius), prec_internal);
    mpfr_set(error, delta, MPFR_RNDU);
    mpfr_set_si(&(zNew->radius), 0, MPFR_RNDU);

    // z_0 = (alpha * x_0) + gamma
    if (mpfr_fma(&(zNew->centre), alpha, &(x->centre), gamma, MPFR_RNDN)) {
        mpfa_error(temp, &(zNew->centre));
        mpfr_add(error, error, temp, MPFR_RNDU);
    }

    // Allocate memory for all possible noise terms in z.
    zNew->nTerms = x->nTerms + 1;
    zNew->symbols = malloc(zNew->nTerms * sizeof(mpfa_uint_t));
    zNew->deviations = malloc(zNew->nTerms * sizeof(mpfr_t));

    for (xTerm = 0, zTerm = 0; xTerm < x->nTerms; xTerm++) {
        zNew->symbols[zTerm] = x->symbols[xTerm];
        mpfr_init2(&(zNew->deviations[zTerm]), prec);

        // z_i = (alpha * x_i)
        if (mpfr_mul(&(zNew->deviations[zTerm]), alpha, &(x->deviations[xTerm]), MPFR_RNDN)) {
            mpfa_error(temp, &(zNew->deviations[zTerm]));
            mpfr_add(error, error, temp, MPFR_RNDU);
        }

        // Store nonzero noise terms.
        if (mpfr_zero_p(&(zNew->deviations[zTerm]))) {
            mpfr_clear(&(zNew->deviations[zTerm]));
        }
        else {
            mpfr_abs(temp, &(zNew->deviations[zTerm]), MPFR_RNDU);
            mpfr_add(&(zNew->radius), &(zNew->radius), temp, MPFR_RNDU);
            zTerm++;
        }
    }

    // Store nonzero numerical error term.
    if (!mpfr_zero_p(error)) {
        zNew->symbols[zTerm] = mpfa_next_sym();
        mpfr_init2(&(zNew->deviations[zTerm]), prec);
        mpfr_set(&(zNew->deviations[zTerm]), error, MPFR_RNDU);
        mpfr_add(&(zNew->radius), &(zNew->radius), &(zNew->deviations[zTerm]), MPFR_RNDU);
        zTerm++;
    }

    // Free noise term memory if number of terms is zero.
    zNew->nTerms = zTerm;
    if (zNew->nTerms == 0) {
        free(zNew->symbols);
        free(zNew->deviations);
    }

    // Clear temp vars, and set z.
    mpfr_clear(temp);
    mpfr_clear(error);
    mpfa_set(z, zNew);
    mpfa_clear(zNew);
}
