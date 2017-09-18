/*
 * set.c -- Set one affine form with the values of another.
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

void mpfa_set (mpfa_ptr z, mpfa_srcptr x)
{
    mpfa_uint_t xTerm, zTerm;
    mpfa_prec_t prec, prec_internal;
    mpfr_t temp, error;

    // Handle trivial cases.
    if (z == x) return;

    // Handle domain violations.
    if (mpfa_none_p(x)) {
        mpfa_set_none(z);
        return;
    }
    if (mpfa_any_p(x)) {
        mpfa_set_any(z);
        return;
    }

    // Init temp vars, and set internal precision.
    prec = mpfa_get_prec(z);
    prec_internal = mpfa_get_internal_prec();
    mpfr_init2(temp, prec_internal);
    mpfr_init2(error, prec_internal);
    mpfr_prec_round(&(z->radius), prec_internal, MPFR_RNDU);
    mpfr_set_si(error, 0, MPFR_RNDU);
    mpfr_set_si(&(z->radius), 0, MPFR_RNDU);

    // z_0 = x_0
    if (mpfr_set(&(z->centre), &(x->centre), MPFR_RNDN)) {
        mpfa_error(temp, &(z->centre));
        mpfr_add(error, error, temp, MPFR_RNDU);
    }

    // Replace existing noise term memory.
    mpfa_clear_terms(z);
    z->nTerms = x->nTerms + 1;
    z->symbols = malloc(z->nTerms * sizeof(mpfa_uint_t));
    z->deviations = malloc(z->nTerms * sizeof(mpfa_t));

    // Copy noise terms over.
    for (xTerm = 0, zTerm = 0; xTerm < x->nTerms; xTerm++) {
        z->symbols[zTerm] = x->symbols[zTerm];
        mpfr_init2(&(z->deviations[zTerm]), prec);

        // z_i = x_i
        if (mpfr_set(&(z->deviations[zTerm]), &(x->deviations[xTerm]), MPFR_RNDN)) {
            mpfa_error(temp, &(z->deviations[zTerm]));
            mpfr_add(error, error, temp, MPFR_RNDU);
        }

        // Store nonzero noise terms.
        if (mpfr_zero_p(&(z->deviations[zTerm]))) {
            mpfr_clear(&(z->deviations[zTerm]));
        }
        else {
            mpfr_abs(temp, &(z->deviations[zTerm]), MPFR_RNDU);
            mpfr_add(&(z->radius), &(z->radius), temp, MPFR_RNDU);
            zTerm++;
        }
    }

    // Store nonzero numerical error term.
    if (!mpfr_zero_p(error)) {
        z->symbols[zTerm] = mpfa_next_sym();
        mpfr_init2(&(z->deviations[zTerm]), prec);
        mpfr_set(&(z->deviations[zTerm]), error, MPFR_RNDU);
        mpfr_add(&(z->radius), &(z->radius), &(z->deviations[zTerm]), MPFR_RNDU);
        zTerm++;
    }

    // Handle domain violations, and resize memory.
    z->nTerms = zTerm;
    if (mpfr_nan_p(&(z->centre)) || mpfr_nan_p(&(z->radius))) {
        mpfa_set_none(z);
    }
    else if (mpfr_inf_p(&(z->centre)) || mpfr_inf_p(&(z->radius))) {
        mpfa_set_any(z);
    }
    else {
        if (z->nTerms == 0) {
            free(z->symbols);
            free(z->deviations);
        }
        else {
            z->symbols = realloc(z->symbols, z->nTerms * sizeof(mpfa_uint_t));
            z->deviations = realloc(z->deviations, z->nTerms * sizeof(mpfr_t));
        }
    }

    // Clear temp vars.
    mpfr_clear(temp);
    mpfr_clear(error);
}
