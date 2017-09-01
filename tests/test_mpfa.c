/*
 * test_mpfa.c -- Test an affine form against another affine form.
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

#include "mpfa-test.h"

mpfa_int_t mpfa_test_cmp_mpfa (mpfa_srcptr x, mpfa_srcptr y)
{
    mpfa_uint_t term;

    // Return 1 if x and y differ.
    if (!mpfr_equal_p(&(x->centre), &(y->centre))) {
        return 1;
    }

    if (!mpfr_equal_p(&(x->radius), &(y->radius))) {
        return 1;
    }

    if (x->nTerms != y->nTerms) {
        return 1;
    }

    for (term = 0; term < x->nTerms; term++) {
        if (x->symbols[term] != y->symbols[term]) {
            return 1;
        }

        if (!mpfr_equal_p(&(x->deviations[term]), &(y->deviations[term]))) {
            return 1;
        }
    }

    // Else return 0.
    return 0;
}

void mpfa_test_rand_mpfa (mpfa_ptr z, enum mpfa_test_rand_mode mode)
{
    mpfa_uint_t zTerm;
    mpfa_prec_t prec, prec_internal;
    mpfr_t temp;

    // Init temp vars and set internal precision.
    prec = mpfa_get_prec(z);
    prec_internal = mpfa_get_internal_prec();
    mpfr_init2(temp, prec_internal);
    mpfr_prec_round(&(z->radius), prec_internal, MPFR_RNDU);
    mpfr_set_ui(&(z->radius), 0, MPFR_RNDN);

    // Set random centre.
    mpfa_test_rand_mpfr(&(z->centre), mode);

    // Clear existing noise terms.
    mpfa_clear_terms(z);

    // Randomly allocate between 0 and 7 noise terms.
    mpfa_test_rand_ui(3);
    z->symbols = malloc(z->nTerms * sizeof(mpfa_uint_t));
    z->deviations = malloc(z->nTerms * sizeof(mpfr_t));

    // Randomly set noise terms.
    for (zTerm = 0; zTerm < z->nTerms; zTerm++) {
        z->symbols[zTerm] = mpfa_next_sym();
        mpfr_init2(&(z->deviations[zTerm]), prec_internal);
        mpfa_test_rand_mpfr(&(z->deviations[zTerm]), mode);

        // Add abs(term) to radius.
        mpfr_abs(temp, &(z->deviations[zTerm]), MPFR_RNDN);
        mpfr_add(&(z->radius), &(z->radius), temp, MPFR_RNDU);
    }

    // Clear temp vars.
    mpfr_clear(temp);
}
