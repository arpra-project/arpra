/*
 * rand_arpra.c -- Generate a random arpra variable.
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

#include "arpra-test.h"

void test_rand_arpra (arpra_range *z,
                      enum test_rand_mode mode_centre,
                      enum test_rand_mode mode_deviations)
{
    arpra_uint zTerm;
    arpra_prec prec, prec_internal;
    arpra_mpfr temp;

    // Initialise vars.
    prec = arpra_get_precision(z);
    prec_internal = arpra_get_internal_precision();
    mpfr_init2(&temp, prec_internal);
    mpfr_set_prec(&(z->radius), prec_internal);
    mpfr_set_ui(&(z->radius), 0, MPFR_RNDN);

    // Set random centre.
    test_rand_mpfr(&(z->centre), mode_centre);

    // Clear existing deviation terms.
    arpra_clear_terms(z);

    // Randomly allocate 0 to 5 deviation terms.
    z->nTerms = gmp_urandomm_ui(test_randstate, 6);
    if (z->nTerms > 0) {
        z->symbols = malloc(z->nTerms * sizeof(arpra_uint));
        z->deviations = malloc(z->nTerms * sizeof(arpra_mpfr));
    }

    // Randomly set deviation terms.
    for (zTerm = 0; zTerm < z->nTerms; zTerm++) {
        z->symbols[zTerm] = arpra_next_symbol();
        mpfr_init2(&(z->deviations[zTerm]), prec);
        test_rand_mpfr(&(z->deviations[zTerm]), mode_deviations);

        // Add abs(term) to radius.
        mpfr_abs(&temp, &(z->deviations[zTerm]), MPFR_RNDU);
        mpfr_add(&(z->radius), &(z->radius), &temp, MPFR_RNDU);
    }

    // Handle domain violations.
    if (mpfr_nan_p(&(z->centre)) || mpfr_nan_p(&(z->radius))) {
        arpra_set_nan(z);
    }
    else if (mpfr_inf_p(&(z->centre)) || mpfr_inf_p(&(z->radius))) {
        arpra_set_inf(z);
    }

    // Clear vars.
    mpfr_clear(&temp);
}
