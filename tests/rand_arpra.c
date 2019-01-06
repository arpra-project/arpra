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

void test_rand_arpra (arpra_range *z, test_rand_mode mode_c, test_rand_mode mode_d)
{
    arpra_uint zTerm;
    arpra_prec prec, prec_internal;
    arpra_mpfr temp1, temp2, error;

    // Initialise vars.
    prec_internal = arpra_get_internal_precision();
    mpfr_init2(&temp1, prec_internal);
    mpfr_init2(&temp2, prec_internal);
    mpfr_init2(&error, prec_internal);
    mpfr_set_prec(&(z->centre), prec_internal);
    mpfr_set_prec(&(z->radius), prec_internal);
    mpfr_set_ui(&error, 0, MPFR_RNDU);
    mpfr_set_ui(&(z->radius), 0, MPFR_RNDU);

    // Set random centre.
    test_rand_mpfr(&(z->centre), z->precision, mode_c);

    // Allocate 0 to 9 terms.
    arpra_clear_terms(z);
    z->nTerms = gmp_urandomm_ui(test_randstate, 10);
    z->symbols = malloc((z->nTerms + 1) * sizeof(arpra_uint));
    z->deviations = malloc((z->nTerms + 1) * sizeof(arpra_mpfr));

    for (zTerm = 0; zTerm < z->nTerms; zTerm++) {
        z->symbols[zTerm] = arpra_next_symbol();
        mpfr_init2(&(z->deviations[zTerm]), prec_internal);
        test_rand_mpfr(&(z->deviations[zTerm]), z->precision, mode_d);
        mpfr_abs(&temp1, &(z->deviations[zTerm]), MPFR_RNDU);
        mpfr_add(&(z->radius), &(z->radius), &temp1, MPFR_RNDU);
    }

    // Round range to target precision.
    mpfr_sub(&temp1, &(z->centre), &(z->radius), MPFR_RNDD);
    mpfr_set(&(z->true_range.left), &temp1, MPFR_RNDD);
    mpfr_sub(&temp1, &temp1, &(z->true_range.left), MPFR_RNDU);
    mpfr_add(&temp2, &(z->centre), &(z->radius), MPFR_RNDU);
    mpfr_set(&(z->true_range.right), &temp2, MPFR_RNDU);
    mpfr_sub(&temp2, &(z->true_range.right), &temp2, MPFR_RNDU);
    mpfr_max(&temp1, &temp1, &temp2, MPFR_RNDU);
    mpfr_add(&error, &error, &temp1, MPFR_RNDU);

    // Store nonzero rounding error term.
    if (!mpfr_zero_p(&error) && mpfr_number_p(&error)) {
        z->nTerms++;
        z->symbols[zTerm] = arpra_next_symbol();
        z->deviations[zTerm] = error;
        mpfr_add(&(z->radius), &(z->radius), &(z->deviations[zTerm]), MPFR_RNDU);
    }
    else {
        mpfr_clear(&error);
    }

    // Clear unused memory
    if (z->nTerms == 0) {
        free(z->symbols);
        free(z->deviations);
    }

    // Handle domain violations.
    if (mpfr_nan_p(&(z->centre)) || mpfr_nan_p(&(z->radius))) {
        arpra_set_nan(z);
    }
    else if (mpfr_inf_p(&(z->centre)) || mpfr_inf_p(&(z->radius))) {
        arpra_set_inf(z);
    }

    // Clear vars.
    mpfr_clear(&temp1);
    mpfr_clear(&temp2);
}
