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
    arpra_mpfr *lo_sum, **lo_sum_ptr;
    arpra_mpfr *hi_sum, **hi_sum_ptr;

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

    // Allocate 5 to 9 terms.
    arpra_clear_terms(z);
    z->nTerms = gmp_urandomm_ui(test_randstate, 5) + 5;
    z->symbols = malloc((z->nTerms + 1) * sizeof(arpra_uint));
    z->deviations = malloc((z->nTerms + 1) * sizeof(arpra_mpfr));
    lo_sum = malloc(z->nTerms * sizeof(arpra_mpfr));
    hi_sum = malloc(z->nTerms * sizeof(arpra_mpfr));
    lo_sum_ptr = malloc((z->nTerms + 1) * sizeof(arpra_mpfr *));
    hi_sum_ptr = malloc((z->nTerms + 1) * sizeof(arpra_mpfr *));

    // Set random deviation terms.
    for (zTerm = 0; zTerm < z->nTerms; zTerm++) {
        z->symbols[zTerm] = arpra_next_symbol();
        mpfr_init2(&(z->deviations[zTerm]), prec_internal);
        test_rand_mpfr(&(z->deviations[zTerm]), z->precision, mode_d);
        mpfr_abs(&temp1, &(z->deviations[zTerm]), MPFR_RNDU);
        mpfr_add(&(z->radius), &(z->radius), &temp1, MPFR_RNDU);

        // Set up true_range bound sums.
        lo_sum[zTerm] = z->deviations[zTerm];
        lo_sum[zTerm]._mpfr_sign = -1;
        lo_sum_ptr[zTerm] = &(lo_sum[zTerm]);
        hi_sum[zTerm] = z->deviations[zTerm];
        hi_sum[zTerm]._mpfr_sign = 1;
        hi_sum_ptr[zTerm] = &(hi_sum[zTerm]);
    }

    // Set up true_range bound sums.
    lo_sum_ptr[zTerm] = &(z->centre);
    hi_sum_ptr[zTerm] = &(z->centre);

    // Round range to target precision.
    if (mpfr_sum(&(z->true_range.left), lo_sum_ptr, (z->nTerms + 1), MPFR_RNDD)) {
        arpra_helper_error_ulp(&temp1, &(z->true_range.left));
    }
    else {
        mpfr_set_ui(&temp1, 0, MPFR_RNDN);
    }
    if (mpfr_sum(&(z->true_range.right), hi_sum_ptr, (z->nTerms + 1), MPFR_RNDU)) {
        arpra_helper_error_ulp(&temp2, &(z->true_range.right));
    }
    else {
        mpfr_set_ui(&temp2, 0, MPFR_RNDN);
    }
    mpfr_max(&error, &temp1, &temp2, MPFR_RNDU);
    mpfr_add(&(z->radius), &(z->radius), &error, MPFR_RNDU);

    // Store nonzero rounding error term.
    if (!mpfr_zero_p(&error) && mpfr_number_p(&error)) {
        z->symbols[zTerm] = arpra_next_symbol();
        z->deviations[zTerm] = error;
        z->nTerms++;
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
    free(lo_sum);
    free(hi_sum);
    free(lo_sum_ptr);
    free(hi_sum_ptr);
}
