/*
 * rand_arpra.c -- Generate a random arpra variable.
 *
 * Copyright 2017-2020 James Paul Turner.
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

void test_rand_arpra (arpra_range *y, test_rand_mode mode_c, test_rand_mode mode_d)
{
    mpfr_t temp1, temp2, error;
    arpra_range yy;
    arpra_prec prec_internal;
    arpra_uint iy;
    mpfr_ptr lo_sum, hi_sum;
    mpfr_ptr *lo_sum_ptr, *hi_sum_ptr;

    // Initialise vars.
    prec_internal = arpra_get_internal_precision();
    mpfr_init2(temp1, prec_internal);
    mpfr_init2(temp2, prec_internal);
    mpfr_init2(error, prec_internal);
    arpra_init2(&yy, y->precision);
    mpfr_set_zero(error, 1);
    mpfr_set_zero(&(yy.radius), 1);

    // y[0] = rand()
    test_rand_mpfr(&(yy.centre), yy.precision, mode_c);

    // Allocate 5 to 9 terms.
    yy.nTerms = gmp_urandomm_ui(test_randstate, 5) + 5;
    yy.symbols = malloc((yy.nTerms + 1) * sizeof(arpra_uint));
    yy.deviations = malloc((yy.nTerms + 1) * sizeof(mpfr_t));
    lo_sum = malloc(yy.nTerms * sizeof(mpfr_t));
    hi_sum = malloc(yy.nTerms * sizeof(mpfr_t));
    lo_sum_ptr = malloc((yy.nTerms + 1) * sizeof(mpfr_ptr));
    hi_sum_ptr = malloc((yy.nTerms + 1) * sizeof(mpfr_ptr));

    for (iy = 0; iy < yy.nTerms; iy++) {
        mpfr_init2(&(yy.deviations[iy]), prec_internal);

        // y[i] = rand()
        yy.symbols[iy] = arpra_helper_next_symbol();
        test_rand_mpfr(&(yy.deviations[iy]), yy.precision, mode_d);

        // Add term to radius.
        mpfr_abs(temp1, &(yy.deviations[iy]), MPFR_RNDU);
        mpfr_add(&(yy.radius), &(yy.radius), temp1, MPFR_RNDU);

        // Set up true_range bound sums.
        lo_sum[iy] = yy.deviations[iy];
        lo_sum[iy]._mpfr_sign = -1;
        lo_sum_ptr[iy] = &(lo_sum[iy]);
        hi_sum[iy] = yy.deviations[iy];
        hi_sum[iy]._mpfr_sign = 1;
        hi_sum_ptr[iy] = &(hi_sum[iy]);
    }

    // Set up true_range bound sums.
    lo_sum_ptr[iy] = &(yy.centre);
    hi_sum_ptr[iy] = &(yy.centre);

    // Compute true_range in working precision.
    mpfr_set_zero(temp1, -1);
    ARPRA_MPFR_RNDERR_SUM(temp1, MPFR_RNDD, &(yy.true_range.left), lo_sum_ptr, (iy + 1));
    mpfr_set_zero(temp2,  1);
    ARPRA_MPFR_RNDERR_SUM(temp2, MPFR_RNDU, &(yy.true_range.right), hi_sum_ptr, (iy + 1));
    mpfr_max(error, temp1, temp2, MPFR_RNDU);

    // Store new deviation term.
    yy.symbols[iy] = arpra_helper_next_symbol();
    yy.deviations[iy] = *error;
    mpfr_add(&(yy.radius), &(yy.radius), &(yy.deviations[iy]), MPFR_RNDU);
    yy.nTerms = iy + 1;

    // Check for NaN and Inf.
    arpra_helper_check_result(&yy);

    // Clear vars.
    mpfr_clear(temp1);
    mpfr_clear(temp2);
    arpra_clear(y);
    *y = yy;
    free(lo_sum);
    free(hi_sum);
    free(lo_sum_ptr);
    free(hi_sum_ptr);
}
