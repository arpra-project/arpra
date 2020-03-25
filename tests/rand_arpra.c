/*
 * rand_arpra.c -- Generate a random arpra range.
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
    mpfr_t error;
    arpra_range yy;
    arpra_prec prec_internal;
    arpra_uint iy;

    // Initialise vars.
    prec_internal = arpra_get_internal_precision();
    mpfr_init2(error, prec_internal);
    arpra_init2(&yy, y->precision);
    mpfr_set_zero(error, 1);

    // y[0] = rand()
    test_rand_mpfr(&(yy.centre), prec_internal, mode_c);

    // Allocate 0 to 9 terms.
    yy.nTerms = gmp_urandomm_ui(test_randstate, 10);
    yy.symbols = malloc((yy.nTerms + 1) * sizeof(arpra_uint));
    yy.deviations = malloc((yy.nTerms + 1) * sizeof(mpfr_t));

    for (iy = 0; iy < yy.nTerms; iy++) {
        mpfr_init2(&(yy.deviations[iy]), prec_internal);

        // y[i] = rand()
        yy.symbols[iy] = arpra_helper_next_symbol();
        test_rand_mpfr(&(yy.deviations[iy]), prec_internal, mode_d);
    }

    // Store new deviation term.
    yy.symbols[iy] = arpra_helper_next_symbol();
    yy.deviations[iy] = *error;
    yy.nTerms = iy + 1;

    // Compute true_range.
    arpra_helper_compute_range(&yy);

    // Check for NaN and Inf.
    arpra_helper_check_result(&yy);

    // Clear vars.
    arpra_clear(y);
    *y = yy;
}

void test_rand_uniform_arpra (arpra_range *y,
                              long int yc_a, long int yc_b,
                              long int yd_a, long int yd_b)
{
    mpfr_t error;
    arpra_range yy;
    arpra_prec prec_internal;
    arpra_uint iy;

    // Initialise vars.
    prec_internal = arpra_get_internal_precision();
    mpfr_init2(error, prec_internal);
    arpra_init2(&yy, y->precision);
    mpfr_set_zero(error, 1);

    // y[0] = rand()
    test_rand_uniform_mpfr(&(yy.centre), yc_a, yc_b);

    // Allocate 0 to 9 terms.
    yy.nTerms = gmp_urandomm_ui(test_randstate, 10);
    yy.symbols = malloc((yy.nTerms + 1) * sizeof(arpra_uint));
    yy.deviations = malloc((yy.nTerms + 1) * sizeof(mpfr_t));

    for (iy = 0; iy < yy.nTerms; iy++) {
        mpfr_init2(&(yy.deviations[iy]), prec_internal);

        // y[i] = rand()
        yy.symbols[iy] = arpra_helper_next_symbol();
        test_rand_uniform_mpfr(&(yy.deviations[iy]), yd_a, yd_b);
    }

    // Store new deviation term.
    yy.symbols[iy] = arpra_helper_next_symbol();
    yy.deviations[iy] = *error;
    yy.nTerms = iy + 1;

    // Compute true_range.
    arpra_helper_compute_range(&yy);

    // Check for NaN and Inf.
    arpra_helper_check_result(&yy);

    // Clear vars.
    arpra_clear(y);
    *y = yy;
}
