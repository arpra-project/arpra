/*
 * rand_mpfr.c -- Generate a non-uniform random MPFR number.
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

void test_rand_mpfr (mpfr_ptr y, arpra_prec prec, test_rand_mode mode)
{
    arpra_uint r_ui;
    mpfr_t r_mpfr;

    if (test_rand_ready) {
        // Check RNG mode.
        switch (mode) {
        case TEST_RAND_MIXED:
            r_ui = gmp_urandomm_ui(test_randstate, 4);
            break; // (-oo <  y  < +oo)

        case TEST_RAND_SMALL_POS:
            r_ui = 0;
            break; // (+0 <=  y  <  +1)

        case TEST_RAND_SMALL_NEG:
            r_ui = 1;
            break; // (-1  <  y  <= -0)

        case TEST_RAND_LARGE_POS:
            r_ui = 2;
            break; // (+1 <=  y  < +oo)

        case TEST_RAND_LARGE_NEG:
            r_ui = 3;
            break; // (-oo <  y  <= -1)

        case TEST_RAND_SMALL:
            r_ui = gmp_urandomm_ui(test_randstate, 2);
            break; // (+0 <= |y| <  +1)

        case TEST_RAND_LARGE:
            r_ui = gmp_urandomm_ui(test_randstate, 2) + 2;
            break; // (+1 <= |y| < +oo)

        case TEST_RAND_POS:
            r_ui = gmp_urandomm_ui(test_randstate, 2) * 2;
            break; // (+0 <=  y  < +oo)

        case TEST_RAND_NEG:
            r_ui = gmp_urandomm_ui(test_randstate, 2) * 2 + 1;
            break; // (-oo <  y  <= -0)

        default:
            fprintf(stderr, "Error: unrecognised RNG mode.\n");
            exit(EXIT_FAILURE);
        }

        // Generate random number.
        mpfr_init2(r_mpfr, prec);
        mpfr_urandom(r_mpfr, test_randstate, MPFR_RNDN);
        if (r_ui == 1) {
            mpfr_neg(r_mpfr, r_mpfr, MPFR_RNDD);
        }
        else if (r_ui >= 2) {
            mpfr_ui_div(r_mpfr, 1, r_mpfr, MPFR_RNDD);
            if (r_ui == 3) {
                mpfr_neg(r_mpfr, r_mpfr, MPFR_RNDD);
            }
        }
        mpfr_set(y, r_mpfr, MPFR_RNDN);
        mpfr_clear(r_mpfr);
    }
    else {
        fprintf(stderr, "Error: RNG is not initialised.\n");
        exit(EXIT_FAILURE);
    }
}

void test_rand_uniform_mpfr (mpfr_ptr y, long int y_a, long int y_b)
{
    mpfr_t temp;

    // Check RNG readiness.
    if (!test_rand_ready) {
        fprintf(stderr, "Error: RNG is not initialised.\n");
        exit(EXIT_FAILURE);
    }

    // Generate random number.
    mpfr_urandom(y, test_randstate, MPFR_RNDN);
    mpfr_init2(temp, (mpfr_get_prec(y) * 2));
    mpfr_set_si(temp, (y_b - y_a), MPFR_RNDN);
    mpfr_mul(temp, temp, y, MPFR_RNDN);
    mpfr_add_si(y, temp, y_a, MPFR_RNDN);
    mpfr_clear(temp);
}
