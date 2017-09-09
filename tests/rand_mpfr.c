/*
 * rand_mpfr.c -- Generate a non-uniform random MPFR number.
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

void test_rand_mpfr (mpfr_ptr z, enum test_rand_mode mode)
{
    mpfa_uint_t r_ui;
    gmp_randstate_t *rand;
    mpfr_t r_mpfr;

    // Check RNG mode.
    if (test_rand_is_init()) {
        switch (mode) {
        case TEST_RAND_MIXED:
            r_ui = test_rand_ui(4);
            break; // (-oo <  z  < +oo)

        case TEST_RAND_SMALL_POS:
            r_ui = 0;
            break; // (+0 <=  z  <  +1)

        case TEST_RAND_SMALL_NEG:
            r_ui = 1;
            break; // (-1  <  z  <= -0)

        case TEST_RAND_LARGE_POS:
            r_ui = 2;
            break; // (+1 <=  z  < +oo)

        case TEST_RAND_LARGE_NEG:
            r_ui = 3;
            break; // (-oo <  z  <= -1)

        case TEST_RAND_SMALL:
            r_ui = test_rand_ui(2);
            break; // (+0 <= |z| <  +1)

        case TEST_RAND_LARGE:
            r_ui = test_rand_ui(2) + 2;
            break; // (+1 <= |z| < +oo)

        case TEST_RAND_POS:
            r_ui = test_rand_ui(2) * 2;
            break; // (+0 <=  z  < +oo)

        case TEST_RAND_NEG:
            r_ui = test_rand_ui(2) * 2 + 1;
            break; // (-oo <  z  <= -0)

        default:
            fprintf(stderr, "Error: unrecognised RNG mode.\n");
            exit(EXIT_FAILURE);
        }

        // Generate random number with 5-bit precision.
        rand = test_randstate_get();
        mpfr_init2(r_mpfr, 5);
        mpfr_urandomb (r_mpfr, *rand);
        if (r_ui == 1) {
            mpfr_neg (r_mpfr, r_mpfr, MPFR_RNDD);
        }
        else if (r_ui >= 2) {
            mpfr_ui_div (r_mpfr, 1, r_mpfr, MPFR_RNDD);
            if (r_ui == 3) {
                mpfr_neg (r_mpfr, r_mpfr, MPFR_RNDD);
            }
        }
        mpfr_set(z, r_mpfr, MPFR_RNDN);
        mpfr_clear(r_mpfr);
    }
    else {
        fprintf(stderr, "Error: RNG is not initialised.\n");
        exit(EXIT_FAILURE);
    }
}