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
    mpfa_uint_t r;
    gmp_randstate_t *rand;

    // Check RNG mode.
    if (test_rand_is_init()) {
        switch (mode) {
        case TEST_RAND_MIXED:
            r = test_rand_ui(4);
            break;

        case TEST_RAND_SMALL_POS:
            r = 0;
            break;

        case TEST_RAND_SMALL_NEG:
            r = 1;
            break;

        case TEST_RAND_LARGE_POS:
            r = 2;
            break;

        case TEST_RAND_LARGE_NEG:
            r = 3;
            break;

        case TEST_RAND_SMALL:
            r = test_rand_ui(2);
            break;

        case TEST_RAND_LARGE:
            r = test_rand_ui(2) + 2;
            break;

        case TEST_RAND_POS:
            r = test_rand_ui(2) * 2;
            break;

        case TEST_RAND_NEG:
            r = test_rand_ui(2) * 2 + 1;
            break;

        default:
            fprintf(stderr, "Error: unrecognised RNG mode.\n");
            exit(EXIT_FAILURE);
        }

        // Generate number.
        rand = test_rand_get();
        mpfr_urandomb (z, *rand);
        if (r == 1) {
            mpfr_neg (z, z, MPFR_RNDD);
        }
        else if (r >= 2) {
            mpfr_ui_div (z, 1, z, MPFR_RNDD);
            if (r == 3) mpfr_neg (z, z, MPFR_RNDD);
        }
    }
    else {
        fprintf(stderr, "Error: RNG is not initialised.\n");
        exit(EXIT_FAILURE);
    }
}
