/*
 * rand_mpfr.c -- Generate a non-uniform random MPFR number.
 *
 * Copyright 2017 James Paul Turner.
 *
 * This file is part of the ArPRA library.
 *
 * The ArPRA library is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * The ArPRA library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with the ArPRA library. If not, see <http://www.gnu.org/licenses/>.
 */

#include "arpra-test.h"

void test_rand_mpfr (mpfr_ptr z, enum test_rand_mode mode)
{
    arpra_uint_t r_ui;
    mpfr_t r_mpfr;

    if (test_rand_ready) {
        // Check RNG mode.
        switch (mode) {
        case TEST_RAND_MIXED:
            r_ui = gmp_urandomm_ui(test_randstate, 4);
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
            r_ui = gmp_urandomm_ui(test_randstate, 2);
            break; // (+0 <= |z| <  +1)

        case TEST_RAND_LARGE:
            r_ui = gmp_urandomm_ui(test_randstate, 2) + 2;
            break; // (+1 <= |z| < +oo)

        case TEST_RAND_POS:
            r_ui = gmp_urandomm_ui(test_randstate, 2) * 2;
            break; // (+0 <=  z  < +oo)

        case TEST_RAND_NEG:
            r_ui = gmp_urandomm_ui(test_randstate, 2) * 2 + 1;
            break; // (-oo <  z  <= -0)

        default:
            fprintf(stderr, "Error: unrecognised RNG mode.\n");
            exit(EXIT_FAILURE);
        }

        // Generate random number with 5-bit precision.
        mpfr_init2(r_mpfr, 5);
        mpfr_urandomb (r_mpfr, test_randstate);
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
