/*
 * t_add.c -- Test the MPFA add function.
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

int main (int argc, char *argv[])
{
    const mpfa_uint_t n_tests = 100000;
    const mpfa_prec_t prec = 53;
    const mpfa_prec_t prec_internal = 128;

    mpfa_uint_t i, total_fail, t1_fail, t2_fail;
    mpfa_t a, b, c;
    mpfr_t a_lo, a_hi;
    mpfr_t b_lo, b_hi;
    mpfr_t c_lo, c_hi;

    // Init test.
    test_rand_init();
    mpfa_set_internal_prec(prec_internal);
    mpfa_init2(a, prec);
    mpfa_init2(b, prec);
    mpfa_init2(c, prec);
    mpfr_init2(a_lo, prec);
    mpfr_init2(a_hi, prec);
    mpfr_init2(b_lo, prec);
    mpfr_init2(b_hi, prec);
    mpfr_init2(c_lo, prec);
    mpfr_init2(c_hi, prec);
    total_fail = 0;

    for (t1_fail = 0, i = 0; i < n_tests; i++) {
        // Set random A.
        test_rand_mpfa(a, TEST_RAND_MIXED);
        mpfa_get_bounds(a_lo, a_hi, a);

        // Set random B.
        test_rand_mpfa(b, TEST_RAND_MIXED);
        mpfa_get_bounds(b_lo, b_hi, b);

        // Randomly share symbols.
        //test_share_syms(a, b, 5);

        // Compare MPFA result with MPFR bound results.
        mpfa_add(c, a, b);
        mpfr_add(c_lo, a_lo, b_lo, MPFR_RNDD);
        mpfr_add(c_hi, a_hi, b_hi, MPFR_RNDU);
        t1_fail += test_compare_bounds(c, c_lo, c_hi);
    }

    total_fail += t1_fail;
    printf("Test one: %llu out of %llu failed.\n", t1_fail, n_tests);

#ifdef WITH_MPFI
    // Compare with the MPFI interval arithmetic library.
    mpfi_t ai, bi, ci;

    // Init test vars, and set working precision.
    mpfi_init2(ai, prec);
    mpfi_init2(bi, prec);
    mpfi_init2(ci, prec);

    for (t2_fail = 0, i = 0; i < n_tests; i++) {
        // Set random A.
        test_rand_mpfa(a, TEST_RAND_MIXED);
        mpfa_get_mpfi(ai, a);
        mpfa_set_mpfi(a, ai);

        // Set random B.
        test_rand_mpfa(b, TEST_RAND_MIXED);
        mpfa_get_mpfi(bi, b);
        mpfa_set_mpfi(b, bi);

        // Randomly share symbols.
        //test_share_syms(a, b, 5);

        // Compare MPFA result with MPFI result.
        mpfa_add(c, a, b);
        mpfi_add(ci, ai, bi);
        t2_fail += test_compare_mpfi(c, ci);
    }

    total_fail += t2_fail;
    printf("Test two: %llu out of %llu failed.\n", t2_fail, n_tests);

    // Cleanup test vars.
    mpfi_clear(ai);
    mpfi_clear(bi);
    mpfi_clear(ci);

#else // WITH_MPFI
    fprintf(stderr,
            "Part of this test uses the MPFI interval arithmetic library.\n"
            "Recompile with MPFI support enabled to run the full test.\n");
#endif // WITH_MPFI

    // Cleanup test.
    mpfa_clear(a);
    mpfa_clear(b);
    mpfa_clear(c);
    mpfr_clear(a_lo);
    mpfr_clear(a_hi);
    mpfr_clear(b_lo);
    mpfr_clear(b_hi);
    mpfr_clear(c_lo);
    mpfr_clear(c_hi);
    test_rand_clear();
    mpfr_free_cache();
    return total_fail > 0;
}
