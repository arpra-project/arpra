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
#ifdef WITH_MPFI
    const mpfa_uint_t n_tests = 100000;
    const mpfa_prec_t prec = 53;
    const mpfa_prec_t prec_internal = 128;

    mpfa_t a_A, b_A, c_A;
    mpfi_t a_I, b_I, c_I, c_AI;
    mpfa_uint_t i, total_fail;
    mpfa_uint_t t1_fail, t2_fail;

    // Init test.
    test_rand_init();
    mpfa_set_internal_prec(prec_internal);
    mpfa_init2(a_A, prec);
    mpfa_init2(b_A, prec);
    mpfa_init2(c_A, prec);
    mpfi_init2(a_I, prec);
    mpfi_init2(b_I, prec);
    mpfi_init2(c_I, prec);
    mpfi_init2(c_AI, prec);
    total_fail = 0;

    for (t1_fail = 0, i = 0; i < n_tests; i++) {
        // Set random A.
        test_rand_mpfa(a_A, TEST_RAND_SMALL);
        mpfa_get_mpfi(a_I, a_A);

        // Set random B.
        test_rand_mpfa(b_A, TEST_RAND_SMALL);
        mpfa_get_mpfi(b_I, b_A);

        // Randomly share symbols.
        //test_share_syms(a_A, b_A, 5);

        // MPFA and MPFI operations.
        mpfa_add(c_A, a_A, b_A);
        mpfi_add(c_I, a_I, b_I);

        // Compare MPFA result with MPFI result.
        mpfa_get_mpfi(c_AI, c_A);
        t1_fail += !mpfi_is_inside(c_I, c_AI);
    }

    total_fail += t1_fail;
    printf("Test one: %llu out of %llu failed.\n", t1_fail, n_tests);

    // Cleanup test.
    mpfa_clear(a_A);
    mpfa_clear(b_A);
    mpfa_clear(c_A);
    mpfi_clear(a_I);
    mpfi_clear(b_I);
    mpfi_clear(c_I);
    mpfi_clear(c_AI);
    test_rand_clear();
    mpfr_free_cache();
    return total_fail > 0;

#else // WITH_MPFI
    fprintf(stderr,
            "This test uses the MPFI interval arithmetic library.\n"
            "Recompile with MPFI support enabled to run this test.\n");
    return 77; // Exit code 77: skip test.
#endif // WITH_MPFI
}
