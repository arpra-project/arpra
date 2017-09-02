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
    const int n_tests = 1000;
    const mpfa_prec_t prec = 53;

    int i, n_failed = 0;
    mpfa_t a_a, b_a, c_a;
    mpfi_t a_i, b_i, c_i;

    // Init test, and set working precision.
    test_rand_init();
    mpfa_init2(a_a, prec);
    mpfa_init2(b_a, prec);
    mpfa_init2(c_a, prec);
    mpfi_init2(a_i, prec);
    mpfi_init2(b_i, prec);
    mpfi_init2(c_i, prec);

    for (i = 0; i < n_tests; i++) {
        // Set random A.
        test_rand_mpfa(a_a, TEST_RAND_MIXED);
        mpfa_get_mpfi(a_i, a_a);

        // Set random B.
        test_rand_mpfa(b_a, TEST_RAND_MIXED);
        mpfa_get_mpfi(b_i, b_a);

        // Compare MPFA result with MPFI result.
        mpfa_add(c_a, a_a, b_a);
        mpfi_add(c_i, a_i, b_i);
        n_failed += test_cmp_mpfi(c_a, c_i);
    }

    printf("Failed %d out of %d.\n", n_failed, n_tests);

    // Cleanup test.
    test_rand_clear();
    mpfa_clear(a_a);
    mpfa_clear(b_a);
    mpfa_clear(c_a);
    mpfi_clear(a_i);
    mpfi_clear(b_i);
    mpfi_clear(c_i);
    mpfr_free_cache();
    return n_failed > 0;

#else // WITH_MPFI
    fprintf(stderr,
            "This test program uses the MPFI interval arithmetic library.\n"
            "Recompile MPFA with MPFI support enabled to run this test.\n");
    return 77; // Skip test code.
#endif // WITH_MPFI
}
