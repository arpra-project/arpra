/*
 * t_div.c -- Test the mpfa_div function.
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
    const mpfa_prec_t prec = 53;
    const mpfa_prec_t prec_internal = 128;
    const mpfa_uint_t n_tests = 100000;

    mpfa_t x_A, y_A, z_A;
    mpfi_t x_I, y_I, z_I, z_AI;
    mpfr_t rdiam_I, rdiam_AI, diff;
    mpfa_uint_t i, n_fail, total_fail;

    // Init test.
    test_rand_init();
    mpfa_set_internal_prec(prec_internal);
    mpfa_init2(x_A, prec);
    mpfa_init2(y_A, prec);
    mpfa_init2(z_A, prec);
    mpfi_init2(x_I, prec);
    mpfi_init2(y_I, prec);
    mpfi_init2(z_I, prec);
    mpfi_init2(z_AI, prec);
    mpfr_init2(rdiam_I, prec_internal);
    mpfr_init2(rdiam_AI, prec_internal);
    mpfr_init2(diff, prec_internal);
    total_fail = 0;

    // Test without shared symbols.
    test_log_init("div_not_shared");
    for (n_fail = 0, i = 0; i < n_tests; i++) {
        // Set random x and y.
        test_rand_mpfa(x_A, TEST_RAND_SMALL);
        test_rand_mpfa(y_A, TEST_RAND_SMALL);
        mpfa_get_mpfi(x_I, x_A);
        mpfa_get_mpfi(y_I, y_A);

        // Compute z with MPFA and MPFI.
        mpfa_div(z_A, x_A, y_A);
        mpfi_div(z_I, x_I, y_I);
        mpfa_get_mpfi(z_AI, z_A);

        // Compute relative diameter difference in z.
        mpfi_diam_rel(rdiam_AI, z_AI);
        mpfi_diam_rel(rdiam_I, z_I);
        mpfr_sub(diff, rdiam_I, rdiam_AI, MPFR_RNDN);

        // Log operands and results.
        if (mpfi_is_inside(z_I, z_AI)) {
            test_log_printf("#%lu: Pass\n", i);
        }
        else if (!mpfi_bounded_p(z_I) && !mpfi_bounded_p(z_AI)) {
            test_log_printf("#%lu: Pass\n", i);
        }
        else {
            test_log_printf("#%lu: Fail\n", i);
            n_fail++;
        }
        test_log_mpfi(x_I, "x");
        test_log_mpfi(y_I, "y");
        test_log_mpfi(z_I, "zI");
        test_log_mpfi(z_AI, "zA");
        test_log_mpfr(diff, "zD");
        test_log_printf("\n");
    }
    total_fail += n_fail;
    printf("Not shared sym: %llu out of %llu failed.\n", n_fail, n_tests);
    test_log_clear();

    // Test with some shared symbols.
    test_log_init("div_some_shared");
    for (n_fail = 0, i = 0; i < n_tests; i++) {
        // Set random x and y.
        test_rand_mpfa(x_A, TEST_RAND_SMALL);
        test_rand_mpfa(y_A, TEST_RAND_SMALL);
        mpfa_get_mpfi(x_I, x_A);
        mpfa_get_mpfi(y_I, y_A);

        // Randomly share symbols.
        test_share_syms(x_A, y_A, 5);

        // Compute z with MPFA and MPFI.
        mpfa_div(z_A, x_A, y_A);
        mpfi_div(z_I, x_I, y_I);
        mpfa_get_mpfi(z_AI, z_A);

        // Compute relative diameter difference in z.
        mpfi_diam_rel(rdiam_AI, z_AI);
        mpfi_diam_rel(rdiam_I, z_I);
        mpfr_sub(diff, rdiam_I, rdiam_AI, MPFR_RNDN);

        // Log operands and results.
        if (mpfi_bounded_p(z_I) == mpfi_bounded_p(z_AI)) {
            test_log_printf("#%lu: Pass\n", i);
        }
        else {
            test_log_printf("#%lu: Fail\n", i);
            n_fail++;
        }
        test_log_mpfi(x_I, "x");
        test_log_mpfi(y_I, "y");
        test_log_mpfi(z_I, "zI");
        test_log_mpfi(z_AI, "zA");
        test_log_mpfr(diff, "zD");
        test_log_printf("\n");
    }
    total_fail += n_fail;
    printf("Some shared sym: %llu out of %llu failed.\n", n_fail, n_tests);
    test_log_clear();

    // Test with all shared symbols.
    test_log_init("div_all_shared");
    for (n_fail = 0, i = 0; i < n_tests; i++) {
        // Set random x and y.
        test_rand_mpfa(x_A, TEST_RAND_SMALL);
        test_rand_mpfa(y_A, TEST_RAND_SMALL);
        mpfa_get_mpfi(x_I, x_A);
        mpfa_get_mpfi(y_I, y_A);

        // Share all symbols.
        test_share_syms(x_A, y_A, 9);

        // Compute z with MPFA and MPFI.
        mpfa_div(z_A, x_A, y_A);
        mpfi_div(z_I, x_I, y_I);
        mpfa_get_mpfi(z_AI, z_A);

        // Compute relative diameter difference in z.
        mpfi_diam_rel(rdiam_AI, z_AI);
        mpfi_diam_rel(rdiam_I, z_I);
        mpfr_sub(diff, rdiam_I, rdiam_AI, MPFR_RNDN);

        // Log operands and results.
        if (mpfi_bounded_p(z_I) == mpfi_bounded_p(z_AI)) {
            test_log_printf("#%lu: Pass\n", i);
        }
        else {
            test_log_printf("#%lu: Fail\n", i);
            n_fail++;
        }
        test_log_mpfi(x_I, "x");
        test_log_mpfi(y_I, "y");
        test_log_mpfi(z_I, "zI");
        test_log_mpfi(z_AI, "zA");
        test_log_mpfr(diff, "zD");
        test_log_printf("\n");
    }
    total_fail += n_fail;
    printf("All shared sym: %llu out of %llu failed.\n", n_fail, n_tests);
    test_log_clear();

    // Cleanup test.
    mpfa_clear(x_A);
    mpfa_clear(y_A);
    mpfa_clear(z_A);
    mpfi_clear(x_I);
    mpfi_clear(y_I);
    mpfi_clear(z_I);
    mpfi_clear(z_AI);
    mpfr_clear(rdiam_I);
    mpfr_clear(rdiam_AI);
    mpfr_clear(diff);
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
