/*
 * t_div.c -- Test the arpra_div function.
 *
 * Copyright 2017-2018 James Paul Turner.
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

int main (int argc, char *argv[])
{
#ifdef WITH_MPFI
    const arpra_prec prec = 53;
    const arpra_prec prec_internal = 128;
    const arpra_uint test_n = 100000;
    arpra_uint i, fail, fail_n;

    // Init test.
    test_fixture_init(prec, prec_internal);
    test_log_init("div");
    test_rand_init();
    fail_n = 0;

    // Run test.
    for (i = 0; i < test_n; i++) {
        fail = 0;
        test_rand_arpra(&x_A, TEST_RAND_MIXED, TEST_RAND_SMALL);
        test_rand_arpra(&y_A, TEST_RAND_MIXED, TEST_RAND_SMALL);

        // Pass criteria (unshared symbols):
        // 1) Arpra x contains 0, Arpra y contains 0 and Arpra z = NaN.
        // 2) Arpra y contains 0 and Arpra z = Inf.
        // 3) Arpra z contains MPFI z.
        // 4) Arpra z unbounded and MPFI z unbounded.
        test_bivariate(arpra_div, mpfi_div);
        if (arpra_has_zero_p(&x_A) && arpra_has_zero_p(&y_A) && arpra_nan_p(&z_A)) {
            test_log_printf("Result (unshared symbols): PASS\n\n");
        }
        else if (arpra_has_zero_p(&y_A) && arpra_inf_p(&z_A)) {
            test_log_printf("Result (unshared symbols): PASS\n\n");
        }
        else if (mpfr_greaterequal_p(&(z_I->left), &(z_AI->left))
                 && mpfr_lessequal_p(&(z_I->right), &(z_AI->right))) {
            test_log_printf("Result (unshared symbols): PASS\n\n");
        }
        else if (!arpra_bounded_p(&z_A) && !mpfi_bounded_p(z_I)) {
            test_log_printf("Result (random shared symbols): PASS\n\n");
        }
        else {
            test_log_printf("Result (unshared symbols): FAIL\n\n");
            fail = 1;
        }

        // Pass criteria (random shared symbols):
        // 1) Arpra x contains 0, Arpra y contains 0 and Arpra z = NaN.
        // 2) Arpra y contains 0 and Arpra z = Inf.
        // 1) bounded(Arpra z) = bounded(MPFI z).
        test_share_rand_syms(&x_A, &y_A);
        test_bivariate(arpra_div, mpfi_div);
        if (arpra_has_zero_p(&x_A) && arpra_has_zero_p(&y_A) && arpra_nan_p(&z_A)) {
            test_log_printf("Result (random shared symbols): PASS\n\n");
        }
        else if (arpra_has_zero_p(&y_A) && arpra_inf_p(&z_A)) {
            test_log_printf("Result (random shared symbols): PASS\n\n");
        }
        else if (arpra_bounded_p(&z_A) == mpfi_bounded_p(z_I)) {
            test_log_printf("Result (random shared symbols): PASS\n\n");
        }
        else {
            test_log_printf("Result (random shared symbols): FAIL\n\n");
            fail = 1;
        }

        // Pass criteria (all shared symbols):
        // 1) Arpra x contains 0, Arpra y contains 0 and Arpra z = NaN.
        // 2) Arpra y contains 0 and Arpra z = Inf.
        // 1) bounded(Arpra z) = bounded(MPFI z).
        test_share_all_syms(&x_A, &y_A);
        test_bivariate(arpra_div, mpfi_div);
        if (arpra_has_zero_p(&x_A) && arpra_has_zero_p(&y_A) && arpra_nan_p(&z_A)) {
            test_log_printf("Result (all shared symbols): PASS\n\n");
        }
        else if (arpra_has_zero_p(&y_A) && arpra_inf_p(&z_A)) {
            test_log_printf("Result (all shared symbols): PASS\n\n");
        }
        else if (arpra_bounded_p(&z_A) == mpfi_bounded_p(z_I)) {
            test_log_printf("Result (all shared symbols): PASS\n\n");
        }
        else {
            test_log_printf("Result (all shared symbols): FAIL\n\n");
            fail = 1;
        }

        if (fail) fail_n++;
    }

    // Cleanup test.
    printf("%lu out of %lu failed.\n", fail_n, test_n);
    test_fixture_clear();
    test_log_clear();
    test_rand_clear();
    return fail_n > 0;

#else // WITH_MPFI
    fprintf(stderr,
            "This test uses the MPFI interval arithmetic library.\n"
            "Recompile with MPFI support enabled to run this test.\n");
    return 77; // Exit code 77: skip test.
#endif // WITH_MPFI
}
