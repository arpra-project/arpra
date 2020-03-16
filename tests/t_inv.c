/*
 * t_inv.c -- Test the arpra_inv function.
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

int main (int argc, char *argv[])
{
    const arpra_prec prec = 24;
    const arpra_prec prec_internal = 256;
    const arpra_uint test_n = 100000;
    unsigned i, fail_n;

    FILE *unshared_log;
    unshared_log = fopen("inv_unshared.log", "w");

    // Init test.
    test_fixture_init(prec, prec_internal);
    test_log_init("inv");
    test_rand_init();
    fail_n = 0;

    // Run test.
    for (i = 0; i < test_n; i++) {
        test_rand_arpra(&x1_A, TEST_RAND_MIXED, TEST_RAND_SMALL);

        // Pass criteria:
        // 1) Arpra x1 contains 0 and Arpra y = Inf.
        // 2) Arpra y contains MPFI y.
        test_univariate(arpra_inv, mpfi_inv);
        if (arpra_has_zero_p(&x1_A) && arpra_inf_p(&y_A)) {
            test_log_printf("Result: PASS\n\n");
        }
        else if (mpfr_greaterequal_p(&(y_I->left), &(y_A.true_range.left))
                 && mpfr_lessequal_p(&(y_I->right), &(y_A.true_range.right))) {
            test_log_printf("Result: PASS\n\n");
        }
        else {
            test_log_printf("Result: FAIL\n\n");
            fail_n++;
        }

        mpfr_out_str(unshared_log, 10, 40, y_A_diam_rel, MPFR_RNDN);
        fputs("\n", unshared_log);

    }

    fclose(unshared_log);

    // Cleanup test.
    printf("%lu out of %lu failed.\n", fail_n, test_n);
    test_fixture_clear();
    test_log_clear();
    test_rand_clear();
    return fail_n > 0;
}
