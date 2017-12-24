/*
 * t_log.c -- Test the arpra_log function.
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

int main (int argc, char *argv[])
{
#ifdef WITH_MPFI
    const arpra_prec_t prec = 53;
    const arpra_prec_t prec_internal = 128;
    const arpra_uint_t test_n = 100000;
    arpra_uint_t i, fail_n;

    // Init test.
    test_fixture_init(prec, prec_internal);
    test_log_init("log");
    test_rand_init();
    fail_n = 0;

    // Run test.
    for (i = 0; i < test_n; i++) {
        test_rand_arpra(x_A, TEST_RAND_POS, TEST_RAND_SMALL);

        // Pass criteria:
        // 1) ArPRA x has negative and ArPRA z = NaN.
        // 2) ArPRA x contains 0 and ArPRA z = Inf.
        // 3) ArPRA z contains MPFI z.
        test_univariate(arpra_log, mpfi_log);
        if (arpra_has_neg_p(x_A) && arpra_nan_p(z_A)) {
            test_log_printf("Result: PASS\n\n");
        }
        else if (arpra_has_zero_p(x_A) && arpra_inf_p(z_A)) {
            test_log_printf("Result: PASS\n\n");
        }
        else if (mpfr_greaterequal_p(&(z_I->left), &(z_AI->left))
                 && mpfr_lessequal_p(&(z_I->right), &(z_AI->right))) {
            test_log_printf("Result: PASS\n\n");
        }
        else {
            test_log_printf("Result: FAIL\n\n");
            fail_n++;
        }
    }

    // Cleanup test.
    printf("%llu out of %llu failed.\n", fail_n, test_n);
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
