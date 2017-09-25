/*
 * univariate_mpfi.c -- Test bivariate MPFA function against MPFI function.
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

int test_univariate_mpfi (
    void (*f_MPFA) (mpfa_ptr z, mpfa_srcptr x),
    int  (*f_MPFI) (mpfi_ptr z, mpfi_srcptr x))
{
    int fail = 0;

    // Set and log random x.
    test_rand_mpfa(x_A, TEST_RAND_SMALL);
    test_log_mpfi(x_I, "x_I");

    // Compute and log z with MPFI.
    mpfa_get_mpfi(x_I, x_A);
    f_MPFI(z_I, x_I);
    test_log_mpfi(z_I, "z_I");

    // Compute z with MPFA.
    f_MPFA(z_A, x_A);
    mpfa_get_mpfi(z_AI, z_A);
    mpfi_diam_rel(rdiam_I, z_I);
    mpfi_diam_rel(rdiam_AI, z_AI);
    mpfr_sub(rdiam_diff, rdiam_I, rdiam_AI, MPFR_RNDN);

    // Log result.
    if ((mpfr_cmp(&(z_I->left), &(z_AI->left)) >= 0)
            || (mpfr_cmp(&(z_I->right), &(z_AI->right)) <= 0)) {
        // Fail if MPFA range does not include MPFI range.
        test_log_printf("Unshared symbol: PASS\n");
    }
    else if (!mpfi_bounded_p(z_I) && !mpfi_bounded_p(z_AI)) {
        test_log_printf("Unshared symbol: PASS\n");
    }
    else {
        test_log_printf("Unshared symbol: FAIL\n");;
        fail++;
    }
    test_log_mpfi(z_AI, "z_A");
    test_log_mpfr(rdiam_diff, "z_D");

    return fail;
}
