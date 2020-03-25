/*
 * univariate.c -- Compare univariate Arpra and MPFI functions.
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

void test_univariate (
    void (*f_arpra) (arpra_range *y, const arpra_range *x1),
    int  (*f_mpfi) (mpfi_ptr y, mpfi_srcptr x1))
{
    // Convert arguments.
    test_log_mpfi(&(x1_A.true_range), "x1  ");

    // Compute y with MPFI.
    f_mpfi(y_I, &(x1_A.true_range));
    test_log_mpfi(y_I, "y_I");

    // Compute y with Arpra.
    f_arpra(&y_A, &x1_A);
    test_log_mpfi(&(y_A.true_range), "y_A");

    // Compute Arpra diameter relative to MPFI diameter.
    mpfr_sub(y_I_diam, &(y_I->right), &(y_I->left), MPFR_RNDU);
    test_log_mpfr(y_I_diam, "y_I_diam    ");
    mpfr_sub(y_A_diam, &(y_A.true_range.right), &(y_A.true_range.left), MPFR_RNDU);
    test_log_mpfr(y_A_diam, "y_A_diam    ");
    mpfr_div(y_A_diam_rel, y_A_diam, y_I_diam, MPFR_RNDN);
    test_log_mpfr(y_A_diam_rel, "y_A_diam_rel");
}
