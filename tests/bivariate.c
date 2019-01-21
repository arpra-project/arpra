/*
 * bivariate.c -- Compare bivariate Arpra and MPFI functions.
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

void test_bivariate (
    void (*f_arpra) (arpra_range *z, const arpra_range *x, const arpra_range *y),
    int  (*f_mpfi) (mpfi_ptr z, mpfi_srcptr x, mpfi_srcptr y))
{
    // Convert arguments.
    test_log_mpfi(&(x_A.true_range), "x  ");
    test_log_mpfi(&(y_A.true_range), "y  ");

    // Compute z with MPFI.
    f_mpfi(&z_I, &(x_A.true_range), &(y_A.true_range));
    test_log_mpfi(&z_I, "z_I");

    // Compute z with Arpra.
    f_arpra(&z_A, &x_A, &y_A);
    test_log_mpfi(&(z_A.true_range), "z_A");

    // Compute relative diameter difference.
    mpfi_diam_rel(&rdiam_I, &z_I);
    mpfi_diam_rel(&rdiam_A, &(z_A.true_range));
    mpfr_sub(&rdiam_diff, &rdiam_I, &rdiam_A, MPFR_RNDN);
    test_log_mpfr(&rdiam_diff, "z_D");
}
