/*
 * bivariate.c -- Compare bivariate MPFA and MPFI functions.
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

void test_bivariate (
    void (*f_MPFA) (mpfa_ptr z, mpfa_srcptr x, mpfa_srcptr y),
    int  (*f_MPFI) (mpfi_ptr z, mpfi_srcptr x, mpfi_srcptr y))
{
    // Convert arguments.
    mpfa_get_mpfi(x_I, x_A);
    mpfa_get_mpfi(y_I, y_A);
    test_log_mpfi(x_I, "x  ");
    test_log_mpfi(y_I, "y  ");

    // Compute z with MPFI.
    f_MPFI(z_I, x_I, y_I);
    test_log_mpfi(z_I, "z_I");

    // Compute z with MPFA (unshared symbols).
    f_MPFA(z_A, x_A, y_A);
    mpfa_get_mpfi(z_AI, z_A);
    test_log_mpfi(z_AI, "z_A");

    // Compute relative diameter difference.
    mpfi_diam_rel(rdiam_I, z_I);
    mpfi_diam_rel(rdiam_AI, z_AI);
    mpfr_sub(rdiam_diff, rdiam_AI, rdiam_I, MPFR_RNDN);
    test_log_mpfr(rdiam_diff, "z_D");
}
