/*
 * ext_mpfr.c -- Arpra extensions to the MPFR library.
 *
 * Copyright 2018-2020 James Paul Turner.
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

#include "arpra-impl.h"

// TODO: USE MPFR 4 FMMA
int arpra_ext_mpfr_fmma (mpfr_ptr y, mpfr_srcptr x1, mpfr_srcptr x2,
                         mpfr_srcptr x3, mpfr_srcptr x4, mpfr_rnd_t rnd)
{
    mpfr_t x1x2, x3x4;
    int ternary;

    // x1 * x2 needs precision prec(x1) + prec(x2) to be exact.
    mpfr_init2(x1x2, (mpfr_get_prec(x1) + mpfr_get_prec(x2)));
    mpfr_mul(x1x2, x1, x2, MPFR_RNDN);

    // x3 * x4 needs precision prec(x3) + prec(x4) to be exact.
    mpfr_init2(x3x4, (mpfr_get_prec(x3) + mpfr_get_prec(x4)));
    mpfr_mul(x3x4, x3, x4, MPFR_RNDN);

    // y = (x1 * x2) + (x3 * x4)
    ternary = mpfr_add(y, x1x2, x3x4, rnd);

    // Clear temp vars.
    mpfr_clear(x1x2);
    mpfr_clear(x3x4);

    return ternary;
}

// TOOD: use MPFR 4 SUM SYNTAX
int arpra_ext_mpfr_fmmaa (mpfr_ptr y, mpfr_srcptr x1, mpfr_srcptr x2,
                          mpfr_srcptr x3, mpfr_srcptr x4, mpfr_srcptr x5, mpfr_rnd_t rnd)
{
    mpfr_t x1x2, x3x4;
    int ternary;

    // x1 * x2 needs precision prec(x1) + prec(x2) to be exact.
    mpfr_init2(x1x2, (mpfr_get_prec(x1) + mpfr_get_prec(x2)));
    mpfr_mul(x1x2, x1, x2, MPFR_RNDN);

    // x3 * x4 needs precision prec(x3) + prec(x4) to be exact.
    mpfr_init2(x3x4, (mpfr_get_prec(x3) + mpfr_get_prec(x4)));
    mpfr_mul(x3x4, x3, x4, MPFR_RNDN);

    // y = (x1 * x2) + (x3 * x4) + (x5)
    //if (mpfr_sum(y, (mpfr_ptr[3]) {x1x2, x3x4, x5}, 3, rnd)) {
    ternary = mpfr_sum(y, (mpfr_ptr[3]) {x1x2, x3x4, (mpfr_ptr) x5}, 3, rnd);

    // Clear temp vars.
    mpfr_clear(x1x2);
    mpfr_clear(x3x4);

    return ternary;
}

int arpra_ext_mpfr_sum (arpra_mpfr *y, arpra_mpfr *x,
                        const arpra_uint n, const mpfr_rnd_t rnd)
{
    arpra_mpfr **buffer_mpfr_ptr;
    arpra_uint i;

    // Save number pointers to buffer.
    buffer_mpfr_ptr = arpra_helper_buffer_mpfr_ptr(n);
    for (i = 0; i < n; i++) {
        buffer_mpfr_ptr[i] = &(x[i]);
    }

    // Sum the numbers.
    return mpfr_sum(y, buffer_mpfr_ptr, n, rnd);
}

int arpra_ext_mpfr_sumabs (arpra_mpfr *y, arpra_mpfr *x,
                           const arpra_uint n, const mpfr_rnd_t rnd)
{
    arpra_mpfr *buffer_mpfr;
    arpra_uint i;

    // Save absolute value numbers to buffer.
    buffer_mpfr = arpra_helper_buffer_mpfr(n);
    for (i = 0; i < n; i++) {
        buffer_mpfr[i] = x[i];
        buffer_mpfr[i]._mpfr_sign = 1;
    }

    // Sum the absolute value numbers.
    return arpra_ext_mpfr_sum(y, buffer_mpfr, n, rnd);
}
