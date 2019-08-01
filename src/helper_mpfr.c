/*
 * helper_mpfr.c -- Compute MPFR functions with rounding error.
 *
 * Copyright 2017-2019 James Paul Turner.
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

/*
 * This function assumes that the MPFR function which computed y returned
 * a nonzero ternary value, and thus y is inexact.
 *
 * Note: exponents of MPFR numbers are one greater than eponents of equivalent
 * IEEE-754 floating-point numbers, since MPFR significands are in [0.5, 1.0).
 */

static void arpra_helper_rnd_err (mpfr_srcptr y, mpfr_rnd_t rnd_mode, mpfr_ptr rnd_err)
{
    mpfr_t temp;
    mpfr_exp_t e;
    mpfr_prec_t p;
    arpra_prec internal_prec;

    // Initialise vars.
    internal_prec = arpra_get_internal_precision();
    mpfr_init2(temp, internal_prec);

    // Was y flushed to zero?
    if (mpfr_zero_p(y)) {
        // Rounding error is nextabove(0).
        mpfr_set_zero(temp, 1);
        mpfr_nextabove(temp);
    }
    else {
        e = mpfr_get_exp(y);
        p = mpfr_get_prec(y);

        // Nearest or directed rounding?
        if ((rnd_mode == MPFR_RNDN) || (rnd_mode == MPFR_RNDNA)) {
            // Rounding error is 0.5 ULP(y) = 2^(e-p-1).
            mpfr_set_si_2exp(temp, 1, (e - p - 1), MPFR_RNDU);
        }
        else {
            // Rounding error is ULP(y) = 2^(e-p).
            mpfr_set_si_2exp(temp, 1, (e - p), MPFR_RNDU);
        }
    }

    // Add rounding error to total.
    mpfr_add(rnd_err, rnd_err, temp, MPFR_RNDU);

    // Clear vars.
    mpfr_clear(temp);
}

void arpra_helper_mpfr_f1 (int (*f) (mpfr_ptr, mpfr_srcptr, mpfr_rnd_t),
                           mpfr_ptr y, mpfr_srcptr x1, mpfr_rnd_t rnd_mode,
                           mpfr_ptr rnd_err)
{
    if (f(y, x1, rnd_mode)) {
        arpra_helper_rnd_err(y, rnd_mode, rnd_err);
    }
}

void arpra_helper_mpfr_f2 (int (*f) (mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mpfr_rnd_t),
                           mpfr_ptr y, mpfr_srcptr x1, mpfr_srcptr x2, mpfr_rnd_t rnd_mode,
                           mpfr_ptr rnd_err)
{
    if (f(y, x1, x2, rnd_mode)) {
        arpra_helper_rnd_err(y, rnd_mode, rnd_err);
    }
}

void arpra_helper_mpfr_f3 (int (*f) (mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mpfr_srcptr, mpfr_rnd_t),
                           mpfr_ptr y, mpfr_srcptr x1, mpfr_srcptr x2, mpfr_srcptr x3, mpfr_rnd_t rnd_mode,
                           mpfr_ptr rnd_err)
{
    if (f(y, x1, x2, x3, rnd_mode)) {
        arpra_helper_rnd_err(y, rnd_mode, rnd_err);
    }
}

void arpra_helper_mpfr_fmma (mpfr_ptr y, mpfr_srcptr x1, mpfr_srcptr x2,
                             mpfr_srcptr a, mpfr_srcptr b,
                             mpfr_ptr rnd_err)
{
    mpfr_t ax1, bx2;

    // a * x1 needs precision prec(a) + prec(x1) to be exact.
    mpfr_init2(ax1, (mpfr_get_prec(a) + mpfr_get_prec(x1)));
    mpfr_mul(ax1, a, x1, MPFR_RNDN);

    // b * x2 needs precision prec(b) + prec(x2) to be exact.
    mpfr_init2(bx2, (mpfr_get_prec(b) + mpfr_get_prec(x2)));
    mpfr_mul(bx2, b, x2, MPFR_RNDN);

    // y = (a * x1) + (b * x2)
    if (mpfr_add(y, ax1, bx2, MPFR_RNDN)) {
        arpra_helper_rnd_err(y, MPFR_RNDN, rnd_err);
    }

    // Clear temp vars.
    mpfr_clear(ax1);
    mpfr_clear(bx2);
}

void arpra_helper_mpfr_fmmaa (mpfr_ptr y, mpfr_srcptr x1, mpfr_srcptr x2,
                              mpfr_srcptr a, mpfr_srcptr b, mpfr_srcptr c,
                              mpfr_ptr rnd_err)
{
    mpfr_t ax1, bx2;

    // a * x1 needs precision prec(a) + prec(x1) to be exact.
    mpfr_init2(ax1, (mpfr_get_prec(a) + mpfr_get_prec(x1)));
    mpfr_mul(ax1, a, x1, MPFR_RNDN);

    // b * x2 needs precision prec(b) + prec(x2) to be exact.
    mpfr_init2(bx2, (mpfr_get_prec(b) + mpfr_get_prec(x2)));
    mpfr_mul(bx2, b, x2, MPFR_RNDN);

    // Newer MPFR 4 syntax
    //if (mpfr_sum(y, (mpfr_ptr[3]) {ax1, bx2, c}, 3, MPFR_RNDN)) {

    // y = (a * x1) + (b * x2) + (c)
    if (mpfr_sum(y, (mpfr_ptr[3]) {ax1, bx2, (mpfr_ptr) c}, 3, MPFR_RNDN)) {
        arpra_helper_rnd_err(y, MPFR_RNDN, rnd_err);
    }

    // Clear temp vars.
    mpfr_clear(ax1);
    mpfr_clear(bx2);
}
