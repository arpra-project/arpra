/*
 * helper_mpfr_rnderr.c -- Compute MPFR functions with rounding error.
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

#include "arpra-impl.h"

/*
 * This function assumes that the MPFR function which computed y returned
 * a nonzero ternary value, and thus y is inexact.
 *
 * Note: exponents of MPFR numbers are one greater than eponents of equivalent
 * IEEE-754 floating-point numbers, since MPFR significands are in [0.5, 1.0).
 */

void arpra_helper_mpfr_rnderr (mpfr_ptr err, mpfr_rnd_t rnd, mpfr_srcptr y)
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
        if ((rnd == MPFR_RNDN) || (rnd == MPFR_RNDNA)) {
            // Rounding error is 0.5 ULP(y) = 2^(e-p-1).
            mpfr_set_si_2exp(temp, 1, (e - p - 1), MPFR_RNDU);
        }
        else {
            // Rounding error is ULP(y) = 2^(e-p).
            mpfr_set_si_2exp(temp, 1, (e - p), MPFR_RNDU);
        }
    }

    // Add rounding error to total.
    mpfr_add(err, err, temp, MPFR_RNDU);

    // Clear vars.
    mpfr_clear(temp);
}




void arpra_helper_term_mul (mpfr_ptr error, mpfr_ptr y, mpfr_srcptr x1, mpfi_srcptr alpha)
{
    mpfr_t temp1, temp2;
    mpfi_t y_range;
    arpra_prec prec_internal;

    // Initialise vars.
    prec_internal = arpra_get_internal_precision();
    mpfr_init2(temp1, prec_internal);
    mpfr_init2(temp2, prec_internal);
    mpfi_init2(y_range, prec_internal);

    // y = (alpha * x1)
    mpfi_mul_fr(y_range, alpha, x1);

    // y mid
    mpfr_div_ui(temp1, &(y_range.left), 2, MPFR_RNDD);
    mpfr_div_ui(temp2, &(y_range.right), 2, MPFR_RNDU);
    mpfr_add(y, temp1, temp2, MPFR_RNDN);

    // y rad
    mpfr_sub(temp1, y, &(y_range.left), MPFR_RNDU);
    mpfr_sub(temp2, &(y_range.right), y, MPFR_RNDU);
    mpfr_max(temp1, temp1, temp2, MPFR_RNDU);
    mpfr_add(error, error, temp1, MPFR_RNDU);

    // Clear temp vars.
    mpfr_clear(temp1);
    mpfr_clear(temp2);
    mpfi_clear(y_range);
}


void arpra_helper_term_fma (mpfr_ptr error, mpfr_ptr y, mpfr_srcptr x1,
                            mpfi_srcptr alpha, mpfi_srcptr gamma)
{
    mpfr_t temp1, temp2;
    mpfi_t y_range, alpha_x1;
    arpra_prec prec_internal;

    // a * b needs precision prec(a) + prec(b) to be exact.

    // Initialise vars.
    prec_internal = arpra_get_internal_precision();
    mpfr_init2(temp1, prec_internal);
    mpfr_init2(temp2, prec_internal);
    mpfi_init2(y_range, prec_internal);
    mpfi_init2(alpha_x1, (mpfi_get_prec(alpha) + mpfr_get_prec(x1)));

    // y = (alpha * x1) + (gamma)
    mpfi_mul_fr(alpha_x1, alpha, x1);
    mpfi_add(y_range, alpha_x1, gamma);

    // y mid
    mpfr_div_ui(temp1, &(y_range.left), 2, MPFR_RNDD);
    mpfr_div_ui(temp2, &(y_range.right), 2, MPFR_RNDU);
    mpfr_add(y, temp1, temp2, MPFR_RNDN);

    // y rad
    mpfr_sub(temp1, y, &(y_range.left), MPFR_RNDU);
    mpfr_sub(temp2, &(y_range.right), y, MPFR_RNDU);
    mpfr_max(temp1, temp1, temp2, MPFR_RNDU);
    mpfr_add(error, error, temp1, MPFR_RNDU);

    // Clear temp vars.
    mpfr_clear(temp1);
    mpfr_clear(temp2);
    mpfi_clear(y_range);
    mpfi_clear(alpha_x1);
}


void arpra_helper_term_fmma (mpfr_ptr error, mpfr_ptr y, mpfr_srcptr x1, mpfr_srcptr x2,
                             mpfi_srcptr alpha, mpfi_srcptr beta)
{
    mpfr_t temp1, temp2;
    mpfi_t y_range, alpha_x1, beta_x2;
    arpra_prec prec_internal;

    // a * b needs precision prec(a) + prec(b) to be exact.

    // Initialise vars.
    prec_internal = arpra_get_internal_precision();
    mpfr_init2(temp1, prec_internal);
    mpfr_init2(temp2, prec_internal);
    mpfi_init2(y_range, prec_internal);
    mpfi_init2(alpha_x1, (mpfi_get_prec(alpha) + mpfr_get_prec(x1)));
    mpfi_init2(beta_x2, (mpfi_get_prec(beta) + mpfr_get_prec(x2)));

    // y = (alpha * x1) + (beta * x2)
    mpfi_mul_fr(alpha_x1, alpha, x1);
    mpfi_mul_fr(beta_x2, beta, x2);
    mpfi_add(y_range, alpha_x1, beta_x2);

    // y mid
    mpfr_div_ui(temp1, &(y_range.left), 2, MPFR_RNDD);
    mpfr_div_ui(temp2, &(y_range.right), 2, MPFR_RNDU);
    mpfr_add(y, temp1, temp2, MPFR_RNDN);

    // y rad
    mpfr_sub(temp1, y, &(y_range.left), MPFR_RNDU);
    mpfr_sub(temp2, &(y_range.right), y, MPFR_RNDU);
    mpfr_max(temp1, temp1, temp2, MPFR_RNDU);
    mpfr_add(error, error, temp1, MPFR_RNDU);

    // Clear temp vars.
    mpfr_clear(temp1);
    mpfr_clear(temp2);
    mpfi_clear(y_range);
    mpfi_clear(alpha_x1);
    mpfi_clear(beta_x2);
}


void arpra_helper_term_fmmaa (mpfr_ptr error, mpfr_ptr y, mpfr_srcptr x1, mpfr_srcptr x2,
                              mpfi_srcptr alpha, mpfi_srcptr beta, mpfi_srcptr gamma)
{
    mpfr_t temp1, temp2;
    mpfi_t y_range, alpha_x1, beta_x2;
    arpra_prec prec_internal;

    // a * b needs precision prec(a) + prec(b) to be exact.

    // Initialise vars.
    prec_internal = arpra_get_internal_precision();
    mpfr_init2(temp1, prec_internal);
    mpfr_init2(temp2, prec_internal);
    mpfi_init2(y_range, prec_internal);
    mpfi_init2(alpha_x1, (mpfi_get_prec(alpha) + mpfr_get_prec(x1)));
    mpfi_init2(beta_x2, (mpfi_get_prec(beta) + mpfr_get_prec(x2)));

    // y = (alpha * x1) + (beta * x2) + (gamma)
    mpfi_mul_fr(alpha_x1, alpha, x1);
    mpfi_mul_fr(beta_x2, beta, x2);
    mpfr_sum(&(y_range.left), (mpfr_ptr[3]) {&(alpha_x1.left), &(beta_x2.left), &(gamma.left)}, 3, MPFR_RNDD);
    mpfr_sum(&(y_range.right), (mpfr_ptr[3]) {&(alpha_x1.right), &(beta_x2.right), &(gamma.right)}, 3, MPFR_RNDU);

    // y mid
    mpfr_div_ui(temp1, &(y_range.left), 2, MPFR_RNDD);
    mpfr_div_ui(temp2, &(y_range.right), 2, MPFR_RNDU);
    mpfr_add(y, temp1, temp2, MPFR_RNDN);

    // y rad
    mpfr_sub(temp1, y, &(y_range.left), MPFR_RNDU);
    mpfr_sub(temp2, &(y_range.right), y, MPFR_RNDU);
    mpfr_max(temp1, temp1, temp2, MPFR_RNDU);
    mpfr_add(error, error, temp1, MPFR_RNDU);

    // Clear temp vars.
    mpfr_clear(temp1);
    mpfr_clear(temp2);
    mpfi_clear(y_range);
    mpfi_clear(alpha_x1);
    mpfi_clear(beta_x2);
}
