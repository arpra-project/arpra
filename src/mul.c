/*
 * mul.c -- Multiply one Arpra range by another.
 *
 * Copyright 2016-2020 James Paul Turner.
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

static arpra_mul_method mul_method = ARPRA_DEFAULT_MUL_METHOD;

arpra_mul_method arpra_get_mul_method ()
{
    return mul_method;
}

void arpra_set_mul_method (arpra_mul_method new_mul_method)
{
    mul_method = new_mul_method;
}

static void mul_err_trivial (mpfr_ptr error, const arpra_range *x1, const arpra_range *x2)
{
    mpfr_t temp;
    arpra_prec prec_internal;

    // Init temp vars.
    prec_internal = arpra_get_internal_precision();
    mpfr_init2(temp, prec_internal);

    // Trivial approximation error is rad(x1) * rad(x2).
    mpfr_mul(temp, &(x1->radius), &(x2->radius), MPFR_RNDU);
    mpfr_add(error, error, temp, MPFR_RNDU);

    // Clear temp vars.
    mpfr_clear(temp);
}

static void mul_err_rump_kashiwagi (mpfr_ptr error, const arpra_range *x1, const arpra_range *x2)
{
    mpfr_t temp;
    mpfr_t x1ix2j, x1jx2i;
    mpfr_t x1ix2i_pos_error, x1ix2i_neg_error;
    arpra_uint x1i_idx, x1j_idx, x2i_idx, x2j_idx;
    arpra_int x1HasNext, x2HasNext;
    arpra_prec prec_internal;

    /*
     * The approximation of the quadratic term of arpra_mul is defined the same as in (26) of
     * S. M. Rump and M. Kashiwagi, Implementation and improvements of affine arithmetic,
     * Nonlinear Theory an Its Applications, IEICE, vol. 6, no. 3, pp. 341-359, 2015.
     */

    // Init temp vars.
    prec_internal = arpra_get_internal_precision();
    mpfr_init2(temp, prec_internal);
    mpfr_init2(x1ix2j, prec_internal);
    mpfr_init2(x1jx2i, prec_internal);
    mpfr_init2(x1ix2i_pos_error, prec_internal);
    mpfr_init2(x1ix2i_neg_error, prec_internal);
    mpfr_set_zero(x1ix2i_pos_error, 1);
    mpfr_set_zero(x1ix2i_neg_error, 1);

    x1i_idx = 0;
    x2i_idx = 0;
    while ((x1i_idx < x1->nTerms) && (x2i_idx < x2->nTerms)) {
        if (x1->symbols[x1i_idx] < x2->symbols[x2i_idx]) {
            for (x2j_idx = x2i_idx; x2j_idx < x2->nTerms; x2j_idx++) {
                // x1 has symbol i, and x2 has symbol j, so error += abs(x1[i] * x2[j])
                mpfr_mul(x1ix2j, &(x1->deviations[x1i_idx]), &(x2->deviations[x2j_idx]), MPFR_RNDA);
                mpfr_abs(x1ix2j, x1ix2j, MPFR_RNDU);
                mpfr_add(error, error, x1ix2j, MPFR_RNDU);
            }
            x1i_idx++;
        }
        else if (x2->symbols[x2i_idx] < x1->symbols[x1i_idx]) {
            for (x1j_idx = x1i_idx; x1j_idx < x1->nTerms; x1j_idx++) {
                // x2 has symbol i, and x1 has symbol j, so error += abs(x1[j] * x2[i])
                mpfr_mul(x1jx2i, &(x1->deviations[x1j_idx]), &(x2->deviations[x2i_idx]), MPFR_RNDA);
                mpfr_abs(x1jx2i, x1jx2i, MPFR_RNDU);
                mpfr_add(error, error, x1jx2i, MPFR_RNDU);
            }
            x2i_idx++;
        }
        else {
            // both x1 and x2 have symbol i, so error += abs(x1[i] * x2[i])
            mpfr_mul(temp, &(x1->deviations[x1i_idx]), &(x2->deviations[x2i_idx]), MPFR_RNDA);
            if (mpfr_sgn(temp) > 0) {
                mpfr_add(x1ix2i_pos_error, x1ix2i_pos_error, temp, MPFR_RNDU);
            }
            else if (mpfr_sgn(temp) < 0) {
                mpfr_sub(x1ix2i_neg_error, x1ix2i_neg_error, temp, MPFR_RNDU);
            }

            x1j_idx = x1i_idx + 1;
            x2j_idx = x2i_idx + 1;
            x1HasNext = x1j_idx < x1->nTerms;
            x2HasNext = x2j_idx < x2->nTerms;
            while (x1HasNext || x2HasNext) {
                if ((!x2HasNext) || (x1HasNext && (x1->symbols[x1i_idx] < x2->symbols[x2i_idx]))) {
                    // both x1 and x2 have symbol i, but only x1 has symbol j, so error += abs(x1[j] * x2[i])
                    mpfr_mul(x1jx2i, &(x1->deviations[x1j_idx]), &(x2->deviations[x2i_idx]), MPFR_RNDA);
                    mpfr_abs(x1jx2i, x1jx2i, MPFR_RNDU);
                    mpfr_add(error, error, x1jx2i, MPFR_RNDU);
                    x1HasNext = ++x1j_idx < x1->nTerms;
                }
                else if ((!x1HasNext) || (x2HasNext && (x2->symbols[x2i_idx] < x1->symbols[x1i_idx]))) {
                    // both x1 and x2 have symbol i, but only x2 has symbol j, so error += abs(x1[i] * x2[j])
                    mpfr_mul(x1ix2j, &(x1->deviations[x1i_idx]), &(x2->deviations[x2j_idx]), MPFR_RNDA);
                    mpfr_abs(x1ix2j, x1ix2j, MPFR_RNDU);
                    mpfr_add(error, error, x1ix2j, MPFR_RNDU);
                    x2HasNext = ++x2j_idx < x2->nTerms;
                }
                else {
                    // both x1 and x2 have symbols i and j, so error += abs(x1[i] * x2[j] + x1[j] * x2[i])
                    mpfr_mul(x1ix2j, &(x1->deviations[x1i_idx]), &(x2->deviations[x2j_idx]), MPFR_RNDU);
                    mpfr_mul(x1jx2i, &(x1->deviations[x1j_idx]), &(x2->deviations[x2i_idx]), MPFR_RNDU);
                    mpfr_add(temp, x1ix2j, x1jx2i, MPFR_RNDU);
                    if (mpfr_sgn(temp) < 0) {
                        mpfr_mul(x1ix2j, &(x1->deviations[x1i_idx]), &(x2->deviations[x2j_idx]), MPFR_RNDD);
                        mpfr_mul(x1jx2i, &(x1->deviations[x1j_idx]), &(x2->deviations[x2i_idx]), MPFR_RNDD);
                        mpfr_add(temp, x1ix2j, x1jx2i, MPFR_RNDD);
                    }
                    mpfr_abs(temp, temp, MPFR_RNDU);
                    mpfr_add(error, error, temp, MPFR_RNDU);
                    x1HasNext = ++x1j_idx < x1->nTerms;
                    x2HasNext = ++x2j_idx < x2->nTerms;
                }
            }
            x1i_idx++;
            x2i_idx++;
        }
    }

    mpfr_max(temp, x1ix2i_pos_error, x1ix2i_neg_error, MPFR_RNDU);
    mpfr_add(error, error, temp, MPFR_RNDU);

    // Clear temp vars.
    mpfr_clear(temp);
    mpfr_clear(x1ix2j);
    mpfr_clear(x1jx2i);
    mpfr_clear(x1ix2i_pos_error);
    mpfr_clear(x1ix2i_neg_error);
}

void arpra_mul (arpra_range *y, const arpra_range *x1, const arpra_range *x2)
{
    mpfi_t ia_range;
    mpfr_t error;
    arpra_range yy;
    arpra_uint i_y, i_x1, i_x2;
    arpra_int x1HasNext, x2HasNext;
    arpra_prec prec_internal;

    // Domain violations:
    // (NaN) * (NaN) = (NaN)
    // (NaN) * (R)   = (NaN)
    // (R)   * (NaN) = (NaN)
    // (Inf) * (Inf) = (NaN)
    // (Inf) * (0)   = (NaN)
    // (0)   * (Inf) = (NaN)
    // (Inf) * (R)   = (Inf)
    // (R)   * (Inf) = (Inf)

    // Handle domain violations.
    if (arpra_nan_p(x1) || arpra_nan_p(x2)) {
        arpra_set_nan(y);
        return;
    }
    if (arpra_inf_p(x1)) {
        if (arpra_has_zero_p(x2)) {
            arpra_set_nan(y);
        }
        else {
            arpra_set_inf(y);
        }
        return;
    }
    if (arpra_inf_p(x2)) {
        if (arpra_has_zero_p(x1)) {
            arpra_set_nan(y);
        }
        else {
            arpra_set_inf(y);
        }
        return;
    }

    // Initialise vars.
    prec_internal = arpra_get_internal_precision();
    mpfi_init2(ia_range, y->precision);
    mpfr_init2(error, prec_internal);
    arpra_init2(&yy, y->precision);
    mpfr_set_zero(error, 1);

    // y[0] = x1[0] * x2[0]
    ARPRA_MPFR_RNDERR_MUL(error, MPFR_RNDN, &(yy.centre), &(x1->centre), &(x2->centre));

    // Allocate memory for deviation terms.
    yy.symbols = malloc((x1->nTerms + x2->nTerms + 1) * sizeof(arpra_uint));
    yy.deviations = malloc((x1->nTerms + x2->nTerms + 1) * sizeof(mpfr_t));

    i_y = 0;
    i_x1 = 0;
    i_x2 = 0;
    x1HasNext = x1->nTerms > 0;
    x2HasNext = x2->nTerms > 0;
    while (x1HasNext || x2HasNext) {
        mpfr_init2(&(yy.deviations[i_y]), prec_internal);

        if ((!x2HasNext) || (x1HasNext && (x1->symbols[i_x1] < x2->symbols[i_x2]))) {
            // y[i] = (x2[0] * x1[i])
            yy.symbols[i_y] = x1->symbols[i_x1];
            ARPRA_MPFR_RNDERR_MUL(error, MPFR_RNDN, &(yy.deviations[i_y]), &(x2->centre), &(x1->deviations[i_x1]));
            x1HasNext = ++i_x1 < x1->nTerms;
        }
        else if ((!x1HasNext) || (x2HasNext && (x2->symbols[i_x2] < x1->symbols[i_x1]))) {
            // y[i] = (x1[0] * x2[i])
            yy.symbols[i_y] = x2->symbols[i_x2];
            ARPRA_MPFR_RNDERR_MUL(error, MPFR_RNDN, &(yy.deviations[i_y]), &(x1->centre), &(x2->deviations[i_x2]));
            x2HasNext = ++i_x2 < x2->nTerms;
        }
        else {
            // y[i] = (x2[0] * x1[i]) + (x1[0] * x2[i])
            yy.symbols[i_y] = x1->symbols[i_x1];
            ARPRA_MPFR_RNDERR_FMMA(error, MPFR_RNDN, &(yy.deviations[i_y]), &(x2->centre), &(x1->deviations[i_x1]), &(x1->centre), &(x2->deviations[i_x2]));
            x1HasNext = ++i_x1 < x1->nTerms;
            x2HasNext = ++i_x2 < x2->nTerms;
        }
        i_y++;
    }

    // Approximation error.
    switch (mul_method) {
    case ARPRA_MUL_TRIVIAL:
        mul_err_trivial(error, x1, x2);
        break;
    case ARPRA_MUL_RUMP_KASHIWAGI:
        mul_err_rump_kashiwagi(error, x1, x2);
        break;
    }

    // Store new deviation term.
    yy.symbols[i_y] = arpra_helper_next_symbol();
    yy.deviations[i_y] = *error;
    yy.nTerms = i_y + 1;

    // MPFI multiplication
    mpfi_mul(ia_range, &(x1->true_range), &(x2->true_range));

    // Compute true_range.
    arpra_helper_compute_range(&yy);

    // Mix with IA range, and trim error term.
    arpra_helper_mix_trim(&yy, ia_range);

    // Check for NaN and Inf.
    arpra_helper_check_result(&yy);

    // Clear vars, and set y.
    mpfi_clear(ia_range);
    arpra_clear(y);
    *y = yy;
}
