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

void arpra_mul (arpra_range *y, const arpra_range *x1, const arpra_range *x2)
{
    mpfi_t ia_range;
    mpfr_t temp1, temp2, error;
    arpra_range yy;
    arpra_prec prec_internal;
    arpra_uint iy, ix1, ix2;
    arpra_int x1HasNext, x2HasNext;

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
    mpfr_init2(temp1, prec_internal + 8);
    mpfr_init2(temp2, prec_internal + 8);
    mpfr_init2(error, prec_internal);
    arpra_init2(&yy, y->precision);
    mpfr_set_zero(error, 1);
    mpfr_set_zero(&(yy.radius), 1);

    // MPFI multiplication
    mpfi_mul(ia_range, &(x1->true_range), &(x2->true_range));

    // y[0] = x1[0] * x2[0]
    arpra_helper_mpfr_rnd_err_f2(error, &mpfr_mul, &(yy.centre), &(x1->centre), &(x2->centre), MPFR_RNDN);

    // Allocate memory for all possible deviation terms.
    yy.symbols = malloc((x1->nTerms + x2->nTerms + 1) * sizeof(arpra_uint));
    yy.deviations = malloc((x1->nTerms + x2->nTerms + 1) * sizeof(mpfr_t));

    iy = 0;
    ix1 = 0;
    ix2 = 0;
    x1HasNext = x1->nTerms > 0;
    x2HasNext = x2->nTerms > 0;
    while (x1HasNext || x2HasNext) {
        mpfr_init2(&(yy.deviations[iy]), prec_internal);

        if ((!x2HasNext) || (x1HasNext && (x1->symbols[ix1] < x2->symbols[ix2]))) {
            // y[i] = (x2[0] * x1[i])
            yy.symbols[iy] = x1->symbols[ix1];
            arpra_helper_mpfr_rnd_err_f2(error, &mpfr_mul, &(yy.deviations[iy]), &(x2->centre), &(x1->deviations[ix1]), MPFR_RNDN);
            x1HasNext = ++ix1 < x1->nTerms;
        }
        else if ((!x1HasNext) || (x2HasNext && (x2->symbols[ix2] < x1->symbols[ix1]))) {
            // y[i] = (x1[0] * x2[i])
            yy.symbols[iy] = x2->symbols[ix2];
            arpra_helper+mpfr_rnd_err_f2(error, &mpfr_mul, &(yy.deviations[iy]), &(x1->centre), &(x2->deviations[ix2]), MPFR_RNDN);
            x2HasNext = ++ix2 < x2->nTerms;
        }
        else {
            // y[i] = (x2[0] * x1[i]) + (x1[0] * x2[i])
            yy.symbols[iy] = x1->symbols[ix1];
            arpra_helper_mpfr_rnd_err_fmma(error, &(yy.deviations[iy]), &(x2->centre), &(x1->deviations[ix1]), &(x1->centre), &(x2->deviations[ix2]), MPFR_RNDN));
            x1HasNext = ++ix1 < x1->nTerms;
            x2HasNext = ++ix2 < x2->nTerms;
        }

        mpfr_abs(temp1, &(yy.deviations[iy]), MPFR_RNDU);
        mpfr_add(&(yy.radius), &(yy.radius), temp1, MPFR_RNDU);
        iy++;
    }

#ifdef ARPRA_TIGHT_MUL
    mpfr_t x1ix2j, x1jx2i, x1ix2i_pos_error, x1ix2i_neg_error;
    arpra_uint x1i_idx, x1j_idx, x2i_idx, x2j_idx;

    // If ARPRA_TIGHT_MUL is defined, then the linear approximation of the quadratic term of
    // arpra_mul (in mul.c) is defined the same as in (26) of:
    // S. M. Rump and M. Kashiwagi, Implementation and improvements of affine arithmetic,
    // Nonlinear Theory an Its Applications, IEICE, vol. 6, no. 3, pp. 341-359, 2015.

    // Init extra temp vars.
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
            mpfr_mul(temp1, &(x1->deviations[x1i_idx]), &(x2->deviations[x2i_idx]), MPFR_RNDA);
            if (mpfr_sgn(temp1) > 0) {
                mpfr_add(x1ix2i_pos_error, x1ix2i_pos_error, temp1, MPFR_RNDU);
            }
            else if (mpfr_sgn(temp1) < 0) {
                mpfr_sub(x1ix2i_neg_error, x1ix2i_neg_error, temp1, MPFR_RNDU);
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
                    mpfr_add(temp1, x1ix2j, x1jx2i, MPFR_RNDU);
                    if (mpfr_sgn(temp1) < 0) {
                        mpfr_mul(x1ix2j, &(x1->deviations[x1i_idx]), &(x2->deviations[x2j_idx]), MPFR_RNDD);
                        mpfr_mul(x1jx2i, &(x1->deviations[x1j_idx]), &(x2->deviations[x2i_idx]), MPFR_RNDD);
                        mpfr_add(temp1, x1ix2j, x1jx2i, MPFR_RNDD);
                    }
                    mpfr_abs(temp1, temp1, MPFR_RNDU);
                    mpfr_add(error, error, temp1, MPFR_RNDU);
                    x1HasNext = ++x1j_idx < x1->nTerms;
                    x2HasNext = ++x2j_idx < x2->nTerms;
                }
            }
            x1i_idx++;
            x2i_idx++;
        }
    }

    mpfr_max(temp1, x1ix2i_pos_error, x1ix2i_neg_error, MPFR_RNDU);
    mpfr_add(error, error, temp1, MPFR_RNDU);

    // Clear extra temp vars.
    mpfr_clear(x1ix2j);
    mpfr_clear(x1jx2i);
    mpfr_clear(x1ix2i_pos_error);
    mpfr_clear(x1ix2i_neg_error);
#else
    // Trivial linear approximation of quadratic term is rad(x1) * rad(x2).
    mpfr_mul(temp1, &(x1->radius), &(x2->radius), MPFR_RNDU);
    mpfr_add(error, error, temp1, MPFR_RNDU);
#endif

    // Store numerical error term.
    yy.symbols[iy] = arpra_helper_next_symbol();
    yy.deviations[iy] = *error;
    mpfr_add(&(yy.radius), &(yy.radius), &(yy.deviations[iy]), MPFR_RNDU);
    yy.nTerms = iy + 1;

    // Compute true_range, and add rounding error.
    arpra_helper_range_rounded(&yy);

#ifdef ARPRA_MIXED_IAAA
    // Intersect AA and IA ranges.
    mpfi_intersect(&(yy.true_range), &(yy.true_range), ia_range);

#ifdef ARPRA_MIXED_TRIMMED_IAAA
    // Trim error term if AA range fully encloses mixed IA/AA range.
    mpfr_sub(temp1, &(yy.centre), &(yy.radius), MPFR_RNDD);
    mpfr_add(temp2, &(yy.centre), &(yy.radius), MPFR_RNDU);
    if (mpfr_less_p(temp1, &(yy.true_range.left))
        && mpfr_greater_p(temp2, &(yy.true_range.right))) {
        mpfr_sub(temp1, &(yy.true_range.left), temp1, MPFR_RNDD);
        mpfr_sub(temp2, temp2, &(yy.true_range.right), MPFR_RNDD);
        mpfr_min(temp1, temp1, temp2, MPFR_RNDD);
        if (mpfr_greater_p(temp1, &(yy.deviations[yy.nTerms - 1]))) {
            mpfr_sub(&(yy.radius), &(yy.radius), &(yy.deviations[yy.nTerms - 1]), MPFR_RNDU);
            mpfr_set_zero(&(yy.deviations[yy.nTerms - 1]), 1);
        }
        else {
            mpfr_sub(&(yy.radius), &(yy.radius), temp1, MPFR_RNDU);
            mpfr_sub(&(yy.deviations[yy.nTerms - 1]), &(yy.deviations[yy.nTerms - 1]), temp1, MPFR_RNDU);
        }
    }
#endif // ARPRA_MIXED_TRIMMED_IAAA
#endif // ARPRA_MIXED_IAAA

    // Check for NaN and Inf.
    arpra_helper_check_result(&yy);

    // Clear vars, and set y.
    mpfi_clear(ia_range);
    mpfr_clear(temp1);
    mpfr_clear(temp2);
    arpra_clear(y);
    *y = yy;
}
