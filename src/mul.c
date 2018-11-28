/*
 * mul.c -- Multiply one Arpra range by another.
 *
 * Copyright 2016-2018 James Paul Turner.
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

void arpra_mul (arpra_range *z, const arpra_range *x, const arpra_range *y)
{
    arpra_uint xTerm, yTerm, zTerm;
    arpra_int xHasNext, yHasNext;
    arpra_mpfr temp, error;
    arpra_mpfi ia_range, temp_range;
    arpra_range zNew;
    arpra_prec prec_internal;

    // Domain violations:
    // NaN  *  NaN  =  NaN
    // NaN  *  r    =  NaN
    // Inf  *  0    =  NaN
    // Inf  *  Inf  =  Inf
    // Inf  *  s    =  Inf
    // s.t. (r in R) and (s in R \ {0})

    // Handle domain violations.
    if (arpra_nan_p(x) || arpra_nan_p(y)) {
        arpra_set_nan(z);
        return;
    }
    if (arpra_inf_p(x)) {
        if (arpra_has_zero_p(y)) {
            arpra_set_nan(z);
        }
        else {
            arpra_set_inf(z);
        }
        return;
    }
    if (arpra_inf_p(y)) {
        if (arpra_has_zero_p(x)) {
            arpra_set_nan(z);
        }
        else {
            arpra_set_inf(z);
        }
        return;
    }

    // Initialise vars.
    prec_internal = arpra_get_internal_precision();
    mpfr_init2(&temp, prec_internal);
    mpfr_init2(&error, prec_internal);
    mpfi_init2(&ia_range, z->precision);
    mpfi_init2(&temp_range, prec_internal);
    arpra_init2(&zNew, z->precision);
    mpfr_set_ui(&error, 0, MPFR_RNDU);
    mpfr_set_ui(&(zNew.radius), 0, MPFR_RNDU);

    // MPFI multiplication
    mpfi_mul(&ia_range, &(x->true_range), &(y->true_range));

    // z_0 = x_0 * y_0
    if (mpfr_mul(&(zNew.centre), &(x->centre), &(y->centre), MPFR_RNDN)) {
        arpra_helper_error_half_ulp(&temp, &(zNew.centre));
        mpfr_add(&error, &error, &temp, MPFR_RNDU);
    }

    // Allocate memory for all possible deviation terms.
    zNew.nTerms = x->nTerms + y->nTerms + 1;
    zNew.symbols = malloc(zNew.nTerms * sizeof(arpra_uint));
    zNew.deviations = malloc(zNew.nTerms * sizeof(arpra_mpfr));

    xTerm = 0;
    yTerm = 0;
    zTerm = 0;
    xHasNext = x->nTerms > 0;
    yHasNext = y->nTerms > 0;
    while (xHasNext || yHasNext) {
        if ((!yHasNext) || (xHasNext && (x->symbols[xTerm] < y->symbols[yTerm]))) {
            zNew.symbols[zTerm] = x->symbols[xTerm];
            mpfr_init2(&(zNew.deviations[zTerm]), prec_internal);

            // z_i = (y_0 * x_i)
            if (mpfr_mul(&(zNew.deviations[zTerm]), &(y->centre), &(x->deviations[xTerm]), MPFR_RNDN)) {
                arpra_helper_error_half_ulp(&temp, &(zNew.deviations[zTerm]));
                mpfr_add(&error, &error, &temp, MPFR_RNDU);
            }

            xHasNext = ++xTerm < x->nTerms;
        }
        else if ((!xHasNext) || (yHasNext && (y->symbols[yTerm] < x->symbols[xTerm]))) {
            zNew.symbols[zTerm] = y->symbols[yTerm];
            mpfr_init2(&(zNew.deviations[zTerm]), prec_internal);

            // z_i = (x_0 * y_i)
            if (mpfr_mul(&(zNew.deviations[zTerm]), &(x->centre), &(y->deviations[yTerm]), MPFR_RNDN)) {
                arpra_helper_error_half_ulp(&temp, &(zNew.deviations[zTerm]));
                mpfr_add(&error, &error, &temp, MPFR_RNDU);
            }

            yHasNext = ++yTerm < y->nTerms;
        }
        else {
            zNew.symbols[zTerm] = x->symbols[xTerm];
            mpfr_init2(&(zNew.deviations[zTerm]), prec_internal);

            // z_i = (y_0 * x_i) + (x_0 * y_i)
            if (arpra_helper_term(&(zNew.deviations[zTerm]), &(x->deviations[xTerm]), &(y->deviations[yTerm]), &(y->centre), &(x->centre), NULL)) {
                arpra_helper_error_half_ulp(&temp, &(zNew.deviations[zTerm]));
                mpfr_add(&error, &error, &temp, MPFR_RNDU);
            }

            xHasNext = ++xTerm < x->nTerms;
            yHasNext = ++yTerm < y->nTerms;
        }

        // Store nonzero deviation terms.
        if (mpfr_zero_p(&(zNew.deviations[zTerm]))) {
            mpfr_clear(&(zNew.deviations[zTerm]));
        }
        else {
            mpfr_abs(&temp, &(zNew.deviations[zTerm]), MPFR_RNDU);
            mpfr_add(&(zNew.radius), &(zNew.radius), &temp, MPFR_RNDU);
            zTerm++;
        }
    }

#ifdef ARPRA_TIGHT_MUL
    // Linear approximation of the quadratic term is defined the same as in (26) of:
    // S. M. Rump and M. Kashiwagi, Implementation and improvements of affine arithmetic,
    // Nonlinear Theory an Its Applications, IEICE, vol. 6, no. 3, pp. 341-359, 2015.

    arpra_uint xi_idx, xj_idx, yi_idx, yj_idx;
    arpra_mpfr xiyj, xjyi, xiyi_pos_error, xiyi_neg_error;

    // Init extra temp vars.
    mpfr_init2(&xiyj, prec_internal);
    mpfr_init2(&xjyi, prec_internal);
    mpfr_init2(&xiyi_pos_error, prec_internal);
    mpfr_init2(&xiyi_neg_error, prec_internal);
    mpfr_set_ui(&xiyi_pos_error, 0, MPFR_RNDN);
    mpfr_set_ui(&xiyi_neg_error, 0, MPFR_RNDN);

    xi_idx = 0;
    yi_idx = 0;
    while ((xi_idx < x->nTerms) && (yi_idx < y->nTerms)) {
        if (x->symbols[xi_idx] < y->symbols[yi_idx]) {
            for (yj_idx = yi_idx; yj_idx < y->nTerms; yj_idx++) {
                // x has symbol i, and y has symbol j, so error += abs(xi * yj)
                mpfr_mul(&xiyj, &(x->deviations[xi_idx]), &(y->deviations[yj_idx]), MPFR_RNDA);
                mpfr_abs(&xiyj, &xiyj, MPFR_RNDU);
                mpfr_add(&error, &error, &xiyj, MPFR_RNDU);
            }
            xi_idx++;
        }
        else if (y->symbols[yi_idx] < x->symbols[xi_idx]) {
            for (xj_idx = xi_idx; xj_idx < x->nTerms; xj_idx++) {
                // y has symbol i, and x has symbol j, so error += abs(xj * yi)
                mpfr_mul(&xjyi, &(x->deviations[xj_idx]), &(y->deviations[yi_idx]), MPFR_RNDA);
                mpfr_abs(&xjyi, &xjyi, MPFR_RNDU);
                mpfr_add(&error, &error, &xjyi, MPFR_RNDU);
            }
            yi_idx++;
        }
        else {
            // both x and y have symbol i, so error += abs(xi * yi)
            mpfr_mul(&temp, &(x->deviations[xi_idx]), &(y->deviations[yi_idx]), MPFR_RNDA);
            if (mpfr_sgn(&temp) > 0) {
                mpfr_add(&xiyi_pos_error, &xiyi_pos_error, &temp, MPFR_RNDU);
            }
            else if (mpfr_sgn(&temp) < 0) {
                mpfr_sub(&xiyi_neg_error, &xiyi_neg_error, &temp, MPFR_RNDU);
            }

            xj_idx = xi_idx + 1;
            yj_idx = yi_idx + 1;
            xHasNext = xj_idx < x->nTerms;
            yHasNext = yj_idx < y->nTerms;
            while (xHasNext || yHasNext) {
                if ((!yHasNext) || (xHasNext && (x->symbols[xi_idx] < y->symbols[yi_idx]))) {
                    // both x and y have symbol i, but only x has symbol j, so error += abs(xj * yi)
                    mpfr_mul(&xjyi, &(x->deviations[xj_idx]), &(y->deviations[yi_idx]), MPFR_RNDA);
                    mpfr_abs(&xjyi, &xjyi, MPFR_RNDU);
                    mpfr_add(&error, &error, &xjyi, MPFR_RNDU);

                    xHasNext = ++xj_idx < x->nTerms;
                }
                else if ((!xHasNext) || (yHasNext && (y->symbols[yi_idx] < x->symbols[xi_idx]))) {
                    // both x and y have symbol i, but only y has symbol j, so error += abs(xi * yj)
                    mpfr_mul(&xiyj, &(x->deviations[xi_idx]), &(y->deviations[yj_idx]), MPFR_RNDA);
                    mpfr_abs(&xiyj, &xiyj, MPFR_RNDU);
                    mpfr_add(&error, &error, &xiyj, MPFR_RNDU);

                    yHasNext = ++yj_idx < y->nTerms;
                }
                else {
                    // both x and y have symbols i and j, so error += abs(xi * yj + xj * yi)
                    mpfr_mul(&xiyj, &(x->deviations[xi_idx]), &(y->deviations[yj_idx]), MPFR_RNDU);
                    mpfr_mul(&xjyi, &(x->deviations[xj_idx]), &(y->deviations[yi_idx]), MPFR_RNDU);
                    mpfr_add(&temp, &xiyj, &xjyi, MPFR_RNDU);
                    if (mpfr_sgn(&temp) < 0) {
                        mpfr_mul(&xiyj, &(x->deviations[xi_idx]), &(y->deviations[yj_idx]), MPFR_RNDD);
                        mpfr_mul(&xjyi, &(x->deviations[xj_idx]), &(y->deviations[yi_idx]), MPFR_RNDD);
                        mpfr_add(&temp, &xiyj, &xjyi, MPFR_RNDD);
                    }
                    mpfr_abs(&temp, &temp, MPFR_RNDU);
                    mpfr_add(&error, &error, &temp, MPFR_RNDU);

                    xHasNext = ++xj_idx < x->nTerms;
                    yHasNext = ++yj_idx < y->nTerms;
                }
            }
            xi_idx++;
            yi_idx++;
        }
    }

    mpfr_max(&temp, &xiyi_pos_error, &xiyi_neg_error, MPFR_RNDU);
    mpfr_add(&error, &error, &temp, MPFR_RNDU);

    // Clear extra temp vars.
    mpfr_clear(&xiyj);
    mpfr_clear(&xjyi);
    mpfr_clear(&xiyi_pos_error);
    mpfr_clear(&xiyi_neg_error);
#else
    // Linear approximation of quadratic term is rad(x) * rad(y).
    mpfr_mul(&temp, &(x->radius), &(y->radius), MPFR_RNDU);
    mpfr_add(&error, &error, &temp, MPFR_RNDU);
#endif

    // Compute target precision rounding error.
    mpfr_sub(&(temp_range.left), &(zNew.centre), &(zNew.radius), MPFR_RNDD);
    mpfr_sub(&(temp_range.left), &(temp_range.left), &error, MPFR_RNDD);
    mpfr_add(&(temp_range.right), &(zNew.centre), &(zNew.radius), MPFR_RNDU);
    mpfr_add(&(temp_range.right), &(temp_range.right), &error, MPFR_RNDU);
    mpfi_set(&(zNew.true_range), &temp_range);

    mpfr_sub(&(temp_range.left), &(temp_range.left), &(zNew.true_range.left), MPFR_RNDU);
    mpfr_sub(&(temp_range.right), &(zNew.true_range.right), &(temp_range.right), MPFR_RNDU);
    mpfr_max(&temp, &(temp_range.left), &(temp_range.right), MPFR_RNDU);
    mpfr_add(&error, &error, &temp, MPFR_RNDU);

    // Store numerical error term.
    zNew.symbols[zTerm] = arpra_next_symbol();
    zNew.deviations[zTerm] = error;
    mpfr_add(&(zNew.radius), &(zNew.radius), &(zNew.deviations[zTerm]), MPFR_RNDU);
    zTerm++;

    // Compute true range.
    mpfr_sub(&(zNew.true_range.left), &(zNew.centre), &(zNew.radius), MPFR_RNDD);
    mpfr_add(&(zNew.true_range.right), &(zNew.centre), &(zNew.radius), MPFR_RNDU);
    mpfi_intersect(&(zNew.true_range), &(zNew.true_range), &ia_range);

    // Handle domain violations, and free unused memory.
    zNew.nTerms = zTerm;
    if (mpfr_nan_p(&(zNew.centre)) || mpfr_nan_p(&(zNew.radius))) {
        arpra_set_nan(&zNew);
    }
    else if (mpfr_inf_p(&(zNew.centre)) || mpfr_inf_p(&(zNew.radius))) {
        arpra_set_inf(&zNew);
    }
    else {
        if (zNew.nTerms == 0) {
            free(zNew.symbols);
            free(zNew.deviations);
        }
    }

    // Clear vars, and set z.
    mpfr_clear(&temp);
    mpfi_clear(&ia_range);
    mpfi_clear(&temp_range);
    arpra_clear(z);
    *z = zNew;
}
