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
    arpra_mpfr temp1, temp2, error;
    arpra_mpfi ia_range;
    arpra_range zNew;
    arpra_prec prec_internal;

    // Domain violations:
    // NaN  *  NaN  =  NaN
    // NaN  *  r    =  NaN
    // Inf  *  0    =  NaN
    // Inf  *  Inf  =  Inf
    // Inf  *  r    =  Inf

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
    mpfr_init2(&temp1, prec_internal + 8);
    mpfr_init2(&temp2, prec_internal + 8);
    mpfr_init2(&error, prec_internal);
    mpfi_init2(&ia_range, z->precision);
    arpra_init2(&zNew, z->precision);
    mpfr_set_ui(&error, 0, MPFR_RNDU);
    mpfr_set_ui(&(zNew.radius), 0, MPFR_RNDU);

    // MPFI multiplication
    mpfi_mul(&ia_range, &(x->true_range), &(y->true_range));

    // z_0 = x_0 * y_0
    if (mpfr_mul(&(zNew.centre), &(x->centre), &(y->centre), MPFR_RNDN)) {
        arpra_helper_error_half_ulp(&temp1, &(zNew.centre));
        mpfr_add(&error, &error, &temp1, MPFR_RNDU);
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
                arpra_helper_error_half_ulp(&temp1, &(zNew.deviations[zTerm]));
                mpfr_add(&error, &error, &temp1, MPFR_RNDU);
            }

            xHasNext = ++xTerm < x->nTerms;
        }
        else if ((!xHasNext) || (yHasNext && (y->symbols[yTerm] < x->symbols[xTerm]))) {
            zNew.symbols[zTerm] = y->symbols[yTerm];
            mpfr_init2(&(zNew.deviations[zTerm]), prec_internal);

            // z_i = (x_0 * y_i)
            if (mpfr_mul(&(zNew.deviations[zTerm]), &(x->centre), &(y->deviations[yTerm]), MPFR_RNDN)) {
                arpra_helper_error_half_ulp(&temp1, &(zNew.deviations[zTerm]));
                mpfr_add(&error, &error, &temp1, MPFR_RNDU);
            }

            yHasNext = ++yTerm < y->nTerms;
        }
        else {
            zNew.symbols[zTerm] = x->symbols[xTerm];
            mpfr_init2(&(zNew.deviations[zTerm]), prec_internal);

            // z_i = (y_0 * x_i) + (x_0 * y_i)
            if (arpra_helper_term(&(zNew.deviations[zTerm]), &(x->deviations[xTerm]), &(y->deviations[yTerm]), &(y->centre), &(x->centre), NULL)) {
                arpra_helper_error_half_ulp(&temp1, &(zNew.deviations[zTerm]));
                mpfr_add(&error, &error, &temp1, MPFR_RNDU);
            }

            xHasNext = ++xTerm < x->nTerms;
            yHasNext = ++yTerm < y->nTerms;
        }

        mpfr_abs(&temp1, &(zNew.deviations[zTerm]), MPFR_RNDU);
        mpfr_add(&(zNew.radius), &(zNew.radius), &temp1, MPFR_RNDU);
        zTerm++;
    }

#ifdef ARPRA_TIGHT_MUL
    arpra_uint xi_idx, xj_idx, yi_idx, yj_idx;
    arpra_mpfr xiyj, xjyi, xiyi_pos_error, xiyi_neg_error;

    // If ARPRA_TIGHT_MUL is defined, then the linear approximation of the quadratic term of
    // arpra_mul (in mul.c) is defined the same as in (26) of:
    // S. M. Rump and M. Kashiwagi, Implementation and improvements of affine arithmetic,
    // Nonlinear Theory an Its Applications, IEICE, vol. 6, no. 3, pp. 341-359, 2015.

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
            mpfr_mul(&temp1, &(x->deviations[xi_idx]), &(y->deviations[yi_idx]), MPFR_RNDA);
            if (mpfr_sgn(&temp1) > 0) {
                mpfr_add(&xiyi_pos_error, &xiyi_pos_error, &temp1, MPFR_RNDU);
            }
            else if (mpfr_sgn(&temp1) < 0) {
                mpfr_sub(&xiyi_neg_error, &xiyi_neg_error, &temp1, MPFR_RNDU);
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
                    mpfr_add(&temp1, &xiyj, &xjyi, MPFR_RNDU);
                    if (mpfr_sgn(&temp1) < 0) {
                        mpfr_mul(&xiyj, &(x->deviations[xi_idx]), &(y->deviations[yj_idx]), MPFR_RNDD);
                        mpfr_mul(&xjyi, &(x->deviations[xj_idx]), &(y->deviations[yi_idx]), MPFR_RNDD);
                        mpfr_add(&temp1, &xiyj, &xjyi, MPFR_RNDD);
                    }
                    mpfr_abs(&temp1, &temp1, MPFR_RNDU);
                    mpfr_add(&error, &error, &temp1, MPFR_RNDU);

                    xHasNext = ++xj_idx < x->nTerms;
                    yHasNext = ++yj_idx < y->nTerms;
                }
            }
            xi_idx++;
            yi_idx++;
        }
    }

    mpfr_max(&temp1, &xiyi_pos_error, &xiyi_neg_error, MPFR_RNDU);
    mpfr_add(&error, &error, &temp1, MPFR_RNDU);

    // Clear extra temp vars.
    mpfr_clear(&xiyj);
    mpfr_clear(&xjyi);
    mpfr_clear(&xiyi_pos_error);
    mpfr_clear(&xiyi_neg_error);
#else
    // Trivial linear approximation of quadratic term is rad(x) * rad(y).
    mpfr_mul(&temp1, &(x->radius), &(y->radius), MPFR_RNDU);
    mpfr_add(&error, &error, &temp1, MPFR_RNDU);
#endif

    // Round range to target precision.
    mpfr_sub(&temp1, &(zNew.centre), &(zNew.radius), MPFR_RNDD);
    mpfr_sub(&temp1, &temp1, &error, MPFR_RNDD);
    mpfr_set(&(zNew.true_range.left), &temp1, MPFR_RNDD);
    mpfr_sub(&temp1, &temp1, &(zNew.true_range.left), MPFR_RNDU);
    mpfr_add(&temp2, &(zNew.centre), &(zNew.radius), MPFR_RNDU);
    mpfr_add(&temp2, &temp2, &error, MPFR_RNDU);
    mpfr_set(&(zNew.true_range.right), &temp2, MPFR_RNDU);
    mpfr_sub(&temp2, &(zNew.true_range.right), &temp2, MPFR_RNDU);
    mpfr_max(&temp1, &temp1, &temp2, MPFR_RNDU);
    mpfr_add(&error, &error, &temp1, MPFR_RNDU);

    // Store numerical error term.
    zNew.symbols[zTerm] = arpra_next_symbol();
    zNew.deviations[zTerm] = error;
    mpfr_add(&(zNew.radius), &(zNew.radius), &(zNew.deviations[zTerm]), MPFR_RNDU);
    zNew.nTerms = zTerm + 1;

#ifdef ARPRA_MIXED_IAAA
#ifdef ARPRA_MIXED_TRIMMED_IAAA
    // Trim error term if Arpra range fully contains IA range.
    if (mpfr_less_p(&(zNew.true_range.left), &(ia_range.left))
        && mpfr_greater_p(&(zNew.true_range.right), &(ia_range.right))) {
        mpfr_sub(&temp1, &(ia_range.left), &(zNew.true_range.left), MPFR_RNDD);
        mpfr_sub(&temp2, &(zNew.true_range.right), &(ia_range.right), MPFR_RNDD);
        mpfr_min(&temp1, &temp1, &temp2, MPFR_RNDD);
        mpfr_sub(&(zNew.deviations[zTerm]), &(zNew.deviations[zTerm]), &temp1, MPFR_RNDU);
        if (mpfr_cmp_ui(&(zNew.deviations[zTerm]), 0) < 0) {
            mpfr_set_ui(&(zNew.deviations[zTerm]), 0, MPFR_RNDN);
        }
    }
#endif // ARPRA_MIXED_TRIMMED_IAAA
    mpfi_intersect(&(zNew.true_range), &(zNew.true_range), &ia_range);
#endif // ARPRA_MIXED_IAAA

    // Handle domain violations.
    if (mpfr_nan_p(&(zNew.centre)) || mpfr_nan_p(&(zNew.radius))) {
        arpra_set_nan(&zNew);
    }
    else if (mpfr_inf_p(&(zNew.centre)) || mpfr_inf_p(&(zNew.radius))) {
        arpra_set_inf(&zNew);
    }

    // Clear vars, and set z.
    mpfr_clear(&temp1);
    mpfr_clear(&temp2);
    mpfi_clear(&ia_range);
    arpra_clear(z);
    *z = zNew;
}
