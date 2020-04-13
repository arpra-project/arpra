/*
 * sum.c -- Sum an array of Arpra ranges.
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

void arpra_sum (arpra_range *y, arpra_range *x, arpra_uint n)
{
    mpfr_t temp1, temp2, error;
    mpfr_ptr *summands;
    arpra_range yy;
    arpra_prec prec_internal;
    arpra_uint i, n_sum;
    arpra_uint i_y, *i_x;
    arpra_uint symbol;
    arpra_int xHasNext;

    // Handle n <= 2 case.
    if (n <= 2) {
        if (n == 2) {
            arpra_add(y, &x[0], &x[1]);
        }
        else if (n == 1) {
            arpra_set(y, &x[0]);
        }
        else {
            arpra_set_nan(y);
        }
        return;
    }

    // Domain violations:
    // (NaN) + ... + (NaN) = (NaN)
    // (NaN) + ... + (R)   = (NaN)
    // (Inf) + ... + (Inf) = (NaN)
    // (Inf) + ... + (R)   = (Inf)

    // Handle domain violations.
    for (i = 0; i < n; i++) {
        if (arpra_nan_p(&x[i])) {
            arpra_set_nan(y);
            return;
        }
    }
    for (i = 0; i < n; i++) {
        if (arpra_inf_p(&x[i])) {
            for (++i; i < n; i++) {
                if (arpra_inf_p(&x[i])) {
                    arpra_set_nan(y);
                    return;
                }
            }
            arpra_set_inf(y);
            return;
        }
    }

    // Initialise vars.
    prec_internal = arpra_get_internal_precision();
    mpfr_init2(temp1, prec_internal + 8);
    mpfr_init2(temp2, prec_internal + 8);
    mpfr_init2(error, prec_internal);
    arpra_init2(&yy, y->precision);
    summands = malloc(n * sizeof(mpfr_ptr));
    i_x = malloc(n * sizeof(arpra_uint));
    mpfr_set_zero(error, 1);

    // Zero term indexes, and fill summand array with centre values.
    i_y = 0;
    for (i = 0; i < n; i++) {
        i_x[i] = 0;
        summands[i] = &(x[i].centre);
    }

    // y[0] = x1[0] + ... + xn[0]
    ARPRA_MPFR_RNDERR_SUM(error, MPFR_RNDN, &(yy.centre), summands, n);

    // Allocate memory for deviation terms.
    yy.nTerms = 1;
    for (i = 0; i < n; i++) {
        yy.nTerms += x[i].nTerms;
    }
    yy.symbols = malloc(yy.nTerms * sizeof(arpra_uint));
    yy.deviations = malloc(yy.nTerms * sizeof(mpfr_t));

    // For all unique symbols in x.
    xHasNext = yy.nTerms > 1;
    while (xHasNext) {
        mpfr_init2(&(yy.deviations[i_y]), prec_internal);
        xHasNext = 0;
        symbol = -1;

        // Find and set the next lowest symbol in y.
        for (i = 0; i < n; i++) {
            if (i_x[i] < x[i].nTerms) {
                if (x[i].symbols[i_x[i]] < symbol) {
                    symbol = x[i].symbols[i_x[i]];
                }
            }
        }
        yy.symbols[i_y] = symbol;

        // For all x with the next symbol:
        for (n_sum = 0, i = 0; i < n; i++) {
            if (x[i].symbols[i_x[i]] == symbol) {
                // Get next deviation pointer of x[i].
                summands[n_sum++] = &(x[i].deviations[i_x[i]]);
                xHasNext += ++i_x[i] < x[i].nTerms;
            }
        }

        // y[i] = x1[i] + ... + xn[i]
        ARPRA_MPFR_RNDERR_SUM(error, MPFR_RNDN, &(yy.deviations[i_y]), summands, n_sum);
        i_y++;
    }

    // Store new deviation term.
    yy.symbols[i_y] = arpra_helper_next_symbol();
    yy.deviations[i_y] = *error;
    yy.nTerms = i_y + 1;

    // Compute true_range.
    arpra_helper_compute_range(&yy);

    // Check for NaN and Inf.
    arpra_helper_check_result(&yy);

    // Clear vars, and set y.
    mpfr_clear(temp1);
    mpfr_clear(temp2);
    arpra_clear(y);
    *y = yy;
    free(summands);
    free(i_x);
}

/*
 * Error bound for recursive summation (any ordering, any n).
 *
 * Rump, S. Error estimation of floating-point summation and dot product.
 * BIT Numer Math (2012) 52:201-220.
 */

void arpra_sum_recursive (arpra_range *y, arpra_range *x, arpra_uint n)
{
    mpfr_t temp1, temp2;
    mpfr_ptr sum_x, *sum_x_ptr;
    arpra_prec prec_internal;
    arpra_uint i;

    // Handle n <= 2 case.
    if (n <= 2) {
        if (n == 2) {
            arpra_add(y, &x[0], &x[1]);
        }
        else if (n == 1) {
            arpra_set(y, &x[0]);
        }
        else {
            arpra_set_nan(y);
        }
        return;
    }

    // Domain violations:
    // (NaN) + ... + (NaN) = (NaN)
    // (NaN) + ... + (R)   = (NaN)
    // (Inf) + ... + (Inf) = (NaN)
    // (Inf) + ... + (R)   = (Inf)

    // Handle domain violations.
    for (i = 0; i < n; i++) {
        if (arpra_nan_p(&x[i])) {
            arpra_set_nan(y);
            return;
        }
    }
    for (i = 0; i < n; i++) {
        if (arpra_inf_p(&x[i])) {
            for (++i; i < n; i++) {
                if (arpra_inf_p(&x[i])) {
                    arpra_set_nan(y);
                    return;
                }
            }
            arpra_set_inf(y);
            return;
        }
    }

    // Initialise vars.
    prec_internal = arpra_get_internal_precision();
    mpfr_init2(temp1, prec_internal);
    mpfr_init2(temp2, prec_internal);
    sum_x = malloc(n * sizeof(mpfr_t));
    sum_x_ptr = malloc(n * sizeof(mpfr_t));

    // Compute |x|.
    for (i = 0; i < n; i++) {
        if (mpfr_cmpabs(&(x[i]->true_range.right), &(x[i]->true_range.left)) >= 0) {
            sum_x[i] = x[i]->true_range.right;
        }
        else {
            sum_x[i] = x[i]->true_range.left;
            sum_x[i]._mpfr_sign = 1;
        }
        sum_x_ptr[i] = &(sum_x[i]);
    }

    // Compute error(sum_recursive(x)) = (n - 1) u sum(|x|).
    mpfr_set_si_2exp(temp1, 1, -y->precision, MPFR_RNDU);
    mpfr_mul_ui(temp1, temp1, (n - 1), MPFR_RNDU);
    mpfr_sum(temp2, sum_x_ptr, n, MPFR_RNDU);
    mpfr_mul(temp1, temp1, temp2, MPFR_RNDU);

    // Sum x, then add recursive sum error.
    arpra_sum(y, x, n);
    arpra_increase(y, y, temp1);

    // Clear vars.
    mpfr_clear(temp1);
    mpfr_clear(temp2);
    free(sum_x);
    free(sum_x_ptr);
}
