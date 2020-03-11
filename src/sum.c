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
    arpra_uint iy, *ix;
    arpra_uint symbol;
    arpra_int xHasNext;

    // Handle n <= 2 cases.
    if (n <= 2) {
        if (n == 2) {
            arpra_add(y, &x[0], &x[1]);
        }
        else if (n == 1) {
            arpra_set_range(y, &x[0]);
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
    ix = malloc(n * sizeof(arpra_uint));
    mpfr_set_zero(error, 1);
    mpfr_set_zero(&(yy.radius), 1);

    // Zero term indexes, and fill summand array with centre values.
    iy = 0;
    for (i = 0; i < n; i++) {
        ix[i] = 0;
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
        mpfr_init2(&(yy.deviations[iy]), prec_internal);
        xHasNext = 0;
        symbol = -1;

        // Find and set the next lowest symbol in y.
        for (i = 0; i < n; i++) {
            if (ix[i] < x[i].nTerms) {
                if (x[i].symbols[ix[i]] < symbol) {
                    symbol = x[i].symbols[ix[i]];
                }
            }
        }
        yy.symbols[iy] = symbol;

        // For all x with the next symbol:
        for (n_sum = 0, i = 0; i < n; i++) {
            if (x[i].symbols[ix[i]] == symbol) {
                // Get next deviation pointer of x[i].
                summands[n_sum++] = &(x[i].deviations[ix[i]]);

                xHasNext += ++ix[i] < x[i].nTerms;
            }
        }

        // y[i] = x1[i] + ... + xn[i]
        ARPRA_MPFR_RNDERR_SUM(error, MPFR_RNDN, &(yy.deviations[iy]), summands, n_sum);

        mpfr_abs(temp1, &(yy.deviations[iy]), MPFR_RNDU);
        mpfr_add(&(yy.radius), &(yy.radius), temp1, MPFR_RNDU);
        iy++;
    }

    // Store new deviation term.
    yy.symbols[iy] = arpra_helper_next_symbol();
    yy.deviations[iy] = *error;
    mpfr_add(&(yy.radius), &(yy.radius), &(yy.deviations[iy]), MPFR_RNDU);
    yy.nTerms = iy + 1;

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
    free(ix);
}

void arpra_sum_recursive (arpra_range *y, arpra_range *x, arpra_uint n)
{
    mpfr_t temp1, temp2, error;
    mpfr_ptr *summands;
    arpra_range yy;
    arpra_prec prec_internal;
    arpra_uint i, n_sum;
    arpra_uint iy, *ix;
    arpra_uint symbol;
    arpra_int xHasNext;

    // Handle n <= 2 cases.
    if (n <= 2) {
        if (n == 2) {
            arpra_add(y, &x[0], &x[1]);
        }
        else if (n == 1) {
            arpra_set_range(y, &x[0]);
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
    ix = malloc(n * sizeof(arpra_uint));
    mpfr_set_zero(error, 1);
    mpfr_set_zero(&(yy.radius), 1);

    // Zero term indexes, and fill summand array with centre values.
    iy = 0;
    for (i = 0; i < n; i++) {
        ix[i] = 0;
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
        mpfr_init2(&(yy.deviations[iy]), prec_internal);
        xHasNext = 0;
        symbol = -1;

        // Find and set the next lowest symbol in y.
        for (i = 0; i < n; i++) {
            if (ix[i] < x[i].nTerms) {
                if (x[i].symbols[ix[i]] < symbol) {
                    symbol = x[i].symbols[ix[i]];
                }
            }
        }
        yy.symbols[iy] = symbol;

        // For all x with the next symbol:
        for (n_sum = 0, i = 0; i < n; i++) {
            if (x[i].symbols[ix[i]] == symbol) {
                // Get next deviation pointer of x[i].
                summands[n_sum++] = &(x[i].deviations[ix[i]]);

                xHasNext += ++ix[i] < x[i].nTerms;
            }
        }

        // y[i] = x1[i] + ... + xn[i]
        ARPRA_MPFR_RNDERR_SUM(error, MPFR_RNDN, &(yy.deviations[iy]), summands, n_sum);

        mpfr_abs(temp1, &(yy.deviations[iy]), MPFR_RNDU);
        mpfr_add(&(yy.radius), &(yy.radius), temp1, MPFR_RNDU);
        iy++;
    }

    /*
     * Error bound for recursive summation (any ordering, any n).
     *
     * Rump, S. Error estimation of floating-point summation and dot product.
     * BIT Numer Math (2012) 52:201-220.
     */

    // Compute |x|.
    mpfr_ptr summands_magnitude;
    summands_magnitude = arpra_helper_buffer_mpfr(n);
    for (i = 0; i < n; i++) {
        mpfr_init2(&(summands_magnitude[i]), prec_internal);

        if (mpfr_sgn(&(x[i].centre)) >= 0) {
            mpfr_add(&(summands_magnitude[i]), &(x[i].centre), &(x[i].radius), MPFR_RNDU);
        }
        else {
            mpfr_sub(&(summands_magnitude[i]), &(x[i].centre), &(x[i].radius), MPFR_RNDD);
            mpfr_abs(&(summands_magnitude[i]), &(summands_magnitude[i]), MPFR_RNDU);
        }
    }

    // Compute error(sum(x)) = (n - 1) u sum(|x|).
    mpfr_set_si_2exp(temp1, 1, -yy.precision, MPFR_RNDU);
    mpfr_mul_ui(temp1, temp1, (n - 1), MPFR_RNDU);
    arpra_ext_mpfr_sum(&(summands_magnitude[0]), summands_magnitude, n, MPFR_RNDU);
    mpfr_mul(&(summands_magnitude[0]), &(summands_magnitude[0]), temp1, MPFR_RNDU);
    mpfr_add(error, error, &(summands_magnitude[0]), MPFR_RNDU);

    for (i = 0; i < n; i++) {
        mpfr_clear(&(summands_magnitude[i]));
    }

    /*
     * END: Error bound for recursive summation.
     */

    // Store new deviation term.
    yy.symbols[iy] = arpra_helper_next_symbol();
    yy.deviations[iy] = *error;
    mpfr_add(&(yy.radius), &(yy.radius), &(yy.deviations[iy]), MPFR_RNDU);
    yy.nTerms = iy + 1;

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
    free(ix);
}
