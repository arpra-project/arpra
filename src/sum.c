/*
 * sum.c -- Sum an array of affine forms.
 *
 * Copyright 2016-2017 James Paul Turner.
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

#include "mpfa-impl.h"

void mpfa_sum (mpfa_ptr z, const mpfa_ptr *x, mpfa_uint_t n)
{
    mpfa_uint_t i, j;
    mpfa_uint_t xSymbol, zTerm;
    mpfa_uint_t *xTerm;
    mpfa_int_t xHasNext;
    mpfr_ptr *summands;
    mpfr_t temp, error;
    mpfa_prec_t prec, prec_internal;
    mpfa_t zNew;

    // Handle trivial cases.
    if (n <= 2) {
        if (n == 2) {
            mpfa_add(z, x[0], x[1]);
            return;
        }
        else if (n == 1) {
            mpfa_set(z, x[0]);
            return;
        }
        else {
            mpfa_set_none(z);
            return;
        }
    }

    // Handle domain violations.
    for (i = 0; i < n; i++) {
        if (mpfa_none_p(x[i])) {
            mpfa_set_none(z);
            return;
        }
    }
    for (i = 0; i < n; i++) {
        if (mpfa_any_p(x[i])) {
            mpfa_set_any(z);
            return;
        }
    }

    // Init temp vars, and set internal precision.
    prec = mpfa_get_prec(z);
    prec_internal = mpfa_get_internal_prec();
    mpfr_init2(temp, prec_internal);
    mpfr_init2(error, prec_internal);
    mpfr_init2(&(zNew->centre), prec);
    mpfr_init2(&(zNew->radius), prec_internal);
    mpfr_set_si(error, 0, MPFR_RNDU);
    mpfr_set_si(&(zNew->radius), 0, MPFR_RNDU);
    xTerm = malloc(n * sizeof(mpfa_uint_t));
    summands = malloc(n * sizeof(mpfr_ptr));

    // Zero term indexes, and fill summand array with centre values.
    zTerm = 0;
    for (i = 0; i < n; i++) {
        xTerm[i] = 0;
        summands[i] = &(x[i]->centre);
    }

    // z_0 = x[1]_0 + ... + x[n]_0
    if (mpfr_sum(&(zNew->centre), summands, n, MPFR_RNDN)) {
        mpfa_error(temp, &(zNew->centre));
        mpfr_add(error, error, temp, MPFR_RNDU);
    }

    // Allocate memory for all possible noise terms.
    zNew->nTerms = 1;
    for (i = 0; i < n; i++) {
        zNew->nTerms += x[i]->nTerms;
    }
    zNew->symbols = malloc(zNew->nTerms * sizeof(mpfa_uint_t));
    zNew->deviations = malloc(zNew->nTerms * sizeof(mpfr_t));

    // For all unique symbols in x.
    xHasNext = zNew->nTerms > 1;
    while (xHasNext) {
        xHasNext = 0;
        xSymbol = -1;

        // Find and set the next symbol in z.
        for (i = 0; i < n; i++) {
            if (xTerm[i] < x[i]->nTerms) {
                if (x[i]->symbols[xTerm[i]] < xSymbol) {
                    xSymbol = x[i]->symbols[xTerm[i]];
                }
            }
        }
        zNew->symbols[zTerm] = xSymbol;
        mpfr_init2(&(zNew->deviations[zTerm]), prec);

        // For all x with the next symbol:
        for (i = 0, j = 0; i < n; i++) {
            if (x[i]->symbols[xTerm[i]] == xSymbol) {
                // Get next deviation pointer of x[i].
                summands[j++] = &(x[i]->deviations[xTerm[i]]);

                // Check for more symbols in x[i].
                if (++xTerm[i] < x[i]->nTerms) {
                    xHasNext = 1;
                }
            }
        }

        // z_i = x[1]_i + ... + x[n]_i
        if (mpfr_sum(&(zNew->deviations[zTerm]), summands, j, MPFR_RNDN)) {
            mpfa_error(temp, &(zNew->deviations[zTerm]));
            mpfr_add(error, error, temp, MPFR_RNDU);
        }

        // Store nonzero noise terms.
        if (mpfr_zero_p(&(zNew->deviations[zTerm]))) {
            mpfr_clear(&(zNew->deviations[zTerm]));
        }
        else {
            mpfr_abs(temp, &(zNew->deviations[zTerm]), MPFR_RNDU);
            mpfr_add(&(zNew->radius), &(zNew->radius), temp, MPFR_RNDU);
            zTerm++;
        }
    }

    // Store nonzero numerical error term.
    if (!mpfr_zero_p(error)) {
        zNew->symbols[zTerm] = mpfa_next_sym();
        mpfr_init2(&(zNew->deviations[zTerm]), prec);
        mpfr_set(&(zNew->deviations[zTerm]), error, MPFR_RNDU);
        mpfr_add(&(zNew->radius), &(zNew->radius), &(zNew->deviations[zTerm]), MPFR_RNDU);
        zTerm++;
    }

    // Handle domain violations, and free unused memory.
    zNew->nTerms = zTerm;
    if (mpfr_nan_p(&(zNew->centre)) || mpfr_nan_p(&(zNew->radius))) {
        mpfa_set_none(zNew);
    }
    else if (mpfr_inf_p(&(zNew->centre)) || mpfr_inf_p(&(zNew->radius))) {
        mpfa_set_any(zNew);
    }
    else {
        if (zNew->nTerms == 0) {
            free(zNew->symbols);
            free(zNew->deviations);
        }
    }

    // Clear temp vars, and set z.
    mpfr_clear(temp);
    mpfr_clear(error);
    mpfa_set(z, zNew);
    mpfa_clear(zNew);
    free(xTerm);
    free(summands);
}
