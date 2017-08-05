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

#include "mpfa.h"
#include <stdlib.h>

void mpfa_sum (mpfa_ptr z, const mpfa_ptr *x, mpfa_uint_t n) {
    mpfa_uint_t i, j;
    mpfa_uint_t xSymbol, zTerm;
    mpfa_uint_t *xTerm;
    mpfa_int_t xHasNext;
    mpfr_ptr *summands;
    mpfr_t temp, error;
    mpfa_prec_t prec, prec_internal;
    mpfa_t zNew;

    if (n < 3) {
        if (n == 2) {
            mpfa_add(z, x[0], x[1]);
        }
        else if (n == 1) {
            mpfa_set(z, x[0]);
        }
        else return;
    }

    prec = mpfa_get_prec(z);
    prec_internal = mpfa_get_internal_prec();
    mpfr_init2(temp, prec_internal);
    mpfr_init2(error, prec_internal);
    mpfr_init2(&(zNew->centre), prec);
    mpfr_init2(&(zNew->radius), prec_internal);
    mpfr_set_si(error, 0, MPFR_RNDN);
    mpfr_set_si(&(zNew->radius), 0, MPFR_RNDN);
    xTerm = malloc(n * sizeof(mpfa_uint_t));
    summands = malloc(n * sizeof(mpfr_ptr));

    zTerm = 0;
    for (i = 0; i < n; i++) {
        xTerm[i] = 0;
        summands[i] = &(x[i]->centre);
    }

    if (mpfr_sum(&(zNew->centre), summands, n, MPFR_RNDN)) {
        mpfa_error(temp, &(zNew->centre));
        mpfr_add(error, error, temp, MPFR_RNDU);
    }

    zNew->nTerms = 1;
    for (i = 0; i < n; i++) {
        zNew->nTerms += x[i]->nTerms;
    }
    zNew->symbols = malloc(zNew->nTerms * sizeof(mpfa_uint_t));
    zNew->deviations = malloc(zNew->nTerms * sizeof(mpfr_t));

    xHasNext = zNew->nTerms > 1;
    while (xHasNext) {
        xHasNext = 0;
        xSymbol = -1;

        for (i = 0; i < n; i++) {
            if (xTerm[i] < x[i]->nTerms) {
                if (x[i]->symbols[xTerm[i]] < xSymbol) {
                    xSymbol = x[i]->symbols[xTerm[i]];
                }
            }
        }

        zNew->symbols[zTerm] = xSymbol;
        mpfr_init2(&(zNew->deviations[zTerm]), prec_internal);

        for (i = 0, j = 0; i < n; i++) {
            if (x[i]->symbols[xTerm[i]] == xSymbol) {
                summands[j++] = &(x[i]->deviations[xTerm[i]]);

                if (mpfr_sum(&(zNew->deviations[zTerm]), summands, j, MPFR_RNDN)) {
                    mpfa_error(temp, &(zNew->deviations[zTerm]));
                    mpfr_add(error, error, temp, MPFR_RNDU);
                }

                if (++xTerm[i] < x[i]->nTerms) {
                    xHasNext = 1;
                }
            }
        }

        if (mpfr_zero_p(&(zNew->deviations[zTerm]))) {
            mpfr_clear(&(zNew->deviations[zTerm]));
        }
        else {
            mpfr_abs(temp, &(zNew->deviations[zTerm]), MPFR_RNDN);
            mpfr_add(&(zNew->radius), &(zNew->radius), temp, MPFR_RNDU);
            zTerm++;
        }
    }

    if (!mpfr_zero_p(error)) {
        zNew->symbols[zTerm] = mpfa_next_sym();
        mpfr_init2(&(zNew->deviations[zTerm]), prec_internal);
        mpfr_set(&(zNew->deviations[zTerm]), error, MPFR_RNDN);
        mpfr_add(&(zNew->radius), &(zNew->radius), error, MPFR_RNDU);
        zTerm++;
    }

    mpfr_prec_round(&(zNew->radius), prec, MPFR_RNDU);

    zNew->nTerms = zTerm;
    if (zNew->nTerms == 0) {
        free(zNew->symbols);
        free(zNew->deviations);
    }

    mpfr_clear(temp);
    mpfr_clear(error);
    mpfa_set(z, zNew);
    mpfa_clear(zNew);
    free(xTerm);
    free(summands);
}
