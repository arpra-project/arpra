/*
 * affine_2.c -- Compute a bivariate affine function of two affine forms.
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

void mpfa_affine_2 (mpfa_ptr z, mpfa_srcptr x, mpfa_srcptr y, mpfr_srcptr alpha, mpfr_srcptr beta, mpfr_srcptr gamma, mpfr_srcptr delta) {
    unsigned xTerm, yTerm, zTerm;
    int xHasNext, yHasNext;
    mpfr_t temp, error;
    mpfr_prec_t prec;
    mpfa_t zNew;

    prec = mpfr_get_prec(&(z->centre));
    mpfr_inits2(prec, temp, error, (mpfr_ptr) NULL);
    mpfa_init2(zNew, prec);
    mpfr_set(error, delta, MPFR_RNDU);
    mpfr_set_si(&(zNew->radius), 0, MPFR_RNDN);

    if (mpfa_term(&(zNew->centre), &(x->centre), &(y->centre), alpha, beta, gamma)) {
        mpfr_mul(temp, &(zNew->centre), &(zNew->u), MPFR_RNDU);
        mpfr_add(error, error, temp, MPFR_RNDU);
    }

    zNew->nTerms = x->nTerms + y->nTerms + 1;
    zNew->symbols = malloc(zNew->nTerms * sizeof(unsigned));
    zNew->deviations = malloc(zNew->nTerms * sizeof(mpfr_t));

    xTerm = 0;
    yTerm = 0;
    zTerm = 0;
    xHasNext = x->nTerms > 0;
    yHasNext = y->nTerms > 0;
    while (xHasNext || yHasNext) {
        if ((!yHasNext) || (xHasNext && (x->symbols[xTerm] < y->symbols[yTerm]))) {
            zNew->symbols[zTerm] = x->symbols[xTerm];
            mpfr_init2(&(zNew->deviations[zTerm]), prec);

            if (mpfr_mul(&(zNew->deviations[zTerm]), alpha, &(x->deviations[xTerm]), MPFR_RNDN)) {
                mpfr_mul(temp, &(zNew->deviations[zTerm]), &(zNew->u), MPFR_RNDU);
                mpfr_add(error, error, temp, MPFR_RNDU);
            }

            xHasNext = ++xTerm < x->nTerms;
        }
        else if ((!xHasNext) || (yHasNext && (y->symbols[yTerm] < x->symbols[xTerm]))) {
            zNew->symbols[zTerm] = y->symbols[yTerm];
            mpfr_init2(&(zNew->deviations[zTerm]), prec);

            if (mpfr_mul(&(zNew->deviations[zTerm]), beta, &(y->deviations[yTerm]), MPFR_RNDN)) {
                mpfr_mul(temp, &(zNew->deviations[zTerm]), &(zNew->u), MPFR_RNDU);
                mpfr_add(error, error, temp, MPFR_RNDU);
            }

            yHasNext = ++yTerm < y->nTerms;
        }
        else {
            zNew->symbols[zTerm] = x->symbols[xTerm];
            mpfr_init2(&(zNew->deviations[zTerm]), prec);

            if (mpfa_term(&(zNew->deviations[zTerm]), &(x->deviations[xTerm]), &(y->deviations[yTerm]), alpha, beta, NULL)) {
                mpfr_mul(temp, &(zNew->deviations[zTerm]), &(zNew->u), MPFR_RNDU);
                mpfr_add(error, error, temp, MPFR_RNDU);
            }

            xHasNext = ++xTerm < x->nTerms;
            yHasNext = ++yTerm < y->nTerms;
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
        mpfr_init2(&(zNew->deviations[zTerm]), prec);
        mpfr_set(&(zNew->deviations[zTerm]), error, MPFR_RNDN);
        mpfr_add(&(zNew->radius), &(zNew->radius), error, MPFR_RNDU);
        zTerm++;
    }

    zNew->nTerms = zTerm;
    if (zNew->nTerms == 0) {
        free(zNew->symbols);
        free(zNew->deviations);
    }

    mpfr_clears(temp, error, (mpfr_ptr) NULL);
    mpfa_set(z, zNew);
    mpfa_clear(zNew);
}
