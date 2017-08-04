/*
 * mul.c -- Multiply one affine form by another.
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

void mpfa_mul (mpfa_ptr z, mpfa_srcptr x, mpfa_srcptr y) {
    mpfa_uint_t xTerm, yTerm, zTerm;
    mpfa_int_t xHasNext, yHasNext;
    mpfr_t temp, error;
    mpfa_prec_t prec, prec_internal;
    mpfa_t zNew;

    prec = mpfr_get_prec(&(z->centre));
    prec_internal = mpfa_get_internal_prec();
    mpfr_init2(temp, prec);
    mpfr_init2(error, prec);
    mpfr_init2(&(zNew->centre), prec);
    mpfr_init2(&(zNew->radius), prec_internal);
    mpfr_set_si(error, 0, MPFR_RNDN);
    mpfr_set_si(&(zNew->radius), 0, MPFR_RNDN);

    if (mpfr_mul(&(zNew->centre), &(x->centre), &(y->centre), MPFR_RNDN)) {
        mpfa_error(temp, &(zNew->centre));
        mpfr_add(error, error, temp, MPFR_RNDU);
    }

    zNew->nTerms = x->nTerms + y->nTerms + 1;
    zNew->symbols = malloc(zNew->nTerms * sizeof(mpfa_uint_t));
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

            if (mpfr_mul(&(zNew->deviations[zTerm]), &(y->centre), &(x->deviations[xTerm]), MPFR_RNDN)) {
                mpfa_error(temp, &(zNew->deviations[zTerm]));
                mpfr_add(error, error, temp, MPFR_RNDU);
            }

            xHasNext = ++xTerm < x->nTerms;
        }
        else if ((!xHasNext) || (yHasNext && (y->symbols[yTerm] < x->symbols[xTerm]))) {
            zNew->symbols[zTerm] = y->symbols[yTerm];
            mpfr_init2(&(zNew->deviations[zTerm]), prec);

            if (mpfr_mul(&(zNew->deviations[zTerm]), &(x->centre), &(y->deviations[yTerm]), MPFR_RNDN)) {
                mpfa_error(temp, &(zNew->deviations[zTerm]));
                mpfr_add(error, error, temp, MPFR_RNDU);
            }

            yHasNext = ++yTerm < y->nTerms;
        }
        else {
            zNew->symbols[zTerm] = x->symbols[xTerm];
            mpfr_init2(&(zNew->deviations[zTerm]), prec);

            if (mpfa_term(&(zNew->deviations[zTerm]), &(x->deviations[xTerm]), &(y->deviations[yTerm]), &(y->centre), &(x->centre), NULL)) {
                mpfa_error(temp, &(zNew->deviations[zTerm]));
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

#ifdef MPFA_TIGHT_MUL
    mpfa_uint_t xNext, yNext;
    mpfr_t xiyiPos, xiyiNeg;

    mpfr_init2(xiyiPos, prec);
    mpfr_init2(xiyiNeg, prec);
    mpfr_set_si(xiyiPos, 0, MPFR_RNDN);
    mpfr_set_si(xiyiNeg, 0, MPFR_RNDN);

    xTerm = 0;
    yTerm = 0;
    while ((xTerm < x->nTerms) && (yTerm < y->nTerms)) {
        if (x->symbols[xTerm] < y->symbols[yTerm]) {
            for (yNext = yTerm; yNext < y->nTerms; yNext++) {
                // x has symbol i, and y has symbol j, so error += abs(xi * yj)
                mpfr_mul(temp, &(x->deviations[xTerm]), &(y->deviations[yNext]), MPFR_RNDA);
                mpfr_abs(temp, temp, MPFR_RNDN);
                mpfr_add(error, error, temp, MPFR_RNDU);
            }
            xTerm++;
        }
        else if (y->symbols[yTerm] < x->symbols[xTerm]) {
            for (xNext = xTerm; xNext < x->nTerms; xNext++) {
                // y has symbol i, and x has symbol j, so error += abs(yi * xj)
                mpfr_mul(temp, &(y->deviations[yTerm]), &(x->deviations[xNext]), MPFR_RNDA);
                mpfr_abs(temp, temp, MPFR_RNDN);
                mpfr_add(error, error, temp, MPFR_RNDU);
            }
            yTerm++;
        }
        else {
            // both x and y have symbol i, so error += abs(xi * yi)
            mpfr_mul(temp, &(x->deviations[xTerm]), &(y->deviations[yTerm]), MPFR_RNDA);
            if (mpfr_sgn(temp) >= 0) {
                mpfr_add(xiyiPos, xiyiPos, temp, MPFR_RNDU);
            }
            else {
                mpfr_sub(xiyiNeg, xiyiNeg, temp, MPFR_RNDU);
            }

            xNext = xTerm + 1;
            yNext = yTerm + 1;
            xHasNext = xNext < x->nTerms;
            yHasNext = yNext < y->nTerms;
            while (xHasNext || yHasNext) {
                if ((!yHasNext) || (xHasNext && (x->symbols[xTerm] < y->symbols[yTerm]))) {
                    // both x and y have symbol i, but only x has symbol j, so error += abs(xj * yi)
                    mpfr_mul(temp, &(y->deviations[yTerm]), &(x->deviations[xNext]), MPFR_RNDA);
                    mpfr_abs(temp, temp, MPFR_RNDN);
                    mpfr_add(error, error, temp, MPFR_RNDU);

                    xHasNext = ++xNext < x->nTerms;
                }
                else if ((!xHasNext) || (yHasNext && (y->symbols[yTerm] < x->symbols[xTerm]))) {
                    // both x and y have symbol i, but only y has symbol j, so error += abs(xi * yj)
                    mpfr_mul(temp, &(x->deviations[xTerm]), &(y->deviations[yNext]), MPFR_RNDA);
                    mpfr_abs(temp, temp, MPFR_RNDN);
                    mpfr_add(error, error, temp, MPFR_RNDU);

                    yHasNext = ++yNext < y->nTerms;
                }
                else {
                    // both x and y have symbols i and j, so error += abs(xi * yj + xj * yi)
                    mpfr_mul(temp, &(x->deviations[xTerm]), &(y->deviations[yNext]), MPFR_RNDU);
                    mpfr_mul(temp, &(x->deviations[xTerm]), &(y->deviations[yNext]), MPFR_RNDU);
                    mpfr_add(temp, temp, temp, MPFR_RNDU);
                    if (mpfr_sgn(temp) < 0) {
                        mpfr_mul(temp, &(x->deviations[xTerm]), &(y->deviations[yNext]), MPFR_RNDD);
                        mpfr_mul(temp, &(x->deviations[xTerm]), &(y->deviations[yNext]), MPFR_RNDD);
                        mpfr_add(temp, temp, temp, MPFR_RNDD);
                    }
                    mpfr_abs(temp, temp, MPFR_RNDN);
                    mpfr_add(error, error, temp, MPFR_RNDU);

                    xHasNext = ++xNext < x->nTerms;
                    yHasNext = ++yNext < y->nTerms;
                }
            }
            xTerm++;
            yTerm++;
        }
    }

    mpfr_max(temp, xiyiPos, xiyiNeg, MPFR_RNDN);
    mpfr_add(error, error, temp, MPFR_RNDU);

    mpfr_clear(xiyiPos);
    mpfr_clear(xiyiNeg);
#else
    mpfr_mul(temp, &(x->radius), &(y->radius), MPFR_RNDU);
    mpfr_add(error, error, temp, MPFR_RNDU);
#endif

    if (!mpfr_zero_p(error)) {
        zNew->symbols[zTerm] = mpfa_next_sym();
        mpfr_init2(&(zNew->deviations[zTerm]), prec);
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
}
