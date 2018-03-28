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
    mpfr_t temp, error;
    arpra_precision prec, prec_internal;
    arpra_range zNew;

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
    prec = arpra_get_precision(z);
    prec_internal = arpra_get_internal_precision();
    mpfr_init2(temp, prec_internal);
    mpfr_init2(error, prec_internal);
    mpfr_init2(&(zNew.centre), prec);
    mpfr_init2(&(zNew.radius), prec_internal);
    mpfr_set_si(error, 0, MPFR_RNDU);
    mpfr_set_si(&(zNew.radius), 0, MPFR_RNDU);

    // z_0 = x_0 * y_0
    if (mpfr_mul(&(zNew.centre), &(x->centre), &(y->centre), MPFR_RNDN)) {
        arpra_error(temp, &(zNew.centre));
        mpfr_add(error, error, temp, MPFR_RNDU);
    }

    // Allocate memory for all possible deviation terms.
    zNew.nTerms = x->nTerms + y->nTerms + 1;
    zNew.symbols = malloc(zNew.nTerms * sizeof(arpra_uint));
    zNew.deviations = malloc(zNew.nTerms * sizeof(mpfr_t));

    xTerm = 0;
    yTerm = 0;
    zTerm = 0;
    xHasNext = x->nTerms > 0;
    yHasNext = y->nTerms > 0;
    while (xHasNext || yHasNext) {
        if ((!yHasNext) || (xHasNext && (x->symbols[xTerm] < y->symbols[yTerm]))) {
            zNew.symbols[zTerm] = x->symbols[xTerm];
            mpfr_init2(&(zNew.deviations[zTerm]), prec);

            // z_i = (y_0 * x_i)
            if (mpfr_mul(&(zNew.deviations[zTerm]), &(y->centre), &(x->deviations[xTerm]), MPFR_RNDN)) {
                arpra_error(temp, &(zNew.deviations[zTerm]));
                mpfr_add(error, error, temp, MPFR_RNDU);
            }

            xHasNext = ++xTerm < x->nTerms;
        }
        else if ((!xHasNext) || (yHasNext && (y->symbols[yTerm] < x->symbols[xTerm]))) {
            zNew.symbols[zTerm] = y->symbols[yTerm];
            mpfr_init2(&(zNew.deviations[zTerm]), prec);

            // z_i = (x_0 * y_i)
            if (mpfr_mul(&(zNew.deviations[zTerm]), &(x->centre), &(y->deviations[yTerm]), MPFR_RNDN)) {
                arpra_error(temp, &(zNew.deviations[zTerm]));
                mpfr_add(error, error, temp, MPFR_RNDU);
            }

            yHasNext = ++yTerm < y->nTerms;
        }
        else {
            zNew.symbols[zTerm] = x->symbols[xTerm];
            mpfr_init2(&(zNew.deviations[zTerm]), prec);

            // z_i = (y_0 * x_i) + (x_0 * y_i)
            if (arpra_term(&(zNew.deviations[zTerm]), &(x->deviations[xTerm]), &(y->deviations[yTerm]), &(y->centre), &(x->centre), NULL)) {
                arpra_error(temp, &(zNew.deviations[zTerm]));
                mpfr_add(error, error, temp, MPFR_RNDU);
            }

            xHasNext = ++xTerm < x->nTerms;
            yHasNext = ++yTerm < y->nTerms;
        }

        // Store nonzero deviation terms.
        if (mpfr_zero_p(&(zNew.deviations[zTerm]))) {
            mpfr_clear(&(zNew.deviations[zTerm]));
        }
        else {
            mpfr_abs(temp, &(zNew.deviations[zTerm]), MPFR_RNDU);
            mpfr_add(&(zNew.radius), &(zNew.radius), temp, MPFR_RNDU);
            zTerm++;
        }
    }

#ifdef ARPRA_TIGHT_MUL
    // Linear approximation of the quadratic term is defined the same as in (26) of:
    // S. M. Rump and M. Kashiwagi, Implementation and improvements of affine arithmetic,
    // Nonlinear Theory an Its Applications, IEICE, vol. 6, no. 3, pp. 341-359, 2015.

    arpra_uint xNext, yNext;
    mpfr_t xiyiPos, xiyiNeg;

    // Init extra temp vars.
    mpfr_init2(xiyiPos, prec_internal);
    mpfr_init2(xiyiNeg, prec_internal);
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

    // Clear extra temp vars.
    mpfr_clear(xiyiPos);
    mpfr_clear(xiyiNeg);
#else
    // Linear approximation of quadratic term is rad(x) * rad(y).
    mpfr_mul(temp, &(x->radius), &(y->radius), MPFR_RNDU);
    mpfr_add(error, error, temp, MPFR_RNDU);
#endif

    // Store nonzero numerical error term.
    if (!mpfr_zero_p(error)) {
        zNew.symbols[zTerm] = arpra_next_symbol();
        mpfr_init2(&(zNew.deviations[zTerm]), prec);
        mpfr_set(&(zNew.deviations[zTerm]), error, MPFR_RNDU);
        mpfr_add(&(zNew.radius), &(zNew.radius), &(zNew.deviations[zTerm]), MPFR_RNDU);
        zTerm++;
    }

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
    mpfr_clear(temp);
    mpfr_clear(error);
    arpra_set(z, &zNew);
    arpra_clear(&zNew);
}
