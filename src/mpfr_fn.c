/*
 * mpfr_fn.c -- Compute Arpra ranges with MPFR functions.
 *
 * Copyright 2019-2020 James Paul Turner.
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

#define ARPRA_MPFR_FN(SIGNATURE, MPFR_CALL)                             \
    void SIGNATURE                                                      \
    {                                                                   \
        mpfr_t error;                                                   \
        arpra_range yy;                                                 \
        arpra_prec prec_internal;                                       \
                                                                        \
        /* Initialise vars. */                                          \
        prec_internal = arpra_get_internal_precision();                 \
        mpfr_init2(error, prec_internal);                               \
        arpra_init2(&yy, y->precision);                                 \
        mpfr_set_zero(error, 1);                                        \
        mpfr_set_zero(&(yy.radius), 1);                                 \
                                                                        \
        /* y[0] = fn(x) */                                              \
        MPFR_CALL;                                                      \
                                                                        \
        /* Allocate memory for deviation terms. */                      \
        yy.symbols = malloc(sizeof(arpra_uint));                        \
        yy.deviations = malloc(sizeof(mpfr_t));                         \
                                                                        \
        /* Store new deviation term. */                                 \
        yy.symbols[0] = arpra_helper_next_symbol();                     \
        yy.deviations[0] = *error;                                      \
        mpfr_add(&(yy.radius), &(yy.radius), &(yy.deviations[0]), MPFR_RNDU); \
        yy.nTerms = 1;                                                  \
                                                                        \
        /* Compute true_range. */                                       \
        arpra_helper_compute_range(&yy);                                \
                                                                        \
        /* Check for NaN and Inf. */                                    \
        arpra_helper_check_result(&yy);                                 \
                                                                        \
        /* Clear vars, and set y. */                                    \
        arpra_clear(y);                                                 \
        *y = yy;                                                        \
    }

#define FN_SIGNATURE(FN_TYPE, ...)                                      \
    arpra_mpfr_##FN_TYPE (int (*fn) (mpfr_ptr y, __VA_ARGS__, mpfr_rnd_t rnd), arpra_range *y, __VA_ARGS__)

#define FN_MPFR_CALL(...)                                               \
    ARPRA_MPFR_RNDERR(error, MPFR_RNDN, fn, &(yy.centre), __VA_ARGS__)

// arpra_mpfr_fn1 (fn, arpra_range *y, mpfr_t x1)
ARPRA_MPFR_FN(FN_SIGNATURE(fn1, mpfr_srcptr x1), FN_MPFR_CALL(x1))

// arpra_mpfr_fn1_ui (fn, arpra_range *y, unsigned long int x1)
ARPRA_MPFR_FN(FN_SIGNATURE(fn1_ui, unsigned long int x1), FN_MPFR_CALL(x1))

// arpra_mpfr_fn1_si (fn, arpra_range *y, long int x1)
ARPRA_MPFR_FN(FN_SIGNATURE(fn1_si, long int x1), FN_MPFR_CALL(x1))

// arpra_mpfr_fn1_d (fn, arpra_range *y, double x1)
ARPRA_MPFR_FN(FN_SIGNATURE(fn1_d, double x1), FN_MPFR_CALL(x1))

// arpra_mpfr_fn2 (fn, arpra_range *y, mpfr_t x1, mpfr_t x2)
ARPRA_MPFR_FN(FN_SIGNATURE(fn2, mpfr_srcptr x1, mpfr_srcptr x2), FN_MPFR_CALL(x1, x2))

// arpra_mpfr_ui_fn2 (fn, arpra_range *y, unsigned long int x1, mpfr_t x2)
ARPRA_MPFR_FN(FN_SIGNATURE(ui_fn2, unsigned long int x1, mpfr_srcptr x2), FN_MPFR_CALL(x1, x2))

// arpra_mpfr_fn2_ui (fn, arpra_range *y, mpfr_t x1, unsigned long int x2)
ARPRA_MPFR_FN(FN_SIGNATURE(fn2_ui, mpfr_srcptr x1, unsigned long int x2), FN_MPFR_CALL(x1, x2))

// arpra_mpfr_si_fn2 (fn, arpra_range *y, long int x1, mpfr_t x2)
ARPRA_MPFR_FN(FN_SIGNATURE(si_fn2, long int x1, mpfr_srcptr x2), FN_MPFR_CALL(x1, x2))

// arpra_mpfr_fn2_si (fn, arpra_range *y, mpfr_t x1, long int x2)
ARPRA_MPFR_FN(FN_SIGNATURE(fn2_si, mpfr_srcptr x1, long int x2), FN_MPFR_CALL(x1, x2))

// arpra_mpfr_d_fn2 (fn, arpra_range *y, double x1, mpfr_t x2)
ARPRA_MPFR_FN(FN_SIGNATURE(d_fn2, double x1, mpfr_srcptr x2), FN_MPFR_CALL(x1, x2))

// arpra_mpfr_fn2_d (fn, arpra_range *y, mpfr_t x1, double x2)
ARPRA_MPFR_FN(FN_SIGNATURE(fn2_d, mpfr_srcptr x1, double x2), FN_MPFR_CALL(x1, x2))

// arpra_mpfr_set_str (arpra_range *y, char *x1, int base)
#define SET_STR_SIGNATURE arpra_mpfr_set_str (arpra_range *y, const char *x1, int base)
#define SET_STR_MPFR_CALL ARPRA_MPFR_RNDERR(error, MPFR_RNDN, mpfr_set_str, &(yy.centre), x1, base)
ARPRA_MPFR_FN(SET_STR_SIGNATURE, SET_STR_MPFR_CALL)
