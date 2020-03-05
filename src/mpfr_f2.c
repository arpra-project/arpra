/*
 * mpfr_f2.c -- Compute Arpra ranges with bivariate MPFR functions.
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



// SET_MPFR.C: DO SET_MPFR AND USE THIS FOR UINT, INT, FLOAT, ETC IN SET_OTHER.C

// WIDEN_MPFR.C, AND USE THIS ALSO FOR UINT INT FLOAT IN WIDEN_OTHER.C


// REPLACE ALL _MPFR WITH _FR, LIKE MPFI DOES

// SET.C SHOULD BE SET_RANGE.C, SET SHOULD BE SET_RANGE
// OTHER SCALAR TYPES SHOULD BE IN SET_NUMBER


// ARPRA_MPFR_F* SHOULD BE HELPER FUNCTIONS - EXPOSE REAL FUNCTIONS INSTEAD


void arpra_mpfr_f2 (int (*f) (mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mpfr_rnd_t),
                    arpra_range *y, mpfr_srcptr x1, mpfr_srcptr x2)
{
    mpfr_t error;
    arpra_range yy;
    arpra_prec prec_internal;

    // Initialise vars.
    prec_internal = arpra_get_internal_precision();
    arpra_init2(&yy, y->precision);
    mpfr_init2(error, prec_internal);
    mpfr_set_zero(error, 1);
    mpfr_set_zero(&(yy.radius), 1);

    // y[0] = f(x1, x2)
    arpra_helper_mpfr_f2(error, f, &(yy.centre), x1, x2, MPFR_RNDN);

    // Allocate memory for deviation terms.
    yy.symbols = malloc(sizeof(arpra_uint));
    yy.deviations = malloc(sizeof(mpfr_t));

    // Store numerical error term.
    yy.symbols[0] = arpra_helper_next_symbol();
    yy.deviations[0] = *error;
    mpfr_add(&(yy.radius), &(yy.radius), &(yy.deviations[0]), MPFR_RNDU);
    yy.nTerms = 1;

    // Compute true_range, and add rounding error.
    arpra_helper_range_rounded(&yy);

    // Check for NaN and Inf.
    arpra_helper_check_result(&yy);

    // Clear vars, and set y.
    arpra_clear(y);
    *y = yy;
}
