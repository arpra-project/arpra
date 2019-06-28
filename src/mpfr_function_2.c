/*
 * mpfr_function_2.c -- Compute Arpra range from bivariate MPFR functions.
 *
 * Copyright 2019 James Paul Turner.
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

void arpra_mpfr_function_2 (arpra_range *z, const arpra_mpfr *x, const arpra_mpfr *y,
                            int (*f) (mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mpfr_rnd_t))
{
    arpra_mpfr error;
    arpra_range zNew;
    arpra_prec prec_internal;

    // Initialise vars.
    prec_internal = arpra_get_internal_precision();
    mpfr_init2(&error, prec_internal);
    arpra_init2(&zNew, z->precision);
    mpfr_set_zero(&error, 1);
    mpfr_set_zero(&(zNew.radius), 1);

    // z_0 = f(x, y)
    if (f(&(zNew.centre), &(x->centre), &(y->centre), MPFR_RNDN)) {
        arpra_helper_error_half_ulp(&error, &(zNew.centre));
    }

    // Allocate memory for deviation term.
    zNew.symbols = malloc(sizeof(arpra_uint));
    zNew.deviations = malloc(sizeof(arpra_mpfr));

    // Store numerical error term.
    zNew.symbols[0] = arpra_helper_next_symbol();
    zNew.deviations[0] = error;
    mpfr_add(&(zNew.radius), &(zNew.radius), &(zNew.deviations[0]), MPFR_RNDU);
    zNew.nTerms = 1;

    // Compute true_range, and add rounding error.
    arpra_helper_range_rounded(&zNew);

    // Check for NaN and Inf.
    arpra_helper_check_result(&zNew);

    // Clear vars, and set z.
    arpra_clear(z);
    *z = zNew;
}
