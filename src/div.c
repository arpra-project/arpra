/*
 * div.c -- Divide one Arpra range by another.
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

void arpra_div (arpra_range *y, const arpra_range *x1, const arpra_range *x2)
{
    mpfi_t ia_range;
    mpfr_t temp1, temp2;
    arpra_range yy;
    arpra_prec prec_internal;

    // Domain violations:
    // (NaN) / (NaN) = (NaN)
    // (NaN) / (R)   = (NaN)
    // (R)   / (NaN) = (NaN)
    // (Inf) / (Inf) = (NaN)
    // (Inf) / (0)   = (NaN)
    // (0)   / (Inf) = (NaN)
    // (0)   / (0)   = (NaN)
    // (Inf) / (R)   = (Inf)
    // (R)   / (Inf) = (Inf)

    // Handle domain violations.
    if (arpra_nan_p(x1) || arpra_nan_p(x2)) {
        arpra_set_nan(y);
        return;
    }
    if (arpra_has_zero_p(x1) && arpra_has_zero_p(x2)) {
        arpra_set_nan(y);
        return;
    }
    if (arpra_inf_p(x1) || arpra_inf_p(x2)) {
        arpra_set_inf(y);
        return;
    }

    // Initialise vars.
    prec_internal = arpra_get_internal_precision();
    mpfi_init2(ia_range, y->precision);
    mpfr_init2(temp1, prec_internal);
    mpfr_init2(temp2, prec_internal);
    arpra_init2(&yy, y->precision);

    // MPFI division
    mpfi_div(ia_range, &(x1->true_range), &(x2->true_range));

    // y = x1 * (1 / x2)
    arpra_inv(&yy, x2);
    arpra_mul(y, x1, &yy);

    // Mix with IA range, and trim error term.
    arpra_helper_mix_trim(y, ia_range);

    // Check for NaN and Inf.
    arpra_helper_check_result(y);

    // Clear vars.
    mpfi_clear(ia_range);
    mpfr_clear(temp1);
    mpfr_clear(temp2);
    arpra_clear(&yy);
}
