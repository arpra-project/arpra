/*
 * add.c -- Add one Arpra range to another.
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

void arpra_add (arpra_range *y, const arpra_range *x1, const arpra_range *x2)
{
    mpfi_t ia_range;
    mpfr_t alpha, beta, gamma, delta;

    // Domain violations:
    // (NaN) + (NaN) = (NaN)
    // (NaN) + (R)   = (NaN)
    // (R)   + (NaN) = (NaN)
    // (Inf) + (Inf) = (NaN)
    // (Inf) + (R)   = (Inf)
    // (R)   + (Inf) = (Inf)

    // Handle domain violations.
    if (arpra_nan_p(x1) || arpra_nan_p(x2)) {
        arpra_set_nan(y);
        return;
    }
    if (arpra_inf_p(x1)) {
        if (arpra_inf_p(x2)) {
            arpra_set_nan(y);
        }
        else {
            arpra_set_inf(y);
        }
        return;
    }
    if (arpra_inf_p(x2)) {
        if (arpra_inf_p(x1)) {
            arpra_set_nan(y);
        }
        else {
            arpra_set_inf(y);
        }
        return;
    }

    // Initialise vars.
    mpfi_init2(ia_range, y->precision);
    mpfr_init2(alpha, 2);
    mpfr_init2(beta, 2);
    mpfr_init2(gamma, 2);
    mpfr_init2(delta, 2);
    mpfr_set_si(alpha, 1, MPFR_RNDN);
    mpfr_set_si(beta, 1, MPFR_RNDN);
    mpfr_set_zero(gamma, 1);
    mpfr_set_zero(delta, 1);

    // MPFI addition
    mpfi_add(ia_range, &(x1->true_range), &(x2->true_range));

    // y = x1 + x2
    arpra_helper_affine_2(y, x1, x2, alpha, beta, gamma, delta);

    // Compute true_range.
    arpra_helper_compute_range(y);

    // Mix with IA range, and trim error term.
    arpra_helper_mix_trim(y, ia_range);

    // Check for NaN and Inf.
    arpra_helper_check_result(y);

    // Clear vars.
    mpfi_clear(ia_range);
    mpfr_clear(alpha);
    mpfr_clear(beta);
    mpfr_clear(gamma);
    mpfr_clear(delta);
}
