/*
 * increase.c -- Append a new numerical error term to an Arpra range.
 *
 * Copyright 2020 James Paul Turner.
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

void arpra_increase (arpra_range *y, const arpra_range *x1, mpfr_srcptr delta)
{
    mpfi_t ia_range;
    mpfr_t alpha, gamma;

    // Domain violations:
    // increase(NaN) = (NaN)
    // increase(Inf) = (Inf)

    // Handle domain violations.
    if (arpra_nan_p(x1)) {
        arpra_set_nan(y);
        return;
    }
    if (arpra_inf_p(x1)) {
        arpra_set_inf(y);
        return;
    }

    // Initialise vars.
    mpfi_init2(ia_range, y->precision);
    mpfr_init2(alpha, 2);
    mpfr_init2(gamma, 2);
    mpfr_set_si(alpha, 1, MPFR_RNDN);
    mpfr_set_zero(gamma, 1);

    // MPFI increase
    mpfi_set(ia_range, &(x1->true_range));
    mpfi_increase(ia_range, delta);

    // y = increase(x1, delta)
    arpra_helper_affine_1(y, x1, alpha, gamma, delta);

    // Compute true_range.
    arpra_helper_compute_range(y);

    // Mix with IA range, and trim error term.
    arpra_helper_mix_trim(y, ia_range);

    // Check for NaN and Inf.
    arpra_helper_check_result(y);

    // Clear vars.
    mpfi_clear(ia_range);
    mpfr_clear(alpha);
    mpfr_clear(gamma);
}
