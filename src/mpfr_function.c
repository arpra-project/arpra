/*
 * mpfr_function.c -- Compute Arpra range from MPFR functions.
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

void arpra_mpfr_function_1 (arpra_range *z, const arpra_range *x,
                            int (*f) (mpfr_ptr, mpfr_srcptr, mpfr_rnd_t))
{
    arpra_prec prec_internal;

    // Initialise vars.
    prec_internal = arpra_get_internal_precision();
    mpfr_set_prec(&(z->centre), prec_internal);
    mpfr_set_prec(&(z->radius), prec_internal);

    // Clear existing deviation terms.
    arpra_clear_terms(z);


    // ============ TODO: RESET RADIUS AND TRUE_RANGE

    // z_0 = f(x_0)
    if (f(&(z->centre), &(x->centre), MPFR_RNDN)) {
        z->nTerms = 1;
        z->symbols = malloc(sizeof(arpra_uint));
        z->deviations = malloc(sizeof(arpra_mpfr));
        z->symbols[0] = arpra_helper_next_symbol();
        mpfr_init2(&(z->deviations[0]), prec_internal);
        arpra_helper_error_half_ulp(&(z->deviations[0]), &(z->centre));
        mpfr_set(&(z->radius), &(z->deviations[0]), MPFR_RNDU);

        // SET TRUE_RANGE

    }
}

void arpra_mpfr_function_2 (arpra_range *z, const arpra_range *x, const arpra_range *y)
{

}
