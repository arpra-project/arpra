/*
 * helper_compute_range.c -- Compute the true_range field of an Arpra range.
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

/*
 * Compute true_range from centre and radius, adding rounding error to the
 * new deviation term.
 */

void arpra_helper_compute_range (arpra_range *y)
{
    mpfr_t temp1, temp2;
    arpra_prec prec_internal;
    arpra_uint iy;

    // Initialise vars.
    prec_internal = arpra_get_internal_precision();
    mpfr_init2(temp1, prec_internal + 8);
    mpfr_init2(temp2, prec_internal + 8);
    iy = y->nTerms - 1;

    // Compute true_range.
    mpfr_sub(temp1, &(y->centre), &(y->radius), MPFR_RNDD);
    mpfr_add(temp2, &(y->centre), &(y->radius), MPFR_RNDU);
    mpfr_set(&(y->true_range.left), temp1, MPFR_RNDD);
    mpfr_set(&(y->true_range.right), temp2, MPFR_RNDU);
    mpfr_sub(temp1, temp1, &(y->true_range.left), MPFR_RNDU);
    mpfr_sub(temp2, &(y->true_range.right), temp2, MPFR_RNDU);
    mpfr_max(temp1, temp1, temp2, MPFR_RNDU);
    mpfr_add(&(y->deviations[iy]), &(y->deviations[iy]), temp1, MPFR_RNDU);
    mpfr_add(&(y->radius), &(y->radius), temp1, MPFR_RNDU);

    // Clear vars.
    mpfr_clear(temp1);
    mpfr_clear(temp2);
}
