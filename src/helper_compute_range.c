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

void arpra_helper_compute_range (arpra_range *z)
{
    arpra_uint zTerm;
    arpra_mpfr temp1, temp2;
    arpra_prec prec_internal;

    // Initialise vars.
    prec_internal = arpra_get_internal_precision();
    mpfr_init2(&temp1, prec_internal + 8);
    mpfr_init2(&temp2, prec_internal + 8);
    zTerm = z->nTerms - 1;

    // Compute true_range.
    mpfr_sub(&temp1, &(z->centre), &(z->radius), MPFR_RNDD);
    mpfr_add(&temp2, &(z->centre), &(z->radius), MPFR_RNDU);
    mpfr_set(&(z->true_range.left), &temp1, MPFR_RNDD);
    mpfr_set(&(z->true_range.right), &temp2, MPFR_RNDU);
    mpfr_sub(&temp1, &temp1, &(z->true_range.left), MPFR_RNDU);
    mpfr_sub(&temp2, &(z->true_range.right), &temp2, MPFR_RNDU);
    mpfr_max(&temp1, &temp1, &temp2, MPFR_RNDU);
    mpfr_add(&(z->deviations[zTerm]), &(z->deviations[zTerm]), &temp1, MPFR_RNDU);
    mpfr_add(&(z->radius), &(z->radius), &temp1, MPFR_RNDU);

    // Clear vars.
    mpfr_clear(&temp1);
    mpfr_clear(&temp2);
}
