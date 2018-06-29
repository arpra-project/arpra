/*
 * precision.c -- Get and set the precision of an Arpra range.
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

arpra_precision arpra_get_precision (const arpra_range *x)
{
    return mpfr_get_prec(&(x->centre));
}

void arpra_set_precision (arpra_range *z, const arpra_precision prec)
{
    arpra_precision internal_prec;

    // Increase internal_prec if < prec.
    internal_prec = arpra_get_internal_precision();
    if (internal_prec < prec) {
        arpra_set_internal_precision(prec);
        internal_prec = prec;
    }

    // Clear existing deviation terms.
    arpra_clear_terms(z);

    // Reset centre and radius precision.
    mpfr_set_prec(&(z->centre), prec);
    mpfr_set_prec(&(z->radius), internal_prec);
}
