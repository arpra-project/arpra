/*
 * prec.c -- Get and set the precision of an arpra_t.
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

arpra_prec_t arpra_get_prec (arpra_srcptr x)
{
    return mpfr_get_prec(&(x->centre));
}

void arpra_set_prec (arpra_ptr z, arpra_prec_t prec)
{
    arpra_prec_t prec_internal;

    // Increase internal precision if < 'prec'.
    prec_internal = arpra_get_internal_prec();
    if (prec_internal < prec) {
        arpra_set_internal_prec(prec);
        prec_internal = prec;
    }

    // Clear existing deviation terms.
    arpra_clear_terms(z);

    // Reset centre and radius with new working precision.
    mpfr_set_prec(&(z->centre), prec);
    mpfr_set_prec(&(z->radius), prec_internal);
}
