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

arpra_prec arpra_get_precision (const arpra_range *x)
{
    return x->precision;
}

void arpra_set_precision (arpra_range *z, const arpra_prec prec)
{
    arpra_prec prec_internal;

    prec_internal = arpra_get_internal_precision();
    mpfr_set_prec(&(z->centre), prec_internal);
    mpfr_set_prec(&(z->radius), prec_internal);
    mpfi_set_prec(&(z->true_range), prec);
    arpra_clear_terms(z);
    z->precision = prec;
}
