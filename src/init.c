/*
 * init.c -- Initialise the memory of Arpra ranges.
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

void arpra_init (arpra_range *y)
{
    arpra_prec prec;

    prec = arpra_get_default_precision();
    arpra_init2(y, prec);
}

void arpra_init2 (arpra_range *y, arpra_prec prec)
{
    arpra_prec prec_internal;

    prec_internal = arpra_get_internal_precision();
    mpfr_init2(&(y->centre), prec_internal);
    mpfr_init2(&(y->radius), prec_internal);
    mpfi_init2(&(y->true_range), prec);
    y->nTerms = 0;
    y->precision = prec;
}
