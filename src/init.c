/*
 * init.c -- Initialise one or more Arpra ranges.
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

void arpra_init (arpra_range *z)
{
    arpra_precision prec, prec_internal;

    // Init centre and radius.
    prec = arpra_get_default_precision();
    prec_internal = arpra_get_internal_precision();
    z->nTerms = 0;
    mpfr_init2(&(z->centre), prec);
    mpfr_init2(&(z->radius), prec_internal);
}

void arpra_inits (arpra_range *z, ...)
{
    va_list arg;

    // Init each argument's centre and radius.
    va_start(arg, z);
    while (z != NULL) {
        arpra_init(z);
        z = (arpra_range *) va_arg(arg, arpra_range *);
    }
    va_end(arg);
}
