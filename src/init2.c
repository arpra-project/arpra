/*
 * init2.c -- Initialise one or more arpra_t, and set their precision.
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

void arpra_init2 (arpra_ptr x, arpra_prec_t prec)
{
    arpra_prec_t prec_internal;

    // Increase internal precision if < 'prec'.
    prec_internal = arpra_get_internal_prec();
    if (prec_internal < prec) {
        arpra_set_internal_prec(prec);
        prec_internal = prec;
    }

    // Init centre and radius with 'prec' precision.
    x->nTerms = 0;
    mpfr_init2(&(x->centre), prec);
    mpfr_init2(&(x->radius), prec_internal);
}

void arpra_inits2 (arpra_prec_t prec, arpra_ptr x, ...)
{
    va_list arg;

    // Init each argument with precision 'prec'.
    va_start(arg, x);
    while (x != NULL) {
        arpra_init2(x, prec);
        x = (arpra_ptr) va_arg(arg, arpra_ptr);
    }
    va_end(arg);
}
