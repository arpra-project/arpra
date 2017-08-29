/*
 * init.c -- Initialise one or more affine forms.
 *
 * Copyright 2016-2017 James Paul Turner.
 *
 * This file is part of the MPFA library.
 *
 * The MPFA library is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * The MPFA library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with the MPFA library. If not, see <http://www.gnu.org/licenses/>.
 */

#include "mpfa-impl.h"

void mpfa_init (mpfa_ptr x) {
    // Init centre and radius with working precision.
    x->nTerms = 0;
    mpfr_init(&(x->centre));
    mpfr_init(&(x->radius));
}

void mpfa_inits (mpfa_ptr x, ...) {
    va_list arg;

    // Init each argument with working precision.
    va_start(arg, x);
    while (x != NULL) {
        mpfa_init(x);
        x = (mpfa_ptr) va_arg(arg, mpfa_ptr);
    }
    va_end(arg);
}
