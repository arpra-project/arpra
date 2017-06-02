/*
 * inits2.c -- Initialise a list of affine forms, and set their precision.
 *
 * Copyright 2017 James Paul Turner.
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

#include "mpfa.h"
#include <stdarg.h>

void mpfa_inits2 (mpfr_prec_t prec, mpfa_ptr x, ...) {
    va_list arg;

    va_start(arg, x);
    while (x != NULL) {
	mpfa_init2(x, prec);
	x = (mpfa_ptr) va_arg(arg, mpfa_ptr);
    }
    va_end(arg);
}
