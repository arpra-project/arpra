/*
 * helper_sumabs.c -- Compute the absolute value sum of n MPFR numbers.
 *
 * Copyright 2018 James Paul Turner.
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

static arpra_uint buffer_size = 0;
static arpra_sign *sign_buffer = NULL;
static arpra_mpfr **term_buffer = NULL;

void arpra_helper_sumabs (arpra_range *sumabs, const arpra_mpfr *x, const arpra_uint n)
{
    // Allocate or double temp buffers, as required.
    if (buffer_size < n) {
        buffer_size = buffer_size ? buffer_size * 2 : ARPRA_DEFAULT_BUFFER_SIZE;
        sign_buffer = realloc(sign_buffer, buffer_size * sizeof(arpra_sign));
        term_buffer = realloc(term_buffer, buffer_size * sizeof(arpra_mpfr *));
    }

}
