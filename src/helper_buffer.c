/*
 * helper_buffer.c -- Access and manage global buffers.
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

// MPFR pointer buffer.
static arpra_mpfr **buffer_mpfr_ptr = NULL;
static arpra_uint buffer_mpfr_ptr_size = 0;

arpra_mpfr **arpra_helper_buffer_mpfr_ptr (arpra_uint n)
{
    // Allocate or resize buffer, as required.
    if (buffer_mpfr_ptr_size < n) {
        buffer_mpfr_ptr_size = ceil(n / ARPRA_BUFFER_RESIZE_FACTOR);
        buffer_mpfr_ptr_size *= ARPRA_BUFFER_RESIZE_FACTOR;
        buffer_mpfr_ptr = realloc(buffer_mpfr_ptr, buffer_mpfr_ptr_size * sizeof(arpra_mpfr *));
    }

    return buffer_mpfr_ptr;
}

// MPFR buffer.
static arpra_mpfr *buffer_mpfr = NULL;
static arpra_uint buffer_mpfr_size = 0;

arpra_mpfr *arpra_helper_buffer_mpfr (arpra_uint n)
{
    // Allocate or resize buffer, as required.
    if (buffer_mpfr_size < n) {
        buffer_mpfr_size = ceil(n / ARPRA_BUFFER_RESIZE_FACTOR);
        buffer_mpfr_size *= ARPRA_BUFFER_RESIZE_FACTOR;
        buffer_mpfr = realloc(buffer_mpfr, buffer_mpfr_size * sizeof(arpra_mpfr));
    }

    return buffer_mpfr;
}

void arpra_clear_buffers ()
{
    // Free MPFR pointer buffer.
    free(buffer_mpfr_ptr);
    buffer_mpfr_ptr = NULL;
    buffer_mpfr_ptr_size = 0;

    // Free MPFR buffer.
    free(buffer_mpfr);
    buffer_mpfr = NULL;
    buffer_mpfr_size = 0;
}
