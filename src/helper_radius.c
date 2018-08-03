/*
 * helper_radius.c -- Compute radiuses and sums of abs value MPFR numbers.
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

static arpra_uint pointer_buffer_size = 0;
static arpra_mpfr **pointer_buffer = NULL;
static arpra_uint sign_buffer_size = 0;
static arpra_sign *sign_buffer = NULL;

/*
 * WARNING: the sign of all MPFR numbers will be permanently overwritten by the
 * arpra_helper_mpfr_sumabs function! Signs should be backed up first, if needed.
 */

int arpra_helper_mpfr_sumabs (arpra_mpfr *sumabs, arpra_mpfr *x,
                              const arpra_uint n, const mpfr_rnd_t round_mode)
{
    arpra_uint i;

    // Allocate or resize pointer_buffer, as required.
    if (pointer_buffer_size < n) {
        pointer_buffer_size = ceil(n / ARPRA_BUFFER_SIZE_FACTOR);
        pointer_buffer_size *= ARPRA_BUFFER_SIZE_FACTOR;
        pointer_buffer = realloc(pointer_buffer, pointer_buffer_size * sizeof(arpra_mpfr *));
    }

    // Save term addresses to pointer_buffer and compute abs.
    for (i = 0; i < n; i++) {
        pointer_buffer[i] = &(x[i]);
        mpfr_abs(pointer_buffer[i], pointer_buffer[i], MPFR_RNDN);
    }

    // Sum the terms pointed to in pointer_buffer.
    return mpfr_sum(sumabs, pointer_buffer, n, round_mode);
}

void arpra_helper_radius (arpra_range *z)
{
    arpra_uint zTerm;

    // Allocate or resize sign_buffer, as required.
    if (sign_buffer_size < z->nTerms) {
        sign_buffer_size = ceil(z->nTerms / ARPRA_BUFFER_SIZE_FACTOR);
        sign_buffer_size *= ARPRA_BUFFER_SIZE_FACTOR;
        sign_buffer = realloc(sign_buffer, sign_buffer_size * sizeof(arpra_sign));
    }

    // Save term signs to sign_buffer.
    for (zTerm = 0; zTerm < z->nTerms; zTerm++) {
        sign_buffer[zTerm] = z->deviations[zTerm]._mpfr_sign;
    }

    // Compute sum of absolute value terms.
    arpra_helper_mpfr_sumabs(&(z->radius), z->deviations, z->nTerms, MPFR_RNDU);

    // Restore term signs from sign_buffer.
    for (zTerm = 0; zTerm < z->nTerms; zTerm++) {
        z->deviations[zTerm]._mpfr_sign = sign_buffer[zTerm];
    }
}
