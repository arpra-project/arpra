/*
 * helper_resize_term_memory.c -- Resize a range's deviation term arrays.
 *
 * Copyright 2019 James Paul Turner.
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

void arpra_helper_resize_term_memory (arpra_range *z, arpra_uint n)
{
    arpra_uint zTerm;

    // Clear existing deviation terms.
    for (zTerm = 0; zTerm < z->nTerms; zTerm++) {
        mpfr_clear(&(z->deviations[zTerm]));
    }
    z->nTerms = 0;

    // Free term memory.
    if ((n == 0) && (z->term_memory > 0)) {
        free(z->symbols);
        free(z->deviations);
        z->term_memory = 0;
    }

    // Allocate term memory.
    else if ((n > 0) && (z->term_memory == 0)) {
        z->symbols = malloc(n * sizeof(arpra_uint));
        z->deviations = malloc(n * sizeof(arpra_mpfr));
        z->term_memory = n;
    }

    // Reallocate term memory.
    else if (n != z->term_memory) {
        z->symbols = realloc(z->symbols, n * sizeof(arpra_uint));
        z->deviations = realloc(z->deviations, n * sizeof(arpra_mpfr));
        z->term_memory = n;
    }
}
