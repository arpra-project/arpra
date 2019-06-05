/*
 * helper_term_arrays.c -- Allocate and free deviation term arrays.
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

void arpra_helper_alloc_term_arrays (arpra_range *z, arpra_uint n)
{
    if (n > 0) {
        z->symbols = malloc(n * sizeof(arpra_uint));
        z->deviations = malloc(n * sizeof(arpra_mpfr));
    }
}

void arpra_helper_free_term_arrays (arpra_range *z)
{
    arpra_uint zTerm;

    // Clear existing deviation terms.
    for (zTerm = 0; zTerm < z->nTerms; zTerm++) {
        mpfr_clear(&(z->deviations[zTerm]));
    }

    if (z->nTerms > 0) {
        free(z->symbols);
        free(z->deviations);
    }
}
