/*
 * clear.c -- Free the memory of one or more Arpra ranges.
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

void arpra_clear (arpra_range *z)
{
    // Clear centre, radius and deviation terms.
    mpfr_clear(&(z->centre));
    mpfr_clear(&(z->radius));
    arpra_clear_terms(z);
}

void arpra_clears (arpra_range *z, ...)
{
    va_list arg;

    // Clear each arguemnt.
    va_start(arg, z);
    while (z != NULL) {
        arpra_clear(z);
        z = (arpra_range *) va_arg(arg, arpra_range *);
    }
    va_end(arg);
}

void arpra_clear_terms (arpra_range *z)
{
    arpra_uint zTerm;

    // Clear ezisting deviation terms.
    if (z->nTerms > 0) {
        for (zTerm = 0; zTerm < z->nTerms; zTerm++) {
            mpfr_clear(&(z->deviations[zTerm]));
        }
        free(z->symbols);
        free(z->deviations);
        z->nTerms = 0;
    }
}
