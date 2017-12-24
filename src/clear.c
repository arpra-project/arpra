/*
 * clear.c -- Free the memory of one or more arpra_t.
 *
 * Copyright 2016-2017 James Paul Turner.
 *
 * This file is part of the ArPRA library.
 *
 * The ArPRA library is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * The ArPRA library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with the ArPRA library. If not, see <http://www.gnu.org/licenses/>.
 */

#include "arpra-impl.h"

void arpra_clear (arpra_ptr x)
{
    // Clear centre, radius and deviation terms.
    mpfr_clear(&(x->centre));
    mpfr_clear(&(x->radius));
    arpra_clear_terms(x);
}

void arpra_clears (arpra_ptr x, ...)
{
    va_list arg;

    // Clear each arguemnt.
    va_start(arg, x);
    while (x != NULL) {
        arpra_clear(x);
        x = (arpra_ptr) va_arg(arg, arpra_ptr);
    }
    va_end(arg);
}

void arpra_clear_terms (arpra_ptr x)
{
    arpra_uint_t xTerm;

    // Clear existing deviation terms.
    if (x->nTerms > 0) {
        for (xTerm = 0; xTerm < x->nTerms; xTerm++) {
            mpfr_clear(&(x->deviations[xTerm]));
        }
        free(x->symbols);
        free(x->deviations);
        x->nTerms = 0;
    }
}
