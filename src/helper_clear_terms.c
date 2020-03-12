/*
 * helper_clear_terms.c -- Clear and free deviation term arrays.
 *
 * Copyright 2019-2020 James Paul Turner.
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

void arpra_helper_clear_terms (arpra_range *y)
{
    arpra_uint iy;

    if (y->nTerms > 0) {
        for (iy = 0; iy < y->nTerms; iy++) {
            mpfr_clear(&(y->deviations[iy]));
        }
        free(y->symbols);
        free(y->deviations);
    }
}
