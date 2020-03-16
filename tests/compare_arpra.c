/*
 * compare_arpra.c -- Check the number of differences in two arpra variables.
 *
 * Copyright 2017-2020 James Paul Turner.
 *
 * This file is part of the  library.
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

#include "arpra-test.h"

int test_compare_arpra (const arpra_range *x1, const arpra_range *x2)
{
    int fail = 0;
    arpra_uint i;

    // Return the number of differences in x1 and x2.
    if (!mpfr_equal_p(&(x1->centre), &(x2->centre))) fail++;
    if (!mpfr_equal_p(&(x1->radius), &(x2->radius))) fail++;
    if (x1->nTerms != x2->nTerms) fail++;
    for (i = 0; i < x1->nTerms; i++) {
        if (x1->symbols[i] != x2->symbols[i]) fail++;
        if (!mpfr_equal_p(&(x1->deviations[i]), &(x2->deviations[i]))) fail++;
    }

    return fail;
}
