/*
 * compare_arpra.c -- Check the number of differences in two arpra variables.
 *
 * Copyright 2017-2018 James Paul Turner.
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

int test_compare_arpra (arpra_srcptr x, arpra_srcptr y)
{
    int fail = 0;
    arpra_uint_t term;

    // Return the number of differences in x and y.
    if (!mpfr_equal_p(&(x->centre), &(y->centre))) fail++;
    if (!mpfr_equal_p(&(x->radius), &(y->radius))) fail++;
    if (x->nTerms != y->nTerms) fail++;
    for (term = 0; term < x->nTerms; term++) {
        if (x->symbols[term] != y->symbols[term]) fail++;
        if (!mpfr_equal_p(&(x->deviations[term]), &(y->deviations[term]))) fail++;
    }

    return fail;
}
