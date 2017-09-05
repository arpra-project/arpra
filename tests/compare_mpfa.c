/*
 * compare_mpfa.c -- Test an MPFA result against another MPFA result.
 *
 * Copyright 2017 James Paul Turner.
 *
 * This file is part of the MPFA library.
 *
 * The MPFA library is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * The MPFA library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with the MPFA library. If not, see <http://www.gnu.org/licenses/>.
 */

#include "mpfa-test.h"

mpfa_int_t test_compare_mpfa (mpfa_srcptr x, mpfa_srcptr y)
{
    mpfa_int_t fail;
    mpfa_uint_t term;

    // Return the number of differences in x and y.
    fail = 0;
    if (!mpfr_equal_p(&(x->centre), &(y->centre))) fail++;
    if (!mpfr_equal_p(&(x->radius), &(y->radius))) fail++;
    if (x->nTerms != y->nTerms) fail++;
    for (term = 0; term < x->nTerms; term++) {
        if (x->symbols[term] != y->symbols[term]) fail++;
        if (!mpfr_equal_p(&(x->deviations[term]), &(y->deviations[term]))) fail++;
    }

    return fail;
}
