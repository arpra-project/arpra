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
    mpfa_uint_t term;

    // Return 1 if x and y differ.
    if (!mpfr_equal_p(&(x->centre), &(y->centre))) {
        return 1;
    }

    if (!mpfr_equal_p(&(x->radius), &(y->radius))) {
        return 1;
    }

    if (x->nTerms != y->nTerms) {
        return 1;
    }

    for (term = 0; term < x->nTerms; term++) {
        if (x->symbols[term] != y->symbols[term]) {
            return 1;
        }

        if (!mpfr_equal_p(&(x->deviations[term]), &(y->deviations[term]))) {
            return 1;
        }
    }

    // Else return 0.
    return 0;
}
