/*
 * test_mpfa.c -- Compare function results against known affine forms.
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

int mpfa_test_mpfa_1 (void (*mpfa_1) (mpfa_ptr z, mpfa_srcptr x),
                      mpfa_ptr expect, mpfa_srcptr x)
{
    int success;
    mpfa_prec_t prec;
    mpfa_t z;

    // Init test var with precision of expected result.
    prec = mpfa_get_prec(expect);
    mpfa_init2(z, prec);

    // Compare test result with expected result.
    mpfa_1(z, x);
    mpfa_test_cmp_mpfa(z, expect);

    // Clear test var and return test result.
    mpfa_clear(z);
    return success;
}

int mpfa_test_mpfa_2 (void (*mpfa_2) (mpfa_ptr z, mpfa_srcptr x, mpfa_srcptr y),
                      mpfa_ptr expect, mpfa_srcptr x, mpfa_srcptr y)
{
    int success;
    mpfa_prec_t prec;
    mpfa_t z;

    // Init test var with precision of expected result.
    prec = mpfa_get_prec(expect);
    mpfa_init2(z, prec);

    // Compare test result with expected result.
    mpfa_2(z, x, y);
    mpfa_test_cmp_mpfa(z, expect);

    // Clear test var and return test result.
    mpfa_clear(z);
    return success;
}

int mpfa_test_cmp_mpfa (mpfa_srcptr got, mpfa_srcptr expect)
{
    mpfa_uint_t term;

    // Compare term count, centre and radius.
    if (got->nTerms != expect->nTerms) {
        return 1;
    }
    if (!mpfr_equal_p(got->centre, expect->centre)) {
        return 1;
    }
    if (!mpfr_equal_p(got->radius, expect->radius)) {
        return 1;
    }

    // Compare all deviation terms.
    for (term = 0; term < got->nTerms; term++) {
        if (got->symbols[term] != expect->symbols[term]) {
            return 1;
        }
        if (!mpfr_equal_p(got->deviations[term], expect->deviations[term])) {
            return 1;
        }
    }

    return 0;
}
