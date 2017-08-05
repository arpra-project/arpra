/*
 * prec.c -- Get and set the precision of an affine form.
 *
 * Copyright 2016-2017 James Paul Turner.
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

#include "mpfa.h"

mpfa_prec_t mpfa_get_prec (mpfa_srcptr x) {
    return mpfr_get_prec(&(x->centre));
}

void mpfa_set_prec (mpfa_ptr x, mpfa_prec_t prec) {
    mpfa_uint_t xTerm;

    if (mpfa_get_internal_prec() < prec) {
        mpfa_set_internal_prec(prec);
    }

    for (xTerm = 0; xTerm < x->nTerms; xTerm++) {
        mpfr_set_prec(&(x->deviations[xTerm]), prec);
    }

    mpfr_set_prec(&(x->centre), prec);
    mpfr_set_prec(&(x->radius), prec);
}
