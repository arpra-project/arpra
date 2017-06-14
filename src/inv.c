/*
 * inv.c -- Compute the inverse of an affine form.
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

#include "mpfa.h"
#include <stdlib.h>

/*
 * This affine inverse function uses a Chebyshev linear approximation.
 */

void mpfa_inv (mpfa_ptr z, mpfa_srcptr x) {
    int sign;
    mpfr_t temp, xa, xb, da, db, du, alpha, gamma, delta;
    mpfr_prec_t prec;

    prec = mpfr_get_prec(&(z->centre));
    mpfr_inits2(prec, temp, xa, xb, da, db, du,
                alpha, gamma, delta, (mpfr_ptr) NULL);

    if (mpfr_zero_p(&(x->radius))) {
        if (mpfr_si_div(temp, 1, &(x->centre), MPFR_RNDN)) {
            mpfa_error(delta, temp);
        }
        else {
            mpfr_set_si(delta, 0, MPFR_RNDN);
        }

        mpfa_set_mpfr_rad(z, temp, delta);
    }
    else {
        mpfr_sub(xa, &(x->centre), &(x->radius), MPFR_RNDD);
        mpfr_add(xb, &(x->centre), &(x->radius), MPFR_RNDU);

        sign = mpfr_sgn(xa);
        if ((sign != mpfr_sgn(xb) || (sign == 0))) {
            if (z->nTerms > 0) {
                unsigned zTerm;
                for (zTerm = 0; zTerm < z->nTerms; zTerm++) {
                    mpfr_clear(&(z->deviations[zTerm]));
                }
                z->nTerms = 0;
                free(z->symbols);
                free(z->deviations);
            }

            // TODO: find a better representation for Inf
            mpfr_set_nan(&(z->centre));
        }
        else {
            if (sign < 0) {
                mpfr_set(temp, xa, MPFR_RNDN);
                mpfr_neg(xa, xb, MPFR_RNDN);
                mpfr_neg(xb, temp, MPFR_RNDN);
            }

            // compute alpha
            mpfr_si_div(alpha, -1, xb, MPFR_RNDN);
            mpfr_div(alpha, alpha, xa, MPFR_RNDN);

            // compute difference (1/a - alpha a)
            mpfr_mul(da, alpha, xa, MPFR_RNDD);
            mpfr_si_div(temp, 1, xa, MPFR_RNDU);
            mpfr_sub(da, temp, da, MPFR_RNDU);

            // compute difference (1/b - alpha b)
            mpfr_mul(db, alpha, xb, MPFR_RNDD);
            mpfr_si_div(temp, 1, xb, MPFR_RNDU);
            mpfr_sub(db, temp, db, MPFR_RNDU);

            mpfr_max(db, da, db, MPFR_RNDN);

            // compute difference (1/u - alpha u)
            mpfr_neg(du, alpha, MPFR_RNDN);
            mpfr_sqrt(du, du, MPFR_RNDD);
            mpfr_mul_si(du, du, 2, MPFR_RNDD);

            // compute gamma
            mpfr_add(gamma, db, du, MPFR_RNDN);
            mpfr_div_si(gamma, gamma, 2, MPFR_RNDN);

            // compute delta
            mpfr_sub(delta, gamma, du, MPFR_RNDU);
            mpfr_sub(temp, db, gamma, MPFR_RNDU);
            mpfr_max(delta, delta, temp, MPFR_RNDN);

            if (sign < 0) {
                mpfr_neg(gamma, gamma, MPFR_RNDN);
            }

            // compute affine approximation
            mpfa_affine_1(z, x, alpha, gamma, delta);
        }
    }

    mpfr_clears(temp, xa, xb, da, db, du,
                alpha, gamma, delta, (mpfr_ptr) NULL);
}
