/*
 * exp.c -- Compute the exponent of an affine form.
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

/*
 * This affine exponential function uses a Chebyshev linear approximation.
 */

void mpfa_exp (mpfa_ptr z, mpfa_srcptr x) {
    mpfr_t temp, xa, xb, da, db, du, alpha, gamma, delta;
    mpfr_prec_t prec;

    prec = mpfr_get_prec(&(z->centre));
    mpfr_inits2(prec, temp, xa, xb, da, db, du,
		alpha, gamma, delta, (mpfr_ptr) NULL);

    if (mpfr_zero_p(&(x->radius))) {
        if (mpfr_exp(temp, &(x->centre), MPFR_RNDN)) {
            mpfr_mul(delta, temp, &(z->u), MPFR_RNDU);
        }
        else {
            mpfr_set_si(delta, 0, MPFR_RNDN);
        }

        mpfa_set_mpfr(z, temp, delta);
    }
    else {
        mpfr_sub(xa, &(x->centre), &(x->radius), MPFR_RNDD);
        mpfr_add(xb, &(x->centre), &(x->radius), MPFR_RNDU);

        // compute alpha
        mpfr_exp(alpha, xb, MPFR_RNDN);
        mpfr_exp(temp, xa, MPFR_RNDN);
        mpfr_sub(alpha, alpha, temp, MPFR_RNDN);
        mpfr_sub(temp, xb, xa, MPFR_RNDN);
        mpfr_div(alpha, alpha, temp, MPFR_RNDN);

        // compute difference (exp(a) - alpha a)
        mpfr_mul(da, alpha, xa, MPFR_RNDD);
        mpfr_exp(temp, xa, MPFR_RNDU);
        mpfr_sub(da, temp, da, MPFR_RNDU);

        // compute difference (exp(b) - alpha b)
        mpfr_mul(db, alpha, xb, MPFR_RNDD);
        mpfr_exp(temp, xb, MPFR_RNDU);
        mpfr_sub(db, temp, db, MPFR_RNDU);

        mpfr_max(da, da, db, MPFR_RNDN);

        // compute difference (exp(u) - alpha u)
        mpfr_log(du, alpha, MPFR_RNDU);
        mpfr_sub_si(du, du, 1, MPFR_RNDU);
        mpfr_mul(du, alpha, du, MPFR_RNDU);
        mpfr_neg(du, du, MPFR_RNDN);

        // compute gamma
        mpfr_add(gamma, da, du, MPFR_RNDN);
        mpfr_div_si(gamma, gamma, 2, MPFR_RNDN);

        // compute delta
        mpfr_sub(delta, du, gamma, MPFR_RNDU);
        mpfr_sub(temp, gamma, da, MPFR_RNDU);
        mpfr_max(delta, delta, temp, MPFR_RNDN);

        // compute affine approximation
        mpfa_affine_1(z, x, alpha, gamma, delta);
    }

    mpfr_clears(temp, xa, xb, da, db, du,
		alpha, gamma, delta, (mpfr_ptr) NULL);
}
