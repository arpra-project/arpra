/*
 * exp.c
 *
 *  Created on: 21 Feb 2017
 *      Author: james
 */

#include "mpfa.h"
#include <assert.h>

/*
 * This affine exponential function uses a Chebyshev linear approximation.
 */

void mpfa_exp (mpfa_ptr z, mpfa_srcptr x) {
	mpfr_t temp, xa, xb, da, db, du, alpha, gamma, delta;
	mpfr_prec_t prec;

	prec = mpfr_get_prec(&(z->centre));

	if (mpfr_zero_p(&(x->radius))) {
		mpfr_inits2(prec, temp, delta, (mpfr_ptr) NULL);

		if (mpfr_exp(temp, &(x->centre), MPFR_RNDN)) {
			assert(mpfr_set_si(delta, (-prec + mpfr_get_exp(temp)), MPFR_RNDN) == 0);
			assert(mpfr_exp2(delta, delta, MPFR_RNDN) == 0);
		}
		else {
			mpfr_set_si(delta, 0, MPFR_RNDN);
		}

		mpfa_set_mpfr(z, temp, delta);

		mpfr_clears(temp, delta, (mpfr_ptr) NULL);
	}
	else {
		mpfr_inits2(prec, temp, xa, xb, da, db, du, alpha, gamma, delta, (mpfr_ptr) NULL);

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

		mpfr_clears(temp, xa, xb, da, db, du, alpha, gamma, delta, (mpfr_ptr) NULL);
	}
}
