/*
 * inv.c
 *
 *  Created on: 14 Feb 2017
 *      Author: jt273
 */

#include "mpfa.h"

/*
 * This affine inverse function uses a Chebyshev linear approximation.
 */

void mpfa_inv (mpfa_ptr z, mpfa_srcptr x) {
	int sign;
	mpfr_t temp, xa, xb, da, db, du, alpha, gamma, delta;
	mpfr_prec_t prec;

	prec = mpfr_get_prec(&(z->centre));
	mpfr_init2(temp, prec);
	mpfr_init2(xa, prec);
	mpfr_init2(xb, prec);
	mpfr_init2(da, prec);
	mpfr_init2(db, prec);
	mpfr_init2(du, prec);
	mpfr_init2(alpha, prec);
	mpfr_init2(gamma, prec);
	mpfr_init2(delta, prec);

	mpfr_sub(xa, &(x->centre), &(x->radius), MPFR_RNDD);
	mpfr_add(xb, &(x->centre), &(x->radius), MPFR_RNDU);
	sign = mpfr_sgn(xa);
	if ((sign != mpfr_sgn(xb) || (sign == 0))) {
		// TODO: handle singularity
		// return ANY_NUM;
	}

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

	mpfr_max(da, da, db, MPFR_RNDN);

	// compute difference (1/u - alpha u)
	mpfr_neg(du, alpha, MPFR_RNDN);
	mpfr_sqrt(du, du, MPFR_RNDD);
	mpfr_mul_si(du, du, 2, MPFR_RNDD);

	// compute gamma
	mpfr_add(gamma, da, du, MPFR_RNDN);
	mpfr_div_si(gamma, gamma, 2, MPFR_RNDN);

	// compute delta
	mpfr_sub(delta, du, gamma, MPFR_RNDU);
	mpfr_sub(temp, gamma, da, MPFR_RNDU);
	mpfr_max(delta, delta, temp, MPFR_RNDN);

	if (sign < 0) {
		mpfr_neg(gamma, gamma, MPFR_RNDN);
	}

	// compute affine approximation
	mpfa_affine_1(z, x, alpha, gamma, delta);

	mpfr_clear(temp);
	mpfr_clear(xa);
	mpfr_clear(xb);
	mpfr_clear(da);
	mpfr_clear(db);
	mpfr_clear(du);
	mpfr_clear(alpha);
	mpfr_clear(gamma);
	mpfr_clear(delta);
}
