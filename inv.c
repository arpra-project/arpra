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
	mpfr_t temp, xl, xh, dl, dh, alpha, gamma, delta;
	mpfr_prec_t prec;

	prec = mpfr_get_prec(&(z->centre));
	mpfr_init2(temp, prec);
	mpfr_init2(xl, prec);
	mpfr_init2(xh, prec);
	mpfr_init2(dl, prec);
	mpfr_init2(dh, prec);
	mpfr_init2(alpha, prec);
	mpfr_init2(gamma, prec);
	mpfr_init2(delta, prec);

	mpfr_sub(xl, &(x->centre), &(x->radius), MPFR_RNDD);
	mpfr_add(xh, &(x->centre), &(x->radius), MPFR_RNDU);
	if (mpfr_sgn(xl) != mpfr_sgn(xh)) {
		// TODO: handle singularity
		// return ANY_NUM;
	}

	// compute alpha
	mpfr_si_div(alpha, -1, xh, MPFR_RNDD);
	mpfr_div(alpha, alpha, xl, MPFR_RNDD);

	// compute max and min differences (1/x - alpha x)
	mpfr_mul(dl, alpha, xl, MPFR_RNDD);
	mpfr_si_div(temp, 1, xl, MPFR_RNDU);
	mpfr_sub(dl, temp, dl, MPFR_RNDU);
	mpfr_mul(dh, alpha, xh, MPFR_RNDD);
	mpfr_si_div(temp, 1, xh, MPFR_RNDU);
	mpfr_sub(dh, temp, dh, MPFR_RNDU);
	mpfr_max(dh, dl, dh, MPFR_RNDN);
	mpfr_neg(dl, alpha, MPFR_RNDN);
	mpfr_sqrt(dl, dl, MPFR_RNDD);
	mpfr_mul_si(dl, dl, 2, MPFR_RNDD);

	// compute gamma
	mpfr_add(gamma, dl, dh, MPFR_RNDN);
	mpfr_div_si(gamma, gamma, 2, MPFR_RNDN);

	// compute delta
	mpfr_sub(delta, dh, gamma, MPFR_RNDU);
	mpfr_sub(temp, gamma, dl, MPFR_RNDU);
	mpfr_max(delta, delta, temp, MPFR_RNDN);

	// compute affine approximation
	mpfa_affine_1(z, x, alpha, gamma, delta);

	mpfr_clear(temp);
	mpfr_clear(xl);
	mpfr_clear(xh);
	mpfr_clear(dl);
	mpfr_clear(dh);
	mpfr_clear(alpha);
	mpfr_clear(gamma);
	mpfr_clear(delta);
}
