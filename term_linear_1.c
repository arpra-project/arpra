/*
 * affine_term_1.c
 *
 *  Created on: 2 Nov 2016
 *      Author: jt273
 */

#include "mpfa.h"

void mpfa_term_linear_1 (mpfr_ptr z, mpfr_srcptr alpha, mpfr_srcptr x, mpfr_srcptr gamma, mpfr_ptr delta) {
	mpfr_t ax, lo, hi, error;
	mpfr_init(ax);
	mpfr_init(lo);
	mpfr_init(hi);
	mpfr_init(error);

	// alpha * x
	mpfr_mul(ax, alpha, x, MPFR_RNDN);
	mpfr_mul(lo, alpha, x, MPFR_RNDD);
	mpfr_mul(hi, alpha, x, MPFR_RNDU);
	mpfr_sub(lo, ax, lo, MPFR_RNDU);
	mpfr_sub(hi, hi, ax, MPFR_RNDU);
	mpfr_max(error, lo, hi, MPFR_RNDU);
	mpfr_add(delta, delta, error, MPFR_RNDU);


	// alpha * x + gamma
	if (gamma) {
		mpfr_add(z, ax, gamma, MPFR_RNDN);
		mpfr_add(lo, ax, gamma, MPFR_RNDD);
		mpfr_add(hi, ax, gamma, MPFR_RNDU);
		mpfr_sub(lo, z, lo, MPFR_RNDU);
		mpfr_sub(hi, hi, z, MPFR_RNDU);
		mpfr_max(error, lo, hi, MPFR_RNDU);
		mpfr_add(delta, delta, error, MPFR_RNDU);
	}
	else {
		mpfr_set(z, ax, MPFR_RNDN);
	}

	mpfr_clear(ax);
	mpfr_clear(lo);
	mpfr_clear(hi);
	mpfr_clear(error);
}
