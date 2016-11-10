/*
 * affine_term_2.c
 *
 *  Created on: 2 Nov 2016
 *      Author: jt273
 */

#include "mpfa.h"

void mpfa_affine_term_2 (mpfr_ptr z, mpfr_srcptr alpha, mpfr_srcptr x, mpfr_srcptr beta, mpfr_srcptr y, mpfr_srcptr gamma, mpfr_ptr delta) {
	mpfr_t ax, by, ax_by, lo, hi, error;
	mpfr_init(ax);
	mpfr_init(by);
	mpfr_init(ax_by);
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

	// beta * y
	mpfr_mul(by, beta, y, MPFR_RNDN);
	mpfr_mul(lo, beta, y, MPFR_RNDD);
	mpfr_mul(hi, beta, y, MPFR_RNDU);
	mpfr_sub(lo, by, lo, MPFR_RNDU);
	mpfr_sub(hi, hi, by, MPFR_RNDU);
	mpfr_max(error, lo, hi, MPFR_RNDU);
	mpfr_add(delta, delta, error, MPFR_RNDU);

	// alpha * x + beta * y
	mpfr_add(ax_by, ax, by, MPFR_RNDN);
	mpfr_add(lo, ax, by, MPFR_RNDD);
	mpfr_add(hi, ax, by, MPFR_RNDU);
	mpfr_sub(lo, ax_by, lo, MPFR_RNDU);
	mpfr_sub(hi, hi, ax_by, MPFR_RNDU);
	mpfr_max(error, lo, hi, MPFR_RNDU);
	mpfr_add(delta, delta, error, MPFR_RNDU);

	// alpha * x + beta + y + gamma
	if (gamma) {
		mpfr_add(z, ax_by, gamma, MPFR_RNDN);
		mpfr_add(lo, ax_by, gamma, MPFR_RNDD);
		mpfr_add(hi, ax_by, gamma, MPFR_RNDU);
		mpfr_sub(lo, z, lo, MPFR_RNDU);
		mpfr_sub(hi, hi, z, MPFR_RNDU);
		mpfr_max(error, lo, hi, MPFR_RNDU);
		mpfr_add(delta, delta, error, MPFR_RNDU);
	}
	else {
		mpfr_set(z, ax_by, MPFR_RNDN);
	}

	mpfr_clear(ax);
	mpfr_clear(by);
	mpfr_clear(ax_by);
	mpfr_clear(lo);
	mpfr_clear(hi);
	mpfr_clear(error);
}
