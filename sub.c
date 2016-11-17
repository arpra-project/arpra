/*
 * sub.c
 *
 *  Created on: 11 Nov 2016
 *      Author: jt273
 */

#include "mpfa.h"

void mpfa_sub (mpfa_ptr z, mpfa_srcptr x, mpfa_srcptr y) {
	mpfr_t alpha, beta, gamma, delta;

	mpfr_init_set_d(alpha, 1.0, MPFR_RNDN);
	mpfr_init_set_d(beta, -1.0, MPFR_RNDN);
	mpfr_init_set_d(gamma, 0.0, MPFR_RNDN);
	mpfr_init_set_d(delta, 0.0, MPFR_RNDN);

	mpfa_affine_2(z, alpha, x, beta, y, gamma, delta);

	mpfr_clear(alpha);
	mpfr_clear(beta);
	mpfr_clear(gamma);
	mpfr_clear(delta);
}
