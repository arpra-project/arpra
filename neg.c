/*
 * neg.c
 *
 *  Created on: 11 Nov 2016
 *      Author: jt273
 */

#include "mpfa.h"

void mpfa_neg (mpfa_ptr z, mpfa_srcptr x) {
	mpfr_t alpha, gamma, delta;

	mpfr_init_set_d(alpha, -1.0, MPFR_RNDN);
	mpfr_init_set_d(gamma, 0.0, MPFR_RNDN);
	mpfr_init_set_d(delta, 0.0, MPFR_RNDN);

	mpfa_affine_1(z, alpha, x, gamma, delta);

	mpfr_clear(alpha);
	mpfr_clear(gamma);
	mpfr_clear(delta);
}
