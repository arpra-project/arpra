/*
 * neg.c
 *
 *  Created on: 11 Nov 2016
 *      Author: jt273
 */

#include "mpfa.h"

void mpfa_neg (mpfa_ptr z, mpfa_srcptr x) {
	mpfr_t alpha, gamma, delta;
	mpfr_prec_t prec;

	prec = mpfr_get_prec(&(z->centre));
	mpfr_init2(alpha, prec);
	mpfr_set_si(alpha, -1, MPFR_RNDN);
	mpfr_init2(gamma, prec);
	mpfr_set_si(gamma, 0, MPFR_RNDN);
	mpfr_init2(delta, prec);
	mpfr_set_si(delta, 0, MPFR_RNDN);

	mpfa_affine_1(z, x, alpha, gamma, delta);

	mpfr_clear(alpha);
	mpfr_clear(gamma);
	mpfr_clear(delta);
}
