/*
 * sub.c
 *
 *  Created on: 11 Nov 2016
 *      Author: jt273
 */

#include "mpfa.h"

void mpfa_sub (mpfa_ptr z, mpfa_srcptr x, mpfa_srcptr y) {
	mpfr_t alpha, beta, gamma, delta;
	mpfr_prec_t prec;

	prec = mpfr_get_prec(&(z->centre));
	mpfr_init2(alpha, prec);
	mpfr_set_si(alpha, 1, MPFR_RNDN);
	mpfr_init2(beta, prec);
	mpfr_set_si(beta, -1, MPFR_RNDN);
	mpfr_init2(gamma, prec);
	mpfr_set_si(gamma, 0, MPFR_RNDN);
	mpfr_init2(delta, prec);
	mpfr_set_si(delta, 0, MPFR_RNDN);

	mpfa_affine_2(z, alpha, x, beta, y, gamma, delta);

	mpfr_clear(alpha);
	mpfr_clear(beta);
	mpfr_clear(gamma);
	mpfr_clear(delta);
}
