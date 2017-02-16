/*
 * inv.c
 *
 *  Created on: 14 Feb 2017
 *      Author: jt273
 */

#include "mpfa.h"

void mpfa_inv (mpfa_ptr z, mpfa_srcptr x) {
	mpfr_t temp1, temp2, xl, xh, alpha, gamma, delta;
	mpfr_prec_t prec;
	int negative;

	prec = mpfr_get_prec(&(z->centre));
	mpfr_init2(temp1, prec);
	mpfr_init2(temp2, prec);
	mpfr_init2(xl, prec);
	mpfr_init2(xh, prec);
	mpfr_init2(alpha, prec);
	mpfr_init2(gamma, prec);
	mpfr_init2(delta, prec);

	mpfr_sub(xl, &(x->centre), &(x->radius), MPFR_RNDD);
	mpfr_add(xh, &(x->centre), &(x->radius), MPFR_RNDU);
	if (mpfr_sgn(xl) != mpfr_sgn(xh)) {
		// TODO: handle singularity
		// return ANY_NUM;
	}
	negative = mpfr_sgn(xl) < 0;
	if (negative) {
		mpfr_swap(xl, xh);
	}

	// compute alpha
	mpfr_si_div(alpha, -1, xh, MPFR_RNDD);
	mpfr_div(alpha, alpha, xl, MPFR_RNDD);

	// compute gamma


	// compute delta


	mpfr_clear(temp1);
	mpfr_clear(temp2);
	mpfr_clear(xl);
	mpfr_clear(xh);
	mpfr_clear(alpha);
	mpfr_clear(gamma);
	mpfr_clear(delta);
}
