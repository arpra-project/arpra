/*
 * term.c
 *
 *  Created on: 20 Mar 2017
 *      Author: james
 */

#include "mpfa.h"
#include <assert.h>

int mpfa_term (mpfr_ptr z, mpfr_srcptr x, mpfr_srcptr y, mpfr_srcptr alpha, mpfr_srcptr beta, mpfr_srcptr gamma) {
	int inexact;
	mpfr_t alpha_x, beta_y;

	mpfr_init2(alpha_x, (mpfr_get_prec(alpha) + mpfr_get_prec(x)));
	assert(mpfr_mul(alpha_x, alpha, x, MPFR_RNDN) == 0);

	mpfr_init2(beta_y, (mpfr_get_prec(beta) + mpfr_get_prec(y)));
	assert(mpfr_mul(beta_y, beta, y, MPFR_RNDN) == 0);

	if (gamma == NULL) {
		inexact = mpfr_add(z, alpha_x, beta_y, MPFR_RNDN);
	}
	else {
		mpfr_ptr sumArray[3];
		sumArray[0] = alpha_x;
		sumArray[1] = beta_y;
		sumArray[2] = (mpfr_ptr) gamma;
		inexact = mpfr_sum(z, sumArray, 3, MPFR_RNDN);
	}

	mpfr_clears(alpha_x, beta_y, (mpfr_ptr) NULL);
	return inexact;
}
