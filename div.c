/*
 * div.c
 *
 *  Created on: 16 Nov 2016
 *      Author: jt273
 */

#include "mpfa.h"
#include <malloc.h>
#include <assert.h>

void mpfa_div (mpfa_ptr z, mpfa_srcptr x, mpfa_srcptr y) {
	unsigned xTerm, yTerm, zTerm;
	int inexact;
	mpfr_prec_t prec;
	mpfr_t u, temp, error, delta;
	mpfa_t z_new;

	prec = mpfr_get_prec(&(z->centre));
	mpfr_init2(u, prec);
	mpfr_init2(temp, prec);
	mpfr_init2(error, prec);
	mpfr_init2(delta, prec);
	mpfr_init2(&(z_new->centre), prec);
	mpfr_init2(&(z_new->radius), prec);
	mpfr_set_d(&(z_new->radius), 0.0, MPFR_RNDN);

	assert(!mpfr_set_si(u, -prec, MPFR_RNDN)); // fails if emax <= log2(prec)
	assert(!mpfr_exp2(u, u, MPFR_RNDN)); // fails if emin > 1-prec

	inexact = mpfr_div(&(z_new->centre), &(x->centre), &(y->centre), MPFR_RNDN);
	if (inexact) {
		mpfr_mul(delta, u, &(z_new->centre), MPFR_RNDU);
	}

	z_new->nTerms = x->nTerms + y->nTerms + 1;
	z_new->symbols = malloc(z_new->nTerms * sizeof(unsigned));
	z_new->deviations = malloc(z_new->nTerms * sizeof(mpfr_t));









	mpfr_clear(u);
	mpfr_clear(temp);
	mpfr_clear(error);
	mpfr_clear(delta);
	mpfa_clear(z);
	*z = *z_new;
}
