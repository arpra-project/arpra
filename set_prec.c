/*
 * set_prec.c
 *
 *  Created on: 22 Sep 2016
 *      Author: jt273
 */

#include "mpfa.h"

void mpfa_set_prec (mpfa_ptr x, mpfr_prec_t prec) {
	unsigned xTerm;
	for (xTerm = 0; xTerm < x->nTerms; xTerm++) {
		mpfr_set_prec(&(x->deviations[xTerm]), prec);
	}
	mpfr_set_prec(&(x->centre), prec);
	mpfr_set_prec(&(x->radius), prec);
}
