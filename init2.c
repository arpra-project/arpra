/*
 * init2.c
 *
 *  Created on: 21 Feb 2017
 *      Author: james
 */

#include "mpfa.h"

void mpfa_init2 (mpfa_ptr x, mpfr_prec_t prec) {
	x->nTerms = 0;
	mpfr_init2(&(x->centre), prec);
	mpfr_init2(&(x->radius), prec);
}
