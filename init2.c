/*
 * init2.c
 *
 *  Created on: 21 Feb 2017
 *      Author: james
 */

#include "mpfa.h"
#include <malloc.h>

void mpfa_init2 (mpfa_ptr x, mpfr_prec_t prec) {
	x->nTerms = 0;
	x->symbols = malloc(0);
	x->deviations = malloc(0);
	mpfr_init2(&(x->centre), prec);
	mpfr_init2(&(x->radius), prec);
	mpfr_set_si(&(x->radius), 0, MPFR_RNDN);
}
