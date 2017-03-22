/*
 * clear.c
 *
 *  Created on: 20 Sep 2016
 *      Author: james
 */

#include "mpfa.h"
#include <malloc.h>

void mpfa_clear (mpfa_ptr x) {
	unsigned xTerm;

	if (x->nTerms > 0) {
		for (xTerm = 0; xTerm < x->nTerms; xTerm++) {
			mpfr_clear(&(x->deviations[xTerm]));
		}
		free(x->symbols);
		free(x->deviations);
	}

	mpfr_clear(&(x->centre));
	mpfr_clear(&(x->radius));
}
