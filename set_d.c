/*
 * set_d.c
 *
 *  Created on: 21 Oct 2016
 *      Author: jt273
 */

#include "mpfa.h"
#include <malloc.h>

void mpfa_set_d (mpfa_ptr z, const double x) {
	if (z->nTerms > 0) {
		unsigned zTerm;
		for (zTerm = 0; zTerm < z->nTerms; zTerm++) {
			mpfr_clear(&(z->deviations[zTerm]));
		}
		z->nTerms = 0;
		free(z->symbols);
		free(z->deviations);
		mpfr_set_si(&(z->radius), 0, MPFR_RNDU);
	}
	mpfr_set_d(&(z->centre), x, MPFR_RNDN);
}
