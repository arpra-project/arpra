/*
 * init_set_d.c
 *
 *  Created on: 21 Oct 2016
 *      Author: jt273
 */

#include "mpfa.h"

int mpfa_init_set_d (mpfa_ptr x, const double newVal) {
	int inexact = 0;

	mpfr_init(&(x->centre));
	inexact += mpfr_set_d(&(x->centre), newVal, MPFR_RNDN);
	mpfr_init(&(x->radius));
	inexact += mpfr_set_d(&(x->radius), 0.0, MPFR_RNDU);

	return inexact;
}
