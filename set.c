/*
 * set.c
 *
 *  Created on: 21 Oct 2016
 *      Author: jt273
 */

#include "mpfa.h"

int mpfa_set (mpfa_ptr x, mpfa_srcptr newVal) {
	int inexact = 0;

	inexact += mpfr_set(&(x->centre), &(newVal->centre), MPFR_RNDN);
	inexact += mpfr_set(&(x->radius), &(newVal->radius), MPFR_RNDU);

	return inexact;
}
