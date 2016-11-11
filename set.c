/*
 * set.c
 *
 *  Created on: 21 Oct 2016
 *      Author: jt273
 */

#include "mpfa.h"

void mpfa_set (mpfa_ptr x, mpfa_srcptr newVal) {
	mpfr_set(&(x->centre), &(newVal->centre), MPFR_RNDN);
	mpfr_set(&(x->radius), &(newVal->radius), MPFR_RNDU);
}
