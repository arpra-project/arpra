/*
 * init_set.c
 *
 *  Created on: 21 Oct 2016
 *      Author: jt273
 */

#include "mpfa.h"

void mpfa_init_set (mpfa_ptr x, mpfa_srcptr newVal) {
	mpfr_init(&(x->centre));
	mpfr_set(&(x->centre), &(newVal->centre), MPFR_RNDN);
	mpfr_init(&(x->radius));
	mpfr_set(&(x->radius), &(newVal->radius), MPFR_RNDU);
}
