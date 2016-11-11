/*
 * init_set_d.c
 *
 *  Created on: 21 Oct 2016
 *      Author: jt273
 */

#include "mpfa.h"

void mpfa_init_set_d (mpfa_ptr x, const double newVal) {
	mpfr_init(&(x->centre));
	mpfr_set_d(&(x->centre), newVal, MPFR_RNDN);
	mpfr_init(&(x->radius));
	mpfr_set_d(&(x->radius), 0.0, MPFR_RNDU);
}
