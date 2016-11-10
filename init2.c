/*
 * init2.c
 *
 *  Created on: 22 Sep 2016
 *      Author: jt273
 */

#include "mpfa.h"

void mpfa_init2 (mpfa_ptr x, mp_prec_t p) {
	mpfr_init2(&(x->centre), p);
	mpfr_init2(&(x->radius), p);
}
