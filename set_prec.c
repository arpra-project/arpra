/*
 * set_prec.c
 *
 *  Created on: 22 Sep 2016
 *      Author: jt273
 */

#include "mpfa.h"

void mpfa_set_prec (mpfa_ptr x, mpfr_prec_t p) {
	mpfr_set_prec(&(x->centre), p);
	mpfr_set_prec(&(x->radius), p);
}
