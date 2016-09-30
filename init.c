/*
 * init.c
 *
 *  Created on: 21 Jul 2016
 *      Author: jt273
 */

#include "mpfa.h"

void mpfa_init (mpfa_ptr x) {
	mpfr_init (&(x->centre));
	mpfr_init (&(x->radius));
}
