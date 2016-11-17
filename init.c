/*
 * init.c
 *
 *  Created on: 21 Jul 2016
 *      Author: jt273
 */

#include "mpfa.h"
#include <malloc.h>

void mpfa_init (mpfa_ptr x) {
	x->symbols = malloc(0);
	x->deviations = malloc(0);
	mpfr_init(&(x->centre));
	mpfr_init(&(x->radius));
}
