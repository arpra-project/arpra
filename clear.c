/*
 * clear.c
 *
 *  Created on: 20 Sep 2016
 *      Author: james
 */

#include "mpfa.h"

void mpfa_clear (mpfa_ptr x) {
	mpfr_clear (&(x->centre));
	mpfr_clear (&(x->radius));
}
