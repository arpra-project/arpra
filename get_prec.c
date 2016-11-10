/*
 * get_prec.c
 *
 *  Created on: 7 Sep 2016
 *      Author: jt273
 */

#include "mpfa.h"

mp_prec_t mpfa_get_prec (mpfa_srcptr x) {
	return mpfr_get_prec(&(x->centre));
}
