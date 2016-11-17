/*
 * div.c
 *
 *  Created on: 16 Nov 2016
 *      Author: jt273
 */

#include "mpfa.h"
#include <malloc.h>

void mpfa_div (mpfa_ptr z, mpfa_srcptr x, mpfa_srcptr y) {
	unsigned xTerm, yTerm, zTerm;
	mpfa_t z_new;
	mpfr_t temp, delta;

	mpfr_init(temp);
	mpfr_init_set_d(delta, 0.0, MPFR_RNDN);
	z_new->nTerms = x->nTerms + y->nTerms + 1;
	z_new->symbols = malloc(z_new->nTerms * sizeof(unsigned));
	z_new->deviations = malloc(z_new->nTerms * sizeof(mpfr_t));
	mpfr_init(&(z_new->centre));


	// change to divide	//mpfa_affine_term_1 (&(z_new->centre), &(x->centre), &(y->centre), NULL, delta);


	mpfr_init_set_d(&(z_new->radius), 0.0, MPFR_RNDN);







	mpfa_clear(z);
	*z = *z_new;
	mpfr_clear(temp);
	mpfr_clear(delta);
}
