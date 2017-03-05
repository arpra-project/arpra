/*
 * affine_2.c
 *
 *  Created on: 2 Nov 2016
 *      Author: jt273
 */

#include "mpfa.h"
#include <malloc.h>
#include <assert.h>

void mpfa_affine_2 (mpfa_ptr z, mpfa_srcptr x, mpfa_srcptr y, mpfr_srcptr alpha, mpfr_srcptr beta, mpfr_ptr gamma, mpfr_ptr delta) {
	unsigned xTerm, yTerm, zTerm;
	int inexact;
	mpfr_t u, temp, error;
	mpfr_prec_t prec;

	prec = mpfr_get_prec(&(z->centre));
	mpfr_init2(u, prec);
	mpfr_init2(temp, prec);
	mpfr_init2(error, prec);
	mpfr_set_si(&(z->radius), 0, MPFR_RNDN);

	assert(mpfr_set_si(u, -prec, MPFR_RNDN) == 0); // fails if emax <= log2(prec)
	assert(mpfr_exp2(u, u, MPFR_RNDN) == 0); // fails if emin > 1-prec

	inexact = mpfr_mul(temp, alpha, &(x->centre), MPFR_RNDN);
	if (inexact) {
		mpfr_mul(error, u, temp, MPFR_RNDU);
		mpfr_add(delta, delta, error, MPFR_RNDU);
	}

	inexact = mpfr_add(&(z->centre), gamma, temp, MPFR_RNDN);
	if (inexact) {
		mpfr_mul(error, u, &(z->centre), MPFR_RNDU);
		mpfr_add(delta, delta, error, MPFR_RNDU);
	}

	inexact = mpfr_mul(temp, beta, &(y->centre), MPFR_RNDN);
	if (inexact) {
		mpfr_mul(error, u, temp, MPFR_RNDU);
		mpfr_add(delta, delta, error, MPFR_RNDU);
	}

	inexact = mpfr_add(&(z->centre), &(z->centre), temp, MPFR_RNDN);
	if (inexact) {
		mpfr_mul(error, u, &(z->centre), MPFR_RNDU);
		mpfr_add(delta, delta, error, MPFR_RNDU);
	}

	if (z->nTerms == 0) {
		z->nTerms = x->nTerms + y->nTerms + 1;
		z->symbols = malloc(z->nTerms * sizeof(unsigned));
		z->deviations = malloc(z->nTerms * sizeof(mpfr_t));
	}
	else {
		for (zTerm = 0; zTerm < z->nTerms; zTerm++) {
			mpfr_clear(&(z->deviations[zTerm]));
		}
		z->nTerms = x->nTerms + y->nTerms + 1;
		z->symbols = realloc(z->symbols, z->nTerms * sizeof(unsigned));
		z->deviations = realloc(z->symbols, z->nTerms * sizeof(mpfr_t));
	}

	for (xTerm = 0, yTerm = 0, zTerm = 0; ((xTerm < x->nTerms) || (yTerm < y->nTerms)); zTerm++) {
		if ((yTerm == y->nTerms) || (x->symbols[xTerm] < y->symbols[yTerm])) {
			z->symbols[zTerm] = x->symbols[xTerm];
			mpfr_init(&(z->deviations[zTerm]));

			inexact = mpfr_mul(&(z->deviations[zTerm]), alpha, &(x->deviations[xTerm]), MPFR_RNDN);
			if (inexact) {
				mpfr_mul(error, u, &(z->deviations[zTerm]), MPFR_RNDU);
				mpfr_add(delta, delta, error, MPFR_RNDU);
			}

			xTerm++;
		}
		else if ((xTerm == x->nTerms) || (y->symbols[yTerm] < x->symbols[xTerm])) {
			z->symbols[zTerm] = y->symbols[yTerm];
			mpfr_init(&(z->deviations[zTerm]));

			inexact = mpfr_mul(&(z->deviations[zTerm]), beta, &(y->deviations[yTerm]), MPFR_RNDN);
			if (inexact) {
				mpfr_mul(error, u, &(z->deviations[zTerm]), MPFR_RNDU);
				mpfr_add(delta, delta, error, MPFR_RNDU);
			}

			yTerm++;
		}
		else {
			z->symbols[zTerm] = x->symbols[xTerm];
			mpfr_init(&(z->deviations[zTerm]));

			inexact = mpfr_mul(&(z->deviations[zTerm]), alpha, &(x->deviations[xTerm]), MPFR_RNDN);
			if (inexact) {
				mpfr_mul(error, u, &(z->deviations[zTerm]), MPFR_RNDU);
				mpfr_add(delta, delta, error, MPFR_RNDU);
			}

			inexact = mpfr_mul(temp, beta, &(y->deviations[yTerm]), MPFR_RNDN);
			if (inexact) {
				mpfr_mul(error, u, temp, MPFR_RNDU);
				mpfr_add(delta, delta, error, MPFR_RNDU);
			}

			inexact = mpfr_add(&(z->deviations[zTerm]), &(z->deviations[zTerm]), temp, MPFR_RNDN);
			if (inexact) {
				mpfr_mul(error, u, &(z->deviations[zTerm]), MPFR_RNDU);
				mpfr_add(delta, delta, error, MPFR_RNDU);
			}

			xTerm++;
			yTerm++;
		}
		mpfr_abs(temp, &(z->deviations[zTerm]), MPFR_RNDN);
		mpfr_add(&(z->radius), &(z->radius), temp, MPFR_RNDU);
	}

	z->nTerms = zTerm + 1;
	z->symbols = realloc(z->symbols, z->nTerms * sizeof(unsigned));
	z->symbols[zTerm] = mpfa_next_sym();
	z->deviations = realloc(z->deviations, z->nTerms * sizeof(mpfr_t));
	mpfr_init_set(&(z->deviations[zTerm]), delta, MPFR_RNDU);
	mpfr_add(&(z->radius), &(z->radius), delta, MPFR_RNDU);

	mpfr_clear(u);
	mpfr_clear(temp);
	mpfr_clear(error);
}
