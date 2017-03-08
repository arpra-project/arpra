/*
 * affine_2.c
 *
 *  Created on: 2 Nov 2016
 *      Author: jt273
 */

#include "mpfa.h"
#include <malloc.h>
#include <assert.h>

void mpfa_affine_2 (mpfa_ptr z, mpfa_srcptr x, mpfa_srcptr y, mpfr_ptr alpha, mpfr_ptr beta, mpfr_ptr gamma, mpfr_ptr delta) {
	unsigned xTerm, yTerm, zTerm;
	int inexact;
	mpfr_t u, temp, error;
	mpfr_prec_t prec;
	mpfa_t zNew;

	prec = mpfr_get_prec(&(z->centre));
	mpfr_init2(u, prec);
	mpfr_init2(temp, prec);
	mpfr_init2(error, prec);
	mpfa_init2(zNew, prec);
	mpfr_set_si(&(zNew->radius), 0, MPFR_RNDN);

	assert(mpfr_set_si(u, -prec, MPFR_RNDN) == 0); // fails if emax <= log2(prec)
	assert(mpfr_exp2(u, u, MPFR_RNDN) == 0); // fails if emin > 1-prec

	inexact = mpfr_mul(temp, alpha, &(x->centre), MPFR_RNDN);
	if (inexact) {
		mpfr_mul(error, u, temp, MPFR_RNDU);
		mpfr_abs(error, error, MPFR_RNDN);
		mpfr_add(delta, delta, error, MPFR_RNDU);
	}

	inexact = mpfr_add(&(zNew->centre), gamma, temp, MPFR_RNDN);
	if (inexact) {
		mpfr_mul(error, u, &(zNew->centre), MPFR_RNDU);
		mpfr_abs(error, error, MPFR_RNDN);
		mpfr_add(delta, delta, error, MPFR_RNDU);
	}

	inexact = mpfr_mul(temp, beta, &(y->centre), MPFR_RNDN);
	if (inexact) {
		mpfr_mul(error, u, temp, MPFR_RNDU);
		mpfr_abs(error, error, MPFR_RNDN);
		mpfr_add(delta, delta, error, MPFR_RNDU);
	}

	inexact = mpfr_add(&(zNew->centre), &(zNew->centre), temp, MPFR_RNDN);
	if (inexact) {
		mpfr_mul(error, u, &(zNew->centre), MPFR_RNDU);
		mpfr_abs(error, error, MPFR_RNDN);
		mpfr_add(delta, delta, error, MPFR_RNDU);
	}

	zNew->nTerms = x->nTerms + y->nTerms + 1;
	zNew->symbols = malloc(zNew->nTerms * sizeof(unsigned));
	zNew->deviations = malloc(zNew->nTerms * sizeof(mpfr_t));

	for (xTerm = 0, yTerm = 0, zTerm = 0; zTerm < (zNew->nTerms - 1); zTerm++) {
		if (yTerm == y->nTerms) {
			for (; xTerm < x->nTerms; xTerm++, zTerm++) {
				zNew->symbols[zTerm] = x->symbols[xTerm];
				mpfr_init2(&(zNew->deviations[zTerm]), prec);

				inexact = mpfr_mul(&(zNew->deviations[zTerm]), alpha, &(x->deviations[xTerm]), MPFR_RNDN);
				if (inexact) {
					mpfr_mul(error, u, &(zNew->deviations[zTerm]), MPFR_RNDU);
					mpfr_abs(error, error, MPFR_RNDN);
					mpfr_add(delta, delta, error, MPFR_RNDU);
				}

				mpfr_abs(temp, &(zNew->deviations[zTerm]), MPFR_RNDN);
				mpfr_add(&(zNew->radius), &(zNew->radius), temp, MPFR_RNDU);
			}
			break;
		}
		if (xTerm == x->nTerms) {
			for (; yTerm < y->nTerms; yTerm++, zTerm++) {
				zNew->symbols[zTerm] = y->symbols[yTerm];
				mpfr_init2(&(zNew->deviations[zTerm]), prec);

				inexact = mpfr_mul(&(zNew->deviations[zTerm]), beta, &(y->deviations[yTerm]), MPFR_RNDN);
				if (inexact) {
					mpfr_mul(error, u, &(zNew->deviations[zTerm]), MPFR_RNDU);
					mpfr_abs(error, error, MPFR_RNDN);
					mpfr_add(delta, delta, error, MPFR_RNDU);
				}

				mpfr_abs(temp, &(zNew->deviations[zTerm]), MPFR_RNDN);
				mpfr_add(&(zNew->radius), &(zNew->radius), temp, MPFR_RNDU);
			}
			break;
		}

		if (x->symbols[xTerm] < y->symbols[yTerm]) {
			zNew->symbols[zTerm] = x->symbols[xTerm];
			mpfr_init2(&(zNew->deviations[zTerm]), prec);

			inexact = mpfr_mul(&(zNew->deviations[zTerm]), alpha, &(x->deviations[xTerm]), MPFR_RNDN);
			if (inexact) {
				mpfr_mul(error, u, &(zNew->deviations[zTerm]), MPFR_RNDU);
				mpfr_abs(error, error, MPFR_RNDN);
				mpfr_add(delta, delta, error, MPFR_RNDU);
			}

			xTerm++;
		}
		else if (y->symbols[yTerm] < x->symbols[xTerm]) {
			zNew->symbols[zTerm] = y->symbols[yTerm];
			mpfr_init2(&(zNew->deviations[zTerm]), prec);

			inexact = mpfr_mul(&(zNew->deviations[zTerm]), beta, &(y->deviations[yTerm]), MPFR_RNDN);
			if (inexact) {
				mpfr_mul(error, u, &(zNew->deviations[zTerm]), MPFR_RNDU);
				mpfr_abs(error, error, MPFR_RNDN);
				mpfr_add(delta, delta, error, MPFR_RNDU);
			}

			yTerm++;
		}
		else {
			zNew->symbols[zTerm] = x->symbols[xTerm];
			mpfr_init2(&(zNew->deviations[zTerm]), prec);

			inexact = mpfr_mul(&(zNew->deviations[zTerm]), alpha, &(x->deviations[xTerm]), MPFR_RNDN);
			if (inexact) {
				mpfr_mul(error, u, &(zNew->deviations[zTerm]), MPFR_RNDU);
				mpfr_abs(error, error, MPFR_RNDN);
				mpfr_add(delta, delta, error, MPFR_RNDU);
			}

			inexact = mpfr_mul(temp, beta, &(y->deviations[yTerm]), MPFR_RNDN);
			if (inexact) {
				mpfr_mul(error, u, temp, MPFR_RNDU);
				mpfr_abs(error, error, MPFR_RNDN);
				mpfr_add(delta, delta, error, MPFR_RNDU);
			}

			inexact = mpfr_add(&(zNew->deviations[zTerm]), &(zNew->deviations[zTerm]), temp, MPFR_RNDN);
			if (inexact) {
				mpfr_mul(error, u, &(zNew->deviations[zTerm]), MPFR_RNDU);
				mpfr_abs(error, error, MPFR_RNDN);
				mpfr_add(delta, delta, error, MPFR_RNDU);
			}

			xTerm++;
			yTerm++;
		}

		mpfr_abs(temp, &(zNew->deviations[zTerm]), MPFR_RNDN);
		mpfr_add(&(zNew->radius), &(zNew->radius), temp, MPFR_RNDU);
	}

	zNew->nTerms = zTerm + 1;
	zNew->symbols = realloc(zNew->symbols, zNew->nTerms * sizeof(unsigned));
	zNew->symbols[zTerm] = mpfa_next_sym();
	zNew->deviations = realloc(zNew->deviations, zNew->nTerms * sizeof(mpfr_t));
	mpfr_init2(&(zNew->deviations[zTerm]), prec);
	mpfr_set(&(zNew->deviations[zTerm]), delta, MPFR_RNDN);
	mpfr_add(&(zNew->radius), &(zNew->radius), delta, MPFR_RNDU);

	mpfr_clear(u);
	mpfr_clear(temp);
	mpfr_clear(error);
	mpfa_set(z, zNew);
	mpfa_clear(zNew);
}
