/*
 * mul.c
 *
 *  Created on: 11 Nov 2016
 *      Author: jt273
 */

#include "mpfa.h"
#include <malloc.h>
#include <assert.h>

void mpfa_mul (mpfa_ptr z, mpfa_srcptr x, mpfa_srcptr y) {
	unsigned xTerm, yTerm, zTerm;
	int inexact;
	mpfr_prec_t prec;
	mpfr_t u, temp, error, delta;
	mpfa_t z_new;

	prec = mpfr_get_prec(&(z->centre));
	mpfr_init2(u, prec);
	mpfr_init2(temp, prec);
	mpfr_init2(error, prec);
	mpfr_init2(delta, prec);
	mpfr_init2(&(z_new->centre), prec);
	mpfr_init2(&(z_new->radius), prec);
	mpfr_set_d(&(z_new->radius), 0.0, MPFR_RNDN);

	assert(!mpfr_set_si(u, -prec, MPFR_RNDN)); // fails if emax <= log2(prec)
	assert(!mpfr_exp2(u, u, MPFR_RNDN)); // fails if emin > 1-prec

	inexact = mpfr_mul(&(z_new->centre), &(x->centre), &(y->centre), MPFR_RNDN);
	if (inexact) {
		mpfr_mul(delta, u, &(z_new->centre), MPFR_RNDU);
	}

	z_new->nTerms = x->nTerms + y->nTerms + 1;
	z_new->symbols = malloc(z_new->nTerms * sizeof(unsigned));
	z_new->deviations = malloc(z_new->nTerms * sizeof(mpfr_t));

	for (xTerm = 0, yTerm = 0, zTerm = 0; ((xTerm < x->nTerms) && (yTerm < y->nTerms)); zTerm++) {
		if ((yTerm == y->nTerms) || (x->symbols[xTerm] < y->symbols[yTerm])) {
			z_new->symbols[zTerm] = x->symbols[xTerm];
			mpfr_init2(&(z_new->deviations[zTerm]), prec);

			inexact = mpfr_mul(&(z_new->deviations[zTerm]), &(y->centre), &(x->deviations[xTerm]), MPFR_RNDN);
			if (inexact) {
				mpfr_mul(error, u, &(z_new->deviations[zTerm]), MPFR_RNDU);
				mpfr_add(delta, delta, error, MPFR_RNDU);
			}

			xTerm++;
		}
		else if ((xTerm == x->nTerms) || (y->symbols[yTerm] < x->symbols[xTerm])) {
			z_new->symbols[zTerm] = y->symbols[yTerm];
			mpfr_init2(&(z_new->deviations[zTerm]), prec);

			inexact = mpfr_mul(&(z_new->deviations[zTerm]), &(x->centre), &(y->deviations[yTerm]), MPFR_RNDN);
			if (inexact) {
				mpfr_mul(error, u, &(z_new->deviations[zTerm]), MPFR_RNDU);
				mpfr_add(delta, delta, error, MPFR_RNDU);
			}

			yTerm++;
		}
		else {
			z_new->symbols[zTerm] = x->symbols[xTerm];
			mpfr_init2(&(z_new->deviations[zTerm]), prec);

			inexact = mpfr_mul(&(z_new->deviations[zTerm]), &(y->centre), &(x->deviations[xTerm]), MPFR_RNDN);
			if (inexact) {
				mpfr_mul(error, u, &(z_new->deviations[zTerm]), MPFR_RNDU);
				mpfr_add(delta, delta, error, MPFR_RNDU);
			}

			inexact = mpfr_mul(temp, &(x->centre), &(y->deviations[yTerm]), MPFR_RNDN);
			if (inexact) {
				mpfr_mul(error, u, temp, MPFR_RNDU);
				mpfr_add(delta, delta, error, MPFR_RNDU);
			}

			inexact = mpfr_add(&(z_new->deviations[zTerm]), &(z_new->deviations[zTerm]), temp, MPFR_RNDN);
			if (inexact) {
				mpfr_mul(error, u, &(z_new->deviations[zTerm]), MPFR_RNDU);
				mpfr_add(delta, delta, error, MPFR_RNDU);
			}

			xTerm++;
			yTerm++;
		}
		mpfr_abs(temp, &(z_new->deviations[zTerm]), MPFR_RNDN);
		mpfr_add(&(z_new->radius), &(z_new->radius), temp, MPFR_RNDU);
	}

	// Trivial estimate of linearisation error
	mpfr_mul(error, &(x->radius), &(y->radius), MPFR_RNDU);
	mpfr_add(delta, delta, error, MPFR_RNDU);

	z_new->nTerms = zTerm + 1;
	z_new->symbols = realloc(z_new->symbols, z_new->nTerms * sizeof(unsigned));
	z_new->symbols[zTerm] = mpfa_next_sym();
	z_new->deviations = realloc(z_new->symbols, z_new->nTerms * sizeof(mpfr_t));
	mpfr_init_set(&(z_new->deviations[zTerm]), delta, MPFR_RNDU);
	mpfr_add(&(z_new->radius), &(z_new->radius), delta, MPFR_RNDU);

	mpfr_clear(u);
	mpfr_clear(temp);
	mpfr_clear(error);
	mpfr_clear(delta);
	mpfa_clear(z);
	*z = *z_new;
}


// Better (Rump and Kashiwagi) estimate of linearisation error

//void mpfa_mul_rump_kashiwagi (mpfa_ptr z, mpfa_srcptr x, mpfa_srcptr y) {
