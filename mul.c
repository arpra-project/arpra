/*
 * mul.c
 *
 *  Created on: 11 Nov 2016
 *      Author: jt273
 */

#include "mpfa.h"
#include <malloc.h>
#include <assert.h>

#define MPFA_TIGHT_MUL

/*
 * If MPFA_TIGHT_MUL is defined, then the linear approximation to the quadratic term of the
 * affine product is defined the same as in (26) of:
 *
 * S. M. Rump and M. Kashiwagi, Implementation and improvements of affine arithmetic,
 * Nonlinear Theory an Its Applications, IEICE, vol. 6, no. 3, pp. 341-359, 2015.
 *
 * Otherwise it is trivially defined as the product of the radii of x and y:
 *
 * \sum^{n}_{i=1} x_{i} \sum^{n}_{i=1} y_{i}
 */

void mpfa_mul (mpfa_ptr z, mpfa_srcptr x, mpfa_srcptr y) {
	unsigned xTerm, yTerm, zTerm;
	int inexact;
	mpfr_t u, temp, error, delta;
	mpfr_prec_t prec;

	prec = mpfr_get_prec(&(z->centre));
	mpfr_init2(u, prec);
	mpfr_init2(temp, prec);
	mpfr_init2(error, prec);
	mpfr_init2(delta, prec);
	mpfr_set_si(&(z->radius), 0, MPFR_RNDN);

#ifdef MPFA_TIGHT_MUL
	unsigned xNext, yNext;
	mpfr_t xiyiPos, xiyiNeg;

	mpfr_init2(xiyiPos, prec);
	mpfr_set_si(xiyiPos, 0, MPFR_RNDN);
	mpfr_init2(xiyiNeg, prec);
	mpfr_set_si(xiyiNeg, 0, MPFR_RNDN);
#endif

	assert(!mpfr_set_si(u, -prec, MPFR_RNDN)); // fails if emax <= log2(prec)
	assert(!mpfr_exp2(u, u, MPFR_RNDN)); // fails if emin > 1-prec

	inexact = mpfr_mul(&(z->centre), &(x->centre), &(y->centre), MPFR_RNDN);
	if (inexact) {
		mpfr_mul(delta, u, &(z->centre), MPFR_RNDU);
	}

	for (zTerm = 0; zTerm < z->nTerms; zTerm++) {
		mpfr_clear(&(z->deviations[zTerm]));
	}
	z->nTerms = x->nTerms + y->nTerms + 1;
	z->symbols = realloc(z->symbols, z->nTerms * sizeof(unsigned));
	z->deviations = realloc(z->symbols, z->nTerms * sizeof(mpfr_t));

	for (xTerm = 0, yTerm = 0, zTerm = 0; ((xTerm < x->nTerms) || (yTerm < y->nTerms)); zTerm++) {
		if ((yTerm == y->nTerms) || (x->symbols[xTerm] < y->symbols[yTerm])) {
			z->symbols[zTerm] = x->symbols[xTerm];
			mpfr_init2(&(z->deviations[zTerm]), prec);

			inexact = mpfr_mul(&(z->deviations[zTerm]), &(y->centre), &(x->deviations[xTerm]), MPFR_RNDN);
			if (inexact) {
				mpfr_mul(error, u, &(z->deviations[zTerm]), MPFR_RNDU);
				mpfr_add(delta, delta, error, MPFR_RNDU);
			}

#ifdef MPFA_TIGHT_MUL
			for (yNext = yTerm; yNext < y->nTerms; yNext++) {
				// x has symbol i, and y has symbol j, so delta += abs(xi * yj)
				mpfr_mul(error, &(x->deviations[xTerm]), &(y->deviations[yNext]), MPFR_RNDA);
				mpfr_abs(error, error, MPFR_RNDN);
				mpfr_add(delta, delta, error, MPFR_RNDU);
			}
#endif

			xTerm++;
		}
		else if ((xTerm == x->nTerms) || (y->symbols[yTerm] < x->symbols[xTerm])) {
			z->symbols[zTerm] = y->symbols[yTerm];
			mpfr_init2(&(z->deviations[zTerm]), prec);

			inexact = mpfr_mul(&(z->deviations[zTerm]), &(x->centre), &(y->deviations[yTerm]), MPFR_RNDN);
			if (inexact) {
				mpfr_mul(error, u, &(z->deviations[zTerm]), MPFR_RNDU);
				mpfr_add(delta, delta, error, MPFR_RNDU);
			}

#ifdef MPFA_TIGHT_MUL
			for (xNext = xTerm; xNext < x->nTerms; xNext++) {
				// y has symbol i, and x has symbol j, so delta += abs(yi * xj)
				mpfr_mul(error, &(y->deviations[yTerm]), &(x->deviations[xNext]), MPFR_RNDA);
				mpfr_abs(error, error, MPFR_RNDN);
				mpfr_add(delta, delta, error, MPFR_RNDU);
			}
#endif

			yTerm++;
		}
		else {
			z->symbols[zTerm] = x->symbols[xTerm];
			mpfr_init2(&(z->deviations[zTerm]), prec);

			inexact = mpfr_mul(&(z->deviations[zTerm]), &(y->centre), &(x->deviations[xTerm]), MPFR_RNDN);
			if (inexact) {
				mpfr_mul(error, u, &(z->deviations[zTerm]), MPFR_RNDU);
				mpfr_add(delta, delta, error, MPFR_RNDU);
			}

			inexact = mpfr_mul(temp, &(x->centre), &(y->deviations[yTerm]), MPFR_RNDN);
			if (inexact) {
				mpfr_mul(error, u, temp, MPFR_RNDU);
				mpfr_add(delta, delta, error, MPFR_RNDU);
			}

			inexact = mpfr_add(&(z->deviations[zTerm]), &(z->deviations[zTerm]), temp, MPFR_RNDN);
			if (inexact) {
				mpfr_mul(error, u, &(z->deviations[zTerm]), MPFR_RNDU);
				mpfr_add(delta, delta, error, MPFR_RNDU);
			}

#ifdef MPFA_TIGHT_MUL
			// both x and y have symbol i, so delta += abs(xi * yi)
			mpfr_mul(error, &(x->deviations[xTerm]), &(y->deviations[yTerm]), MPFR_RNDA);
			if (mpfr_sgn(error) >= 0) {
				mpfr_add(xiyiPos, xiyiPos, error, MPFR_RNDU);
			}
			else {
				mpfr_sub(xiyiNeg, xiyiNeg, error, MPFR_RNDU);
			}

			for (xNext = (xTerm + 1), yNext = (yTerm + 1); ((xNext < x->nTerms) || (yNext < y->nTerms)); ) {
				if ((yNext == y->nTerms) || (x->symbols[xNext] < y->symbols[yNext])) {
					// both x and y have symbol i, but only x has symbol j, so delta += abs(xj * yi)
					mpfr_mul(error, &(y->deviations[yTerm]), &(x->deviations[xNext]), MPFR_RNDA);
					mpfr_abs(error, error, MPFR_RNDN);
					mpfr_add(delta, delta, error, MPFR_RNDU);

					xNext++;
				}
				else if ((xNext == x->nTerms) || (y->symbols[yNext] < x->symbols[xNext])) {
					// both x and y have symbol i, but only y has symbol j, so delta += abs(xi * yj)
					mpfr_mul(error, &(x->deviations[xTerm]), &(y->deviations[yNext]), MPFR_RNDA);
					mpfr_abs(error, error, MPFR_RNDN);
					mpfr_add(delta, delta, error, MPFR_RNDU);

					yNext++;
				}
				else {
					// both x and y have symbols i and j, so delta += abs(xi * yj + xj * yi)
					mpfr_mul(error, &(x->deviations[xTerm]), &(y->deviations[yNext]), MPFR_RNDU);
					mpfr_mul(temp, &(x->deviations[xTerm]), &(y->deviations[yNext]), MPFR_RNDU);
					mpfr_add(error, error, temp, MPFR_RNDU);
					if (mpfr_sgn(error) < 0) {
						mpfr_mul(error, &(x->deviations[xTerm]), &(y->deviations[yNext]), MPFR_RNDD);
						mpfr_mul(temp, &(x->deviations[xTerm]), &(y->deviations[yNext]), MPFR_RNDD);
						mpfr_add(error, error, temp, MPFR_RNDD);
					}
					mpfr_abs(error, error, MPFR_RNDN);
					mpfr_add(delta, delta, error, MPFR_RNDU);

					xNext++;
					yNext++;
				}
			}
#endif

			xTerm++;
			yTerm++;
		}
		mpfr_abs(temp, &(z->deviations[zTerm]), MPFR_RNDN);
		mpfr_add(&(z->radius), &(z->radius), temp, MPFR_RNDU);
	}

#ifdef MPFA_TIGHT_MUL
	mpfr_max(error, xiyiPos, xiyiNeg, MPFR_RNDN);
	mpfr_add(delta, delta, error, MPFR_RNDU);

	mpfr_clear(xiyiPos);
	mpfr_clear(xiyiNeg);
#else
	mpfr_mul(error, &(x->radius), &(y->radius), MPFR_RNDU);
	mpfr_add(delta, delta, error, MPFR_RNDU);
#endif

	z->nTerms = zTerm + 1;
	z->symbols = realloc(z->symbols, z->nTerms * sizeof(unsigned));
	z->symbols[zTerm] = mpfa_next_sym();
	z->deviations = realloc(z->symbols, z->nTerms * sizeof(mpfr_t));
	mpfr_init_set(&(z->deviations[zTerm]), delta, MPFR_RNDU);
	mpfr_add(&(z->radius), &(z->radius), delta, MPFR_RNDU);

	mpfr_clear(u);
	mpfr_clear(temp);
	mpfr_clear(error);
	mpfr_clear(delta);
}
