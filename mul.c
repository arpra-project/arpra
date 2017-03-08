/*
 * mul.c
 *
 *  Created on: 11 Nov 2016
 *      Author: jt273
 */

#include "mpfa.h"
#include <malloc.h>
#include <assert.h>

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

//#define MPFA_TIGHT_MUL

void mpfa_mul (mpfa_ptr z, mpfa_srcptr x, mpfa_srcptr y) {
	unsigned xTerm, yTerm, zTerm;
	int inexact;
	mpfr_t u, temp, error, delta;
	mpfr_prec_t prec;
	mpfa_t zNew;

	prec = mpfr_get_prec(&(z->centre));
	mpfr_init2(u, prec);
	mpfr_init2(temp, prec);
	mpfr_init2(error, prec);
	mpfr_init2(delta, prec);
	mpfa_init2(zNew, prec);
	mpfr_set_si(&(zNew->radius), 0, MPFR_RNDN);

	assert(mpfr_set_si(u, -prec, MPFR_RNDN) == 0); // fails if emax <= log2(prec)
	assert(mpfr_exp2(u, u, MPFR_RNDN) == 0); // fails if emin > 1-prec

	inexact = mpfr_mul(&(zNew->centre), &(x->centre), &(y->centre), MPFR_RNDN);
	if (inexact) {
		mpfr_mul(delta, u, &(zNew->centre), MPFR_RNDA);
		mpfr_abs(delta, delta, MPFR_RNDN);
	}

	zNew->nTerms = x->nTerms + y->nTerms + 1;
	zNew->symbols = malloc(zNew->nTerms * sizeof(unsigned));
	zNew->deviations = malloc(zNew->nTerms * sizeof(mpfr_t));

	for (xTerm = 0, yTerm = 0, zTerm = 0; zTerm < (zNew->nTerms - 1); zTerm++) {
		if (yTerm == y->nTerms) {
			for (; xTerm < x->nTerms; xTerm++, zTerm++) {
				zNew->symbols[zTerm] = x->symbols[xTerm];
				mpfr_init2(&(zNew->deviations[zTerm]), prec);

				inexact = mpfr_mul(&(zNew->deviations[zTerm]), &(y->centre), &(x->deviations[xTerm]), MPFR_RNDN);
				if (inexact) {
					mpfr_mul(error, u, &(zNew->deviations[zTerm]), MPFR_RNDA);
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

				inexact = mpfr_mul(&(zNew->deviations[zTerm]), &(x->centre), &(y->deviations[yTerm]), MPFR_RNDN);
				if (inexact) {
					mpfr_mul(error, u, &(zNew->deviations[zTerm]), MPFR_RNDA);
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

			inexact = mpfr_mul(&(zNew->deviations[zTerm]), &(y->centre), &(x->deviations[xTerm]), MPFR_RNDN);
			if (inexact) {
				mpfr_mul(error, u, &(zNew->deviations[zTerm]), MPFR_RNDA);
				mpfr_abs(error, error, MPFR_RNDN);
				mpfr_add(delta, delta, error, MPFR_RNDU);
			}

			xTerm++;
		}
		else if (y->symbols[yTerm] < x->symbols[xTerm]) {
			zNew->symbols[zTerm] = y->symbols[yTerm];
			mpfr_init2(&(zNew->deviations[zTerm]), prec);

			inexact = mpfr_mul(&(zNew->deviations[zTerm]), &(x->centre), &(y->deviations[yTerm]), MPFR_RNDN);
			if (inexact) {
				mpfr_mul(error, u, &(zNew->deviations[zTerm]), MPFR_RNDA);
				mpfr_abs(error, error, MPFR_RNDN);
				mpfr_add(delta, delta, error, MPFR_RNDU);
			}

			yTerm++;
		}
		else {
			zNew->symbols[zTerm] = x->symbols[xTerm];
			mpfr_init2(&(zNew->deviations[zTerm]), prec);

			inexact = mpfr_mul(&(zNew->deviations[zTerm]), &(y->centre), &(x->deviations[xTerm]), MPFR_RNDN);
			if (inexact) {
				mpfr_mul(error, u, &(zNew->deviations[zTerm]), MPFR_RNDA);
				mpfr_abs(error, error, MPFR_RNDN);
				mpfr_add(delta, delta, error, MPFR_RNDU);
			}

			inexact = mpfr_mul(temp, &(x->centre), &(y->deviations[yTerm]), MPFR_RNDN);
			if (inexact) {
				mpfr_mul(error, u, temp, MPFR_RNDA);
				mpfr_abs(error, error, MPFR_RNDN);
				mpfr_add(delta, delta, error, MPFR_RNDU);
			}

			inexact = mpfr_add(&(zNew->deviations[zTerm]), &(zNew->deviations[zTerm]), temp, MPFR_RNDN);
			if (inexact) {
				mpfr_mul(error, u, &(zNew->deviations[zTerm]), MPFR_RNDA);
				mpfr_abs(error, error, MPFR_RNDN);
				mpfr_add(delta, delta, error, MPFR_RNDU);
			}

			xTerm++;
			yTerm++;
		}

		mpfr_abs(temp, &(zNew->deviations[zTerm]), MPFR_RNDN);
		mpfr_add(&(zNew->radius), &(zNew->radius), temp, MPFR_RNDU);
	}

#ifdef MPFA_TIGHT_MUL
	unsigned xNext, yNext;
	mpfr_t xiyiPos, xiyiNeg;

	mpfr_init2(xiyiPos, prec);
	mpfr_set_si(xiyiPos, 0, MPFR_RNDN);
	mpfr_init2(xiyiNeg, prec);
	mpfr_set_si(xiyiNeg, 0, MPFR_RNDN);

	xTerm = 0;
	yTerm = 0;
	while (1) {
		if (yTerm == y->nTerms) break;
		if (xTerm == x->nTerms) break;

		if (x->symbols[xTerm] < y->symbols[yTerm]) {
			for (yNext = yTerm; yNext < y->nTerms; yNext++) {
				// x has symbol i, and y has symbol j, so delta += abs(xi * yj)
				mpfr_mul(error, &(x->deviations[xTerm]), &(y->deviations[yNext]), MPFR_RNDA);
				mpfr_abs(error, error, MPFR_RNDN);
				mpfr_add(delta, delta, error, MPFR_RNDU);
			}
			xTerm++;
		}
		else if (y->symbols[yTerm] < x->symbols[xTerm]) {
			for (xNext = xTerm; xNext < x->nTerms; xNext++) {
				// y has symbol i, and x has symbol j, so delta += abs(yi * xj)
				mpfr_mul(error, &(y->deviations[yTerm]), &(x->deviations[xNext]), MPFR_RNDA);
				mpfr_abs(error, error, MPFR_RNDN);
				mpfr_add(delta, delta, error, MPFR_RNDU);
			}
			yTerm++;
		}
		else {
			// both x and y have symbol i, so delta += abs(xi * yi)
			mpfr_mul(error, &(x->deviations[xTerm]), &(y->deviations[yTerm]), MPFR_RNDA);
			if (mpfr_sgn(error) >= 0) {
				mpfr_add(xiyiPos, xiyiPos, error, MPFR_RNDU);
			}
			else {
				mpfr_sub(xiyiNeg, xiyiNeg, error, MPFR_RNDU);
			}

			xNext = xTerm + 1;
			yNext = yTerm + 1;
			while (1) {
				if (yNext == y->nTerms) {
					for (; xNext < x->nTerms; xNext++) {
						// both x and y have symbol i, but only x has symbol j, so delta += abs(xj * yi)
						mpfr_mul(error, &(y->deviations[yTerm]), &(x->deviations[xNext]), MPFR_RNDA);
						mpfr_abs(error, error, MPFR_RNDN);
						mpfr_add(delta, delta, error, MPFR_RNDU);
					}
					break;
				}
				if (xNext == x->nTerms) {
					for (; yNext < y->nTerms; yNext++) {
						// both x and y have symbol i, but only y has symbol j, so delta += abs(xi * yj)
						mpfr_mul(error, &(x->deviations[xTerm]), &(y->deviations[yNext]), MPFR_RNDA);
						mpfr_abs(error, error, MPFR_RNDN);
						mpfr_add(delta, delta, error, MPFR_RNDU);
					}
					break;
				}

				if (x->symbols[xNext] < y->symbols[yNext]) {
					// both x and y have symbol i, but only x has symbol j, so delta += abs(xj * yi)
					mpfr_mul(error, &(y->deviations[yTerm]), &(x->deviations[xNext]), MPFR_RNDA);
					mpfr_abs(error, error, MPFR_RNDN);
					mpfr_add(delta, delta, error, MPFR_RNDU);

					xNext++;
				}
				else if (y->symbols[yNext] < x->symbols[xNext]) {
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
			xTerm++;
			yTerm++;
		}
	}

	mpfr_max(error, xiyiPos, xiyiNeg, MPFR_RNDN);
	mpfr_add(delta, delta, error, MPFR_RNDU);

	mpfr_clear(xiyiPos);
	mpfr_clear(xiyiNeg);
#else
	mpfr_mul(error, &(x->radius), &(y->radius), MPFR_RNDU);
	mpfr_add(delta, delta, error, MPFR_RNDU);
#endif

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
	mpfr_clear(delta);
	mpfa_set(z, zNew);
	mpfa_clear(zNew);
}
