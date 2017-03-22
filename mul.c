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
	int xHasNext, yHasNext;
	mpfr_t temp, error;
	mpfr_prec_t prec;
	mpfa_t zNew;

	prec = mpfr_get_prec(&(z->centre));
	mpfr_inits2(prec, temp, error, (mpfr_ptr) NULL);
	mpfa_init2(zNew, prec);
	mpfr_set_si(error, 0, MPFR_RNDN);
	mpfr_set_si(&(zNew->radius), 0, MPFR_RNDN);

	if (mpfr_mul(&(zNew->centre), &(x->centre), &(y->centre), MPFR_RNDN)) {
		assert(mpfr_set_si(temp, (-prec + mpfr_get_exp(&(zNew->centre))), MPFR_RNDN) == 0);
		assert(mpfr_exp2(temp, temp, MPFR_RNDN) == 0);
		mpfr_add(error, error, temp, MPFR_RNDU);
	}

	zNew->nTerms = x->nTerms + y->nTerms + 1;
	zNew->symbols = malloc(zNew->nTerms * sizeof(unsigned));
	zNew->deviations = malloc(zNew->nTerms * sizeof(mpfr_t));

	xTerm = 0; yTerm = 0; zTerm = 0;
	xHasNext = x->nTerms > 0; yHasNext = y->nTerms > 0;
	while (xHasNext || yHasNext) {
		if ((!yHasNext) || (xHasNext && (x->symbols[xTerm] < y->symbols[yTerm]))) {
			zNew->symbols[zTerm] = x->symbols[xTerm];
			mpfr_init2(&(zNew->deviations[zTerm]), prec);

			if (mpfr_mul(&(zNew->deviations[zTerm]), &(y->centre), &(x->deviations[xTerm]), MPFR_RNDN)) {
				assert(mpfr_set_si(temp, (-prec + mpfr_get_exp(&(zNew->deviations[zTerm]))), MPFR_RNDN) == 0);
				assert(mpfr_exp2(temp, temp, MPFR_RNDN) == 0);
				mpfr_add(error, error, temp, MPFR_RNDU);
			}

			xHasNext = ++xTerm < x->nTerms;
		}
		else if ((!xHasNext) || (yHasNext && (y->symbols[yTerm] < x->symbols[xTerm]))) {
			zNew->symbols[zTerm] = y->symbols[yTerm];
			mpfr_init2(&(zNew->deviations[zTerm]), prec);

			if (mpfr_mul(&(zNew->deviations[zTerm]), &(x->centre), &(y->deviations[yTerm]), MPFR_RNDN)) {
				assert(mpfr_set_si(temp, (-prec + mpfr_get_exp(&(zNew->deviations[zTerm]))), MPFR_RNDN) == 0);
				assert(mpfr_exp2(temp, temp, MPFR_RNDN) == 0);
				mpfr_add(error, error, temp, MPFR_RNDU);
			}

			yHasNext = ++yTerm < y->nTerms;
		}
		else {
			zNew->symbols[zTerm] = x->symbols[xTerm];
			mpfr_init2(&(zNew->deviations[zTerm]), prec);

			if (mpfa_term(&(zNew->deviations[zTerm]), &(x->deviations[xTerm]), &(y->deviations[yTerm]), &(y->centre), &(x->centre), NULL)) {
				assert(mpfr_set_si(temp, (-prec + mpfr_get_exp(&(zNew->deviations[zTerm]))), MPFR_RNDN) == 0);
				assert(mpfr_exp2(temp, temp, MPFR_RNDN) == 0);
				mpfr_add(error, error, temp, MPFR_RNDU);
			}

			xHasNext = ++xTerm < x->nTerms;
			yHasNext = ++yTerm < y->nTerms;
		}

		if (mpfr_zero_p(&(zNew->deviations[zTerm]))) {
			mpfr_clear(&(zNew->deviations[zTerm]));
		}
		else {
			mpfr_abs(temp, &(zNew->deviations[zTerm]), MPFR_RNDN);
			mpfr_add(&(zNew->radius), &(zNew->radius), temp, MPFR_RNDU);
			zTerm++;
		}
	}

#ifdef MPFA_TIGHT_MUL
	unsigned xNext, yNext;
	mpfr_t xiyiPos, xiyiNeg;

	mpfr_init2(xiyiPos, prec);
	mpfr_set_si(xiyiPos, 0, MPFR_RNDN);
	mpfr_init2(xiyiNeg, prec);
	mpfr_set_si(xiyiNeg, 0, MPFR_RNDN);

	xTerm = 0; yTerm = 0;
	while ((xTerm < x->nTerms) && (yTerm < y->nTerms)) {
		if (x->symbols[xTerm] < y->symbols[yTerm]) {
			for (yNext = yTerm; yNext < y->nTerms; yNext++) {
				// x has symbol i, and y has symbol j, so error += abs(xi * yj)
				mpfr_mul(temp, &(x->deviations[xTerm]), &(y->deviations[yNext]), MPFR_RNDA);
				mpfr_abs(temp, temp, MPFR_RNDN);
				mpfr_add(error, error, temp, MPFR_RNDU);
			}
			xTerm++;
		}
		else if (y->symbols[yTerm] < x->symbols[xTerm]) {
			for (xNext = xTerm; xNext < x->nTerms; xNext++) {
				// y has symbol i, and x has symbol j, so error += abs(yi * xj)
				mpfr_mul(temp, &(y->deviations[yTerm]), &(x->deviations[xNext]), MPFR_RNDA);
				mpfr_abs(temp, temp, MPFR_RNDN);
				mpfr_add(error, error, temp, MPFR_RNDU);
			}
			yTerm++;
		}
		else {
			// both x and y have symbol i, so error += abs(xi * yi)
			mpfr_mul(temp, &(x->deviations[xTerm]), &(y->deviations[yTerm]), MPFR_RNDA);
			if (mpfr_sgn(temp) >= 0) {
				mpfr_add(xiyiPos, xiyiPos, temp, MPFR_RNDU);
			}
			else {
				mpfr_sub(xiyiNeg, xiyiNeg, temp, MPFR_RNDU);
			}

			xNext = xTerm + 1; yNext = yTerm + 1;
			xHasNext = xNext < x->nTerms; yHasNext = yNext < y->nTerms;
			while (xHasNext || yHasNext) {
				if ((!yHasNext) || (xHasNext && (x->symbols[xTerm] < y->symbols[yTerm]))) {
					// both x and y have symbol i, but only x has symbol j, so error += abs(xj * yi)
					mpfr_mul(temp, &(y->deviations[yTerm]), &(x->deviations[xNext]), MPFR_RNDA);
					mpfr_abs(temp, temp, MPFR_RNDN);
					mpfr_add(error, error, temp, MPFR_RNDU);

					xHasNext = ++xNext < x->nTerms;
				}
				else if ((!xHasNext) || (yHasNext && (y->symbols[yTerm] < x->symbols[xTerm]))) {
					// both x and y have symbol i, but only y has symbol j, so error += abs(xi * yj)
					mpfr_mul(temp, &(x->deviations[xTerm]), &(y->deviations[yNext]), MPFR_RNDA);
					mpfr_abs(temp, temp, MPFR_RNDN);
					mpfr_add(error, error, temp, MPFR_RNDU);

					yHasNext = ++yNext < y->nTerms;
				}
				else {
					// both x and y have symbols i and j, so error += abs(xi * yj + xj * yi)
					mpfr_mul(temp, &(x->deviations[xTerm]), &(y->deviations[yNext]), MPFR_RNDU);
					mpfr_mul(temp, &(x->deviations[xTerm]), &(y->deviations[yNext]), MPFR_RNDU);
					mpfr_add(temp, temp, temp, MPFR_RNDU);
					if (mpfr_sgn(temp) < 0) {
						mpfr_mul(temp, &(x->deviations[xTerm]), &(y->deviations[yNext]), MPFR_RNDD);
						mpfr_mul(temp, &(x->deviations[xTerm]), &(y->deviations[yNext]), MPFR_RNDD);
						mpfr_add(temp, temp, temp, MPFR_RNDD);
					}
					mpfr_abs(temp, temp, MPFR_RNDN);
					mpfr_add(error, error, temp, MPFR_RNDU);

					xHasNext = ++xNext < x->nTerms;
					yHasNext = ++yNext < y->nTerms;
				}
			}
			xTerm++;
			yTerm++;
		}
	}

	mpfr_max(temp, xiyiPos, xiyiNeg, MPFR_RNDN);
	mpfr_add(error, error, temp, MPFR_RNDU);

	mpfr_clear(xiyiPos);
	mpfr_clear(xiyiNeg);
#else
	mpfr_mul(temp, &(x->radius), &(y->radius), MPFR_RNDU);
	mpfr_add(error, error, temp, MPFR_RNDU);
#endif

	if (!mpfr_zero_p(error)) {
		zNew->symbols[zTerm] = mpfa_next_sym();
		mpfr_init2(&(zNew->deviations[zTerm]), prec);
		mpfr_set(&(zNew->deviations[zTerm]), error, MPFR_RNDN);
		mpfr_add(&(zNew->radius), &(zNew->radius), error, MPFR_RNDU);
		zTerm++;
	}

	zNew->nTerms = zTerm;
	if (zNew->nTerms == 0) {
		free(zNew->symbols);
		free(zNew->deviations);
	}
	else {
		zNew->symbols = realloc(zNew->symbols, zNew->nTerms * sizeof(unsigned));
		zNew->deviations = realloc(zNew->deviations, zNew->nTerms * sizeof(mpfr_t));
	}

	mpfr_clears(temp, error, (mpfr_ptr) NULL);
	mpfa_set(z, zNew);
	mpfa_clear(zNew);
}
