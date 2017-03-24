/*
 * affine_1.c -- Compute a univariate affine function of an affine form.
 *
 * Copyright 2016-2017 James Paul Turner.
 *
 * This file is part of the MPFA library.
 *
 * The MPFA library is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * The MPFA library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with the MPFA library. If not, see <http://www.gnu.org/licenses/>.
 */

#include "mpfa.h"
#include <malloc.h>
#include <assert.h>

void mpfa_affine_1 (mpfa_ptr z, mpfa_srcptr x, mpfr_srcptr alpha, mpfr_srcptr gamma, mpfr_srcptr delta) {
	unsigned xTerm, zTerm;
	mpfr_t temp, error;
	mpfr_prec_t prec;
	mpfa_t zNew;

	prec = mpfr_get_prec(&(z->centre));
	mpfr_inits2(prec, temp, error, (mpfr_ptr) NULL);
	mpfa_init2(zNew, prec);
	mpfr_set_si(error, 0, MPFR_RNDN);
	mpfr_set_si(&(zNew->radius), 0, MPFR_RNDN);

	if (mpfr_fma(&(zNew->centre), alpha, &(x->centre), gamma, MPFR_RNDN)) {
		assert(mpfr_set_si(temp, (-prec + mpfr_get_exp(&(zNew->centre))), MPFR_RNDN) == 0);
		assert(mpfr_exp2(temp, temp, MPFR_RNDN) == 0);
		mpfr_add(error, error, temp, MPFR_RNDU);
	}

	zNew->nTerms = x->nTerms + 1;
	zNew->symbols = malloc(zNew->nTerms * sizeof(unsigned));
	zNew->deviations = malloc(zNew->nTerms * sizeof(mpfr_t));

	for (xTerm = 0, zTerm = 0; xTerm < x->nTerms; xTerm++) {
		zNew->symbols[zTerm] = x->symbols[xTerm];
		mpfr_init2(&(zNew->deviations[zTerm]), prec);

		if (mpfr_mul(&(zNew->deviations[zTerm]), alpha, &(x->deviations[xTerm]), MPFR_RNDN)) {
			assert(mpfr_set_si(temp, (-prec + mpfr_get_exp(&(zNew->deviations[zTerm]))), MPFR_RNDN) == 0);
			assert(mpfr_exp2(temp, temp, MPFR_RNDN) == 0);
			mpfr_add(error, error, temp, MPFR_RNDU);
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

	mpfr_add(error, error, delta, MPFR_RNDU);
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
