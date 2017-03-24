/*
 * set_d.c -- Set an affine form using a double-precision float.
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

void mpfa_set_d (mpfa_ptr z, const double centre, const double radius) {
	unsigned zTerm;
	mpfr_prec_t prec;

	prec = mpfr_get_prec(&(z->centre));

	if (mpfr_set_d(&(z->centre), centre, MPFR_RNDN)) {
		assert(mpfr_set_si(&(z->radius), (-prec + mpfr_get_exp(&(z->centre))), MPFR_RNDN) == 0);
		assert(mpfr_exp2(&(z->radius), &(z->radius), MPFR_RNDN) == 0);
		mpfr_add_d(&(z->radius), &(z->radius), radius, MPFR_RNDU);
	}
	else {
		mpfr_set_d(&(z->radius), radius, MPFR_RNDU);
	}

	if (mpfr_zero_p(&(z->radius))) {
		if (z->nTerms > 0) {
			for (zTerm = 0; zTerm < z->nTerms; zTerm++) {
				mpfr_clear(&(z->deviations[zTerm]));
			}
			z->nTerms = 0;
			free(z->symbols);
			free(z->deviations);
		}
	}
	else {
		if (z->nTerms == 0) {
			z->symbols = malloc(sizeof(unsigned));
			z->deviations = malloc(sizeof(mpfa_t));
			mpfr_init2(&(z->deviations[0]), prec);
		}
		else if (z->nTerms >= 2) {
			for (zTerm = 1; zTerm < z->nTerms; zTerm++) {
				mpfr_clear(&(z->deviations[zTerm]));
			}
			z->symbols = realloc(z->symbols, sizeof(unsigned));
			z->deviations = realloc(z->deviations, sizeof(mpfa_t));
		}
		z->nTerms = 1;
		z->symbols[0] = mpfa_next_sym();
		mpfr_set_d(&(z->deviations[0]), radius, MPFR_RNDN);
	}
}
