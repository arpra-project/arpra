/*
 * clear.c -- Free the memory of an affine form.
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
#include <stdlib.h>

void mpfa_clear (mpfa_ptr x) {
	unsigned xTerm;

	if (x->nTerms > 0) {
		for (xTerm = 0; xTerm < x->nTerms; xTerm++) {
			mpfr_clear(&(x->deviations[xTerm]));
		}
		free(x->symbols);
		free(x->deviations);
	}

	mpfr_clear(&(x->centre));
	mpfr_clear(&(x->radius));
}
