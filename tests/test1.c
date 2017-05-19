/*
 * test1.c -- Temporary testing, until a full test suite is implemented.
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

#include <stdio.h>
#include <assert.h>
#include "mpfa.h"


int main (int argc, char *argv[]) {

	// MPFA testing ============================================

	mpfa_t a, b, c;

	mpfa_init(a);
	mpfa_init(b);
	mpfa_init(c);

//	mpfr_t temp, error;
//	mpfr_inits(temp, error, (mpfr_ptr) NULL);
//	assert(!mpfr_mul_si(temp, &(a->u), (-1ul), MPFR_RNDU));
//	assert(!mpfr_si_sub(error, 1, temp, MPFR_RNDD));
//	assert(!mpfr_div(error, temp, error, MPFR_RNDU));
//	mpfr_clears(temp, error, (mpfr_ptr) NULL);
//	return 0;

	mpfa_set_d(a, 0.2, 0.0);
	mpfa_set_d(b, 0.6, 0.0);
	//mpfa_set_d(a, 3.323234235732578935, 0.0);
	//mpfa_set_d(b, 3.332536490684068738946, 0.0);

	mpfa_add(c, a, b);
	//printf("centre: "); mpfr_out_str (stdout, 10, 100, &(c->centre), MPFR_RNDN); putchar('\n');
	//printf("radius: "); mpfr_out_str (stdout, 10, 100, &(c->radius), MPFR_RNDN); putchar('\n');

	mpfa_mul(c, a, b);
	//printf("centre: "); mpfr_out_str (stdout, 10, 100, &(c->centre), MPFR_RNDN); putchar('\n');
	//printf("radius: "); mpfr_out_str (stdout, 10, 100, &(c->radius), MPFR_RNDN); putchar('\n');

	mpfa_div(c, a, b);
	//printf("centre: "); mpfr_out_str (stdout, 10, 100, &(c->centre), MPFR_RNDN); putchar('\n');
	//printf("radius: "); mpfr_out_str (stdout, 10, 100, &(c->radius), MPFR_RNDN); putchar('\n');

	mpfa_exp(c, a);
	//printf("centre: "); mpfr_out_str (stdout, 10, 100, &(c->centre), MPFR_RNDN); putchar('\n');
	//printf("radius: "); mpfr_out_str (stdout, 10, 100, &(c->radius), MPFR_RNDN); putchar('\n');

	mpfa_inv(c, a);
	//printf("centre: "); mpfr_out_str (stdout, 10, 100, &(c->centre), MPFR_RNDN); putchar('\n');
	//printf("radius: "); mpfr_out_str (stdout, 10, 100, &(c->radius), MPFR_RNDN); putchar('\n');



	//*
	printf("centre: "); mpfr_out_str (stdout, 10, 100, &(a->centre), MPFR_RNDN); putchar('\n');
	printf("radius: "); mpfr_out_str (stdout, 10, 100, &(a->radius), MPFR_RNDN); putchar('\n');

	int i;
	for (i = 0; i < 289; i++) {
		//mpfa_mul(a, a, b);
		//mpfa_mul(b, a, b);
		mpfa_div(a, a, b);
		printf("centre: "); mpfr_out_str (stdout, 10, 100, &(a->centre), MPFR_RNDN); putchar('\n');
		printf("radius: "); mpfr_out_str (stdout, 10, 100, &(a->radius), MPFR_RNDN); putchar('\n');
	}
	//*/

	mpfa_clear(a);
	mpfa_clear(b);
	mpfa_clear(c);




/* MPFR testing ============================================

	mpfr_t t;

	mpfr_init2 (t, 24);
	mpfr_set_d (t, 1.19998490810394287109375, MPFR_RNDN);

	mpfr_out_str (stdout, 10, 100, t, MPFR_RNDN);
	putchar ('\n');
	mpfr_out_str (stdout, 2, 0, t, MPFR_RNDN);
	putchar ('\n');

	mpfr_clear (t);
//*/

	mpfr_free_cache();

	return 0;
}
