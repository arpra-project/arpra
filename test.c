/*
 * test.c
 *
 *  Created on: 21 Jul 2016
 *      Author: jt273
 */

#include <stdio.h>
#include "mpfa.h"


int main (int argc, char *argv[]) {

	// MPFA testing ============================================

	mpfa_t a, b, c;

	mpfa_init(a);
	mpfa_init(b);
	mpfa_init(c);

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
	for (i = 0; i < 10; i++) {
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
