/*
 * test.c
 *
 *  Created on: 21 Jul 2016
 *      Author: jt273
 */

#include <stdio.h>
#include "mpfa.h"

#include <malloc.h>

int main (int argc, char *argv[]) {

/*
	mpfa_t c;
	c->deviations = malloc(1 * sizeof(mpfr_t));
	c->deviations = realloc(c->deviations, 1 * sizeof(mpfr_t));
	free(c->deviations);
*/

/*
	__mpfr_struct *c;
	c = malloc(1 * sizeof(mpfr_t));
	c = realloc(c, 1 * sizeof(mpfr_t));
	free(c);
*/



///*
	mpfa_t a, b, c;

	mpfa_init(a);
	mpfa_init(b);
	mpfa_init(c);

	mpfa_set_d(a, 0.2);
	mpfa_set_d(b, 0.6);

	//mpfa_add(c, a, b);
	//mpfr_out_str (stdout, 10, 100, &(c->centre), MPFR_RNDN);
	//putchar ('\n');

	//mpfa_mul(c, a, b);
	//mpfr_out_str (stdout, 10, 100, &(c->centre), MPFR_RNDN);
	//putchar ('\n');



	//c->symbols = realloc(c->symbols, 2 * sizeof(unsigned));
	//c->deviations = realloc(c->deviations, 2 * sizeof(mpfr_t));




	mpfa_div(c, a, b);
	mpfr_out_str (stdout, 10, 100, &(c->centre), MPFR_RNDN);
	putchar ('\n');







/*
	if (c->nTerms > 0) {
		unsigned xTerm;
		for (xTerm = 0; xTerm < c->nTerms; xTerm++) {
			mpfr_clear(&(c->deviations[xTerm]));
		}
		free(c->symbols);
		free(c->deviations);
	}
	mpfr_clear(&(c->centre));
	mpfr_clear(&(c->radius));
*/


	mpfa_clear(a);
	mpfa_clear(b);
	mpfa_clear(c);
//*/










	// MPFR testing ============================================

	/*
	mpfr_t t;

	mpfr_init2 (t, 24);
	mpfr_set_d (t, 1.19998490810394287109375, MPFR_RNDN);

	mpfr_out_str (stdout, 10, 100, t, MPFR_RNDN);
	putchar ('\n');
	mpfr_out_str (stdout, 2, 0, t, MPFR_RNDN);
	putchar ('\n');

	mpfr_clear (t);
	 */

    return 0;
}
