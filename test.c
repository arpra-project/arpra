/*
 * test.c
 *
 *  Created on: 21 Jul 2016
 *      Author: jt273
 */

#include <stdio.h>

#include "mpfa.h"

int main (int argc, char *argv[]) {
	mpfr_t t;

	mpfr_init2 (t, 24);
	mpfr_set_d (t, 1.19998490810394287109375, MPFR_RNDN);

	mpfr_out_str (stdout, 10, 100, t, MPFR_RNDN);
	putchar ('\n');
	mpfr_out_str (stdout, 2, 0, t, MPFR_RNDN);
	putchar ('\n');

	mpfr_clear (t);

    return 0;
}
