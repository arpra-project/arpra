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

#include "mpfa-test.h"

int main (int argc, char *argv[])
{
    mpfa_t a, b, c;

    mpfa_init(a);
    mpfa_init(b);
    mpfa_init(c);

    mpfa_set_d(a, 0.2);
    mpfa_set_d(b, 0.6);
    //mpfa_set_d(a, 3.323234235732578935);
    //mpfa_set_d(b, 3.332536490684068738946);

#ifdef WITH_MPFI
    mpfi_t p;
    mpfi_init(p);
    mpfi_set_d(p, 2.0);
    mpfa_test_cmp_mpfi(a, p);
#endif // WITH_MPFI

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

    mpfa_clear(a);
    mpfa_clear(b);
    mpfa_clear(c);

    mpfr_free_cache();

    return 0;
}
