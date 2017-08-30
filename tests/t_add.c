/*
 * t_add.c -- Test the MPFA add function.
 *
 * Copyright 2017 James Paul Turner.
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
    mpfa_t a, b;

    mpfa_test_rand_init();
    mpfa_init(a);
    mpfa_init(b)

#ifdef WITH_MPFI
    mpfi_t a_i, b_i;

    mpfi_init(a_i);
    mpfi_init(b_i);
#endif // WITH_MPFI







    mpfa_test_rand_clear();
    mpfa_clear(a);
    mpfa_clear(b);

#ifdef WITH_MPFI
    mpfi_init(a_i);
    mpfi_init(b_i);
#endif // WITH_MPFI

    mpfr_clear_cache();
}
