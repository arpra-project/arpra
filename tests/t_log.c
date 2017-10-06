/*
 * t_log.c -- Test the mpfa_log function.
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
#ifdef WITH_MPFI
    const mpfa_prec_t prec = 53;
    const mpfa_prec_t prec_internal = 128;
    const mpfa_uint_t test_n = 100000;
    mpfa_uint_t i, fail_n;

    // Init test.
    test_fixture_init(prec, prec_internal);
    test_rand_init();
    test_log_init("log");
    fail_n = 0;

    // Start test.
    for (i = 0; i < test_n; i++) {
        test_rand_mpfa(x_A, TEST_RAND_POS, TEST_RAND_SMALL);
        if (test_univariate_mpfi(mpfa_log, mpfi_log)) {
            fail_n++;
        }
    }

    // Cleanup test.
    printf("%llu out of %llu failed.\n", fail_n, test_n);
    test_fixture_clear();
    test_rand_clear();
    test_log_clear();
    return fail_n > 0;

#else // WITH_MPFI
    fprintf(stderr,
            "This test uses the MPFI interval arithmetic library.\n"
            "Recompile with MPFI support enabled to run this test.\n");
    return 77; // Exit code 77: skip test.
#endif // WITH_MPFI
}
