/*
 * test.c -- Setup and teardown routines for testing.
 *
 * Copyright 2017 James Paul Turner.
 *
 * This code is adapted from the GNU MPFR library test suite.
 * Copyright 2001-2017 Free Software Foundation, Inc.
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
#include <time.h>

static gmp_randstate_t mpfa_test_randstate;
static char mpfa_test_initialised = 0;

void mpfa_test_start () {
    unsigned long int seed;
    char *environment_seed;

    if (mpfa_test_initialised) {
        fprintf(stderr, "Error: test is alreay initialised.\n");
        exit(1);
    }

    gmp_randinit_default(mpfa_test_randstate);
    mpfa_test_initialised = 1;

    environment_seed = getenv("MPFA_TEST_RAND_SEED");
    if (environment_seed != NULL) {
        seed = strtoul(environment_seed, NULL, 10);
        gmp_randseed_ui(mpfa_test_randstate, seed);
        printf("Seeding with MPFA_TEST_RAND_SEED=%lu.\n", seed);
    }
    else {
#ifdef HAVE_CLOCK_GETTIME
        struct timespec t;
        clock_gettime(CLOCK_REALTIME, &t);
        seed = t.tv_sec + t.tv_nsec;
#else // Else use stdlib clock.
        time(&seed);
#endif
        gmp_randseed_ui(mpfa_test_randstate, seed);
        printf("Seeding with %lu.\n", seed);
    }
}

void mpfa_test_end () {
    if (mpfa_test_initialised) {
        gmp_randclear(mpfa_test_randstate);
        mpfa_test_initialised = 0;
    }
    else {
        fprintf(stderr, "Error: test is not initialised.\n");
        exit(1);
    }
}
