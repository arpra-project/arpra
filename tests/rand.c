/*
 * rand.c -- Get, check, initialise and clear the RNG.
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

static gmp_randstate_t test_randstate;
static int test_rand_initialised = 0;

void test_rand_init ()
{
    unsigned long int seed;
    char *environment_seed;

    if (test_rand_initialised) {
        fprintf(stderr, "Error: RNG is alreay initialised.\n");
        exit(EXIT_FAILURE);
    }

    gmp_randinit_default(test_randstate);
    test_rand_initialised = 1;

    environment_seed = getenv("MPFA_TEST_RAND_SEED");
    if (environment_seed != NULL) {
        seed = strtoul(environment_seed, NULL, 10);
        gmp_randseed_ui(test_randstate, seed);
        printf("Seeding with MPFA_TEST_RAND_SEED=%lu.\n", seed);
    }
    else {
#ifdef HAVE_CLOCK_GETTIME
        // Seed with clock_gettime.
        struct timespec t;
        clock_gettime(CLOCK_REALTIME, &t);
        seed = t.tv_sec + t.tv_nsec;
#else
        // Else seed with stdlib clock.
        time(&seed);
#endif
        gmp_randseed_ui(test_randstate, seed);
        printf("Seeding with %lu.\n", seed);
    }
}

void test_rand_clear ()
{
    if (test_rand_initialised) {
        gmp_randclear(test_randstate);
        test_rand_initialised = 0;
    }
    else {
        fprintf(stderr, "Error: RNG is not initialised.\n");
        exit(EXIT_FAILURE);
    }
}

int test_rand_is_init ()
{
    return test_rand_initialised;
}

gmp_randstate_t *test_rand_get ()
{
    return &test_randstate;
}
