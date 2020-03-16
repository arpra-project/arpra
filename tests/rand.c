/*
 * rand.c -- Intialise and clear the RNG.
 *
 * Copyright 2017-2020 James Paul Turner.
 *
 * This file is part of the Arpra library.
 *
 * The Arpra library is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * The Arpra library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with the Arpra library. If not, see <http://www.gnu.org/licenses/>.
 */

#include "arpra-test.h"

int test_rand_ready = 0;
gmp_randstate_t test_randstate;

void test_rand_init ()
{
    unsigned long int seed;
    char *environment_seed;

    // Ensure that we do not double-initialise.
    if (!test_rand_ready) {
        test_rand_ready = 1;

        // Try to get seed from environment.
        gmp_randinit_default(test_randstate);
        environment_seed = getenv("ARPRA_TEST_RAND_SEED");
        if (environment_seed != NULL) {
            seed = strtoul(environment_seed, NULL, 10);
            gmp_randseed_ui(test_randstate, seed);
            printf("Seeding with ARPRA_TEST_RAND_SEED=%lu.\n", seed);
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
    else {
        fprintf(stderr, "Error: RNG is alreay initialised.\n");
        exit(EXIT_FAILURE);
    }
}

void test_rand_clear ()
{
    // Ensure that we do not double-clear.
    if (test_rand_ready) {
        test_rand_ready = 0;
        gmp_randclear(test_randstate);
    }
    else {
        fprintf(stderr, "Error: RNG is not initialised.\n");
        exit(EXIT_FAILURE);
    }
}
