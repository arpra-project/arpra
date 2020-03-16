/*
 * logfile.c -- Initialise and clear the logfile.
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

int test_log_ready = 0;
FILE *test_log;

void test_log_init (const char *test_name)
{
    char filename[100];

    // Ensure that we do not double-initialise.
    if (!test_log_ready) {
        test_log_ready = 1;

        strcpy(filename, test_name);
        strcat(filename, ".log");
        test_log = fopen(filename, "w");

        // Check that logfile opened OK.
        if (test_log == NULL) {
            fprintf(stderr, "Error: logfile cannot be opened.\n");
            exit(EXIT_FAILURE);
        }
    }
    else {
        fprintf(stderr, "Error: logfile is alreay initialised.\n");
        exit(EXIT_FAILURE);
    }
}

void test_log_clear ()
{
    // Ensure that we do not double-clear.
    if (test_log_ready) {
        test_log_ready = 0;

        // Check that logfile closed OK.
        if (fclose(test_log)) {
            fprintf(stderr, "Error: logfile cannot be closed.\n");
            exit(EXIT_FAILURE);
        }
    }
    else {
        fprintf(stderr, "Error: logfile is not initialised.\n");
        exit(EXIT_FAILURE);
    }
}
