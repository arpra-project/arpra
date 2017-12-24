/*
 * logfile_mpfr.c -- Record an MPFR result in the logfile.
 *
 * Copyright 2017 James Paul Turner.
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

void test_log_mpfr (mpfr_srcptr x, const char *var_name)
{
    if (test_log_ready) {
        // Write variable name and value to logfile.
        fprintf(test_log, "%s: ", var_name);
        mpfr_out_str(test_log, 10, 80, x, MPFR_RNDN);
        fputs("\n", test_log);
    }
    else {
        fprintf(stderr, "Error: Logfile is not initialised.\n");
        exit(EXIT_FAILURE);
    }
}
