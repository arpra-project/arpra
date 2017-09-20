/*
 * logfile_mpfi.c -- Record an MPFI result in the logfile.
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

void test_log_mpfi (mpfi_srcptr x, const char *var_name)
{
    if (test_log_ready) {
        // Record variable name and bounds in logfile.
        fprintf(test_log, "%s: ", var_name);
        mpfr_out_str(test_log, 10, 80, &(x->left), MPFR_RNDN);
        fputs(" ", test_log);
        mpfr_out_str(test_log, 10, 80, &(x->right), MPFR_RNDN);
        fputs("\n", test_log);
    }
    else {
        fprintf(stderr, "Error: Logfile is not initialised.\n");
        exit(EXIT_FAILURE);
    }
}
