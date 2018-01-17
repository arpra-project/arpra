/*
 * logfile_printf.c -- Record formatted text in the logfile.
 *
 * Copyright 2017-2018 James Paul Turner.
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

void test_log_printf (const char *format, ...)
{
    va_list arg;

    if (test_log_ready) {
        // Write formatted text to logfile.
        va_start(arg, format);
        vfprintf(test_log, format, arg);
        va_end(arg);
    }
    else {
        fprintf(stderr, "Error: Logfile is not initialised.\n");
        exit(EXIT_FAILURE);
    }
}
