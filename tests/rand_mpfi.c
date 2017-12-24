/*
 * rand_mpfi.c -- Generate a random MPFI interval.
 *
 * Copyright 2017 James Paul Turner.
 *
 * This file is part of the ArPRA library.
 *
 * The ArPRA library is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * The ArPRA library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with the ArPRA library. If not, see <http://www.gnu.org/licenses/>.
 */

#include "arpra-test.h"

void test_rand_mpfi (mpfi_ptr z,
                     enum test_rand_mode mode_low,
                     enum test_rand_mode mode_high)
{
    // Set random lower and upper bound.
    test_rand_mpfr(&(z->left), mode_low);
    test_rand_mpfr(&(z->right), mode_high);
}
