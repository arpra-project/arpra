/*
 * set_special.c -- Set special values.
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

#include "arpra-impl.h"

void arpra_set_nan (struct arpra_range *z)
{
    // Clear existing deviation terms.
    arpra_clear_terms(z);

    // Set centre to NaN and radius to NaN.
    mpfr_set_nan(&(z->centre));
    mpfr_set_nan(&(z->radius));
}

void arpra_set_inf (struct arpra_range *z)
{
    // Clear existing deviation terms.
    arpra_clear_terms(z);

    // Set centre to +0 and radius to +Inf.
    mpfr_set_zero(&(z->centre), 1);
    mpfr_set_inf(&(z->radius), 1);
}

void arpra_set_zero (struct arpra_range *z)
{
    // Clear existing deviation terms.
    arpra_clear_terms(z);

    // Set centre to +0 and radius to +0.
    mpfr_set_zero(&(z->centre), 1);
    mpfr_set_zero(&(z->radius), 1);
}
