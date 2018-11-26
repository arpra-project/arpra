/*
 * div.c -- Divide one Arpra range by another.
 *
 * Copyright 2016-2018 James Paul Turner.
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

/*
 * For now we just multiply the numerator with the reciprocal of the denominator.
 */

void arpra_div (arpra_range *z, const arpra_range *x, const arpra_range *y)
{
    arpra_range zNew;
    arpra_mpfi z_range;

    // Initialise vars.
    arpra_init2(&zNew, z->precision);
    mpfi_init2(&z_range, z->precision);

    // MPFI division
    mpfi_div(&z_range, &(x->true_range), &(y->true_range));

    // z = x * (1 / y)
    arpra_inv(&zNew, y);
    arpra_mul(z, x, &zNew);

    // Compute true range.
    mpfi_intersect(&(z->true_range), &(z->true_range), &z_range);

    // Clear vars.
    arpra_clear(&zNew);
    mpfi_clear(&z_range);
}
