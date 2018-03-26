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

void arpra_div (struct arpra_range *z, const struct arpra_range *x, const struct arpra_range *y)
{
    arpra_precision prec;
    struct arpra_range zNew;

    // Init temp z.
    prec = arpra_get_precision(z);
    arpra_init2(&zNew, prec);

    // z = x * (1 / y)
    arpra_inv(&zNew, y);
    arpra_mul(z, x, &zNew);

    // Clear temp z.
    arpra_clear(&zNew);
}
