/*
 * div.c -- Divide one arpra_t by another.
 *
 * Copyright 2016-2017 James Paul Turner.
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

#include "arpra-impl.h"

/*
 * For now we just multiply the numerator with the reciprocal of the denominator.
 */

void arpra_div (arpra_ptr z, arpra_srcptr x, arpra_srcptr y)
{
    arpra_t zNew;
    arpra_prec_t prec;

    // Init temp z.
    prec = arpra_get_prec(z);
    arpra_init2(zNew, prec);

    // z = x * (1 / y)
    arpra_inv(zNew, y);
    arpra_mul(z, x, zNew);

    // Clear temp z.
    arpra_clear(zNew);
}
