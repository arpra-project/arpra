/*
 * div.c -- Divide one affine form by another.
 *
 * Copyright 2016-2017 James Paul Turner.
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

#include "mpfa.h"

/*
 * For now we just multiply the numerator with the reciprocal of the denominator.
 *
 * TODO: find a better way to do mpfa_div
 */

void mpfa_div (mpfa_ptr z, mpfa_srcptr x, mpfa_srcptr y) {
    mpfa_t zNew;
    mpfa_prec_t prec;

    prec = mpfa_get_prec(z);
    mpfa_init2(zNew, prec);

    mpfa_inv(zNew, y);
    mpfa_mul(z, x, zNew);

    mpfa_clear(zNew);
}
