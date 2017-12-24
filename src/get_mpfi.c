/*
 * get_mpfi.c -- Get an MPFI interval from an arpra_t.
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

#include "arpra-impl.h"

void arpra_get_mpfi (mpfi_ptr z, arpra_srcptr x)
{
    // z_lo = x_0 - rad(x)
    mpfr_sub(&(z->left), &(x->centre), &(x->radius), MPFR_RNDD);

    // z_hi = x_0 + rad(x)
    mpfr_add(&(z->right), &(x->centre), &(x->radius), MPFR_RNDU);
}
