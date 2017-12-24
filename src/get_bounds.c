/*
 * get_bounds.c -- Get the lower and upper bound of an arpra_t.
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

#include "arpra-impl.h"

void arpra_get_bounds (mpfr_ptr lo, mpfr_ptr hi, arpra_srcptr x)
{
    // lo = x_0 - rad(x)
    mpfr_sub(lo, &(x->centre), &(x->radius), MPFR_RNDD);

    // hi = x_0 + rad(x)
    mpfr_add(hi, &(x->centre), &(x->radius), MPFR_RNDU);
}
