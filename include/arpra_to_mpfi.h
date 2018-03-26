/*
 * arpra_to_mpfi.h -- Arpra public header file for MPFI support.
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

#ifndef ARPRA_TO_MPFI_H
#define ARPRA_TO_MPFI_H

#include <mpfi.h>

#ifdef __cplusplus
extern "C" {
#endif

// Get and set MPFI intervals.
void arpra_get_mpfi (mpfi_ptr z, const struct arpra_range *x);
void arpra_set_mpfi (struct arpra_range *z, mpfi_srcptr x);

#ifdef __cplusplus
}
#endif

#endif // ARPRA_TO_MPFI_H
