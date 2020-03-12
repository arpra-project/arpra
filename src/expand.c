/*
 * expand.c -- Add to the last deviation term of an Arpra range.
 *
 * Copyright 2020 James Paul Turner.
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


// expand.c: arpra_expand_mpfr, alias arpra_expand
// ARG SHOULD BE NON-NEGATIVE
// SAME FOR OTHER TYPES (using mpfr_add_<type>)

#define ARPRA_EXPAND_FN(SIGNATURE)              \
    void SIGNATURE                              \
    {                                           \
                                                \
    }


#define FN_SIGNATURE(FN_TYPE, ...)                                      \
    arpra_expand_##FN_TYPE (arpra_range *y, __VA_ARGS__)

#define FN_MPFR_CALL(...)                                               \
    ARPRA_MPFR_RNDERR(error, MPFR_RNDN, fn, &(yy.centre), __VA_ARGS__)
