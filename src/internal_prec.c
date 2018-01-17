/*
 * internal_prec.c -- Get and set the precision used internally by Arpra.
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

static arpra_prec_t arpra_internal_prec = ARPRA_DEFAULT_INTERNAL_PREC;

arpra_prec_t arpra_get_internal_prec ()
{
    return arpra_internal_prec;
}

void arpra_set_internal_prec (arpra_prec_t prec)
{
    arpra_internal_prec = prec;
}
