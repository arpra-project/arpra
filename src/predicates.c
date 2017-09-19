/*
 * predicates.c -- Predicates on affine forms.
 *
 * Copyright 2017 James Paul Turner.
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

#include "mpfa-impl.h"

int mpfa_nan_p (mpfa_srcptr x)
{
    return mpfr_nan_p(&(x->centre)) || mpfr_nan_p(&(x->radius));
}

int mpfa_inf_p (mpfa_srcptr x)
{
    return !mpfr_nan_p(&(x->centre)) && mpfr_inf_p(&(x->radius));
}

int mpfa_zero_p (mpfa_srcptr x)
{
    return mpfr_zero_p(&(x->centre)) && mpfr_zero_p(&(x->radius));
}

int mpfa_has_zero_p (mpfa_srcptr x)
{
    return mpfr_cmpabs(&(x->centre), &(x->radius)) <= 0;
}
