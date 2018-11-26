/*
 * predicates.c -- Predicates on Arpra ranges.
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

#include "arpra-impl.h"

int arpra_nan_p (const arpra_range *x)
{
    return mpfi_nan_p(&(x->true_range));
}

int arpra_inf_p (const arpra_range *x)
{
    return mpfi_inf_p(&(x->true_range));
}

int arpra_bounded_p (const arpra_range *x)
{
    return mpfi_bounded_p(&(x->true_range));
}

int arpra_zero_p (const arpra_range *x)
{
    return mpfr_zero_p(&(x->true_range.left)) && mpfr_zero_p(&(x->true_range.right));
}

int arpra_has_zero_p (const arpra_range *x)
{
    return !mpfi_nan_p(&(x->true_range))
           && (mpfr_sgn(&(x->true_range.left)) <= 0) && (mpfr_sgn(&(x->true_range.right)) >= 0);
}

int arpra_has_pos_p (const arpra_range *x)
{
    return mpfr_sgn(&(x->true_range.right)) > 0;
}

int arpra_has_neg_p (const arpra_range *x)
{
    return mpfr_sgn(&(x->true_range.left)) < 0;
}
