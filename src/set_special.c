/*
 * set_special.c -- Set special values.
 *
 * Copyright 2017-2020 James Paul Turner.
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

void arpra_set_nan (arpra_range *y)
{
    arpra_prec prec_internal;

    // Initialise vars.
    prec_internal = arpra_get_internal_precision();
    mpfr_set_prec(&(y->centre), prec_internal);
    mpfr_set_prec(&(y->radius), prec_internal);
    arpra_helper_clear_terms(y);

    // y[0] = NaN
    mpfr_set_nan(&(y->centre));

    // Allocate memory for deviation terms.
    y->symbols = malloc(sizeof(arpra_uint));
    y->deviations = malloc(sizeof(mpfr_t));

    // Store new deviation term.
    y->symbols[0] = arpra_helper_next_symbol();
    mpfr_init2(&(y->deviations[0]), prec_internal);
    mpfr_set_nan(&(y->deviations[0]));
    mpfr_set_nan(&(y->radius));
    y->nTerms = 1;

    // Compute true_range.
    mpfr_set_nan(&(y->true_range.left));
    mpfr_set_nan(&(y->true_range.right));
}

void arpra_set_inf (arpra_range *y)
{
    arpra_prec prec_internal;

    // Initialise vars.
    prec_internal = arpra_get_internal_precision();
    mpfr_set_prec(&(y->centre), prec_internal);
    mpfr_set_prec(&(y->radius), prec_internal);
    arpra_helper_clear_terms(y);

    // y[0] = Inf
    mpfr_set_zero(&(y->centre), 1);

    // Allocate memory for deviation terms.
    y->symbols = malloc(sizeof(arpra_uint));
    y->deviations = malloc(sizeof(mpfr_t));

    // Store new deviation term.
    y->symbols[0] = arpra_helper_next_symbol();
    mpfr_init2(&(y->deviations[0]), prec_internal);
    mpfr_set_inf(&(y->deviations[0]), 1);
    mpfr_set_inf(&(y->radius), 1);
    y->nTerms = 1;

    // Compute true_range.
    mpfr_set_inf(&(y->true_range.left), -1);
    mpfr_set_inf(&(y->true_range.right), 1);
}

void arpra_set_zero (arpra_range *y)
{
    arpra_prec prec_internal;

    // Initialise vars.
    prec_internal = arpra_get_internal_precision();
    mpfr_set_prec(&(y->centre), prec_internal);
    mpfr_set_prec(&(y->radius), prec_internal);
    arpra_helper_clear_terms(y);

    // y[0] = 0
    mpfr_set_zero(&(y->centre), 1);

    // Allocate memory for deviation terms.
    y->symbols = malloc(sizeof(arpra_uint));
    y->deviations = malloc(sizeof(mpfr_t));

    // Store new deviation term.
    y->symbols[0] = arpra_helper_next_symbol();
    mpfr_init2(&(y->deviations[0]), prec_internal);
    mpfr_set_zero(&(y->deviations[0]), 1);
    mpfr_set_zero(&(y->radius), 1);
    y->nTerms = 1;

    // Compute true_range.
    mpfr_set_zero(&(y->true_range.left), -1);
    mpfr_set_zero(&(y->true_range.right), 1);
}
