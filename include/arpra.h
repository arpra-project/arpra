/*
 * arpra.h -- Arpra main public header.
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

#ifndef ARPRA_H
#define ARPRA_H

#include <mpfr.h>

// Primitive type definitions.
typedef long int arpra_int;
typedef unsigned long int arpra_uint;
typedef mpfr_prec_t arpra_precision;

// The Arpra range type.
struct arpra_range
{
    __mpfr_struct centre;
    __mpfr_struct radius;
    __mpfr_struct *deviations;
    arpra_uint *symbols;
    arpra_uint nTerms;
};

#ifdef __cplusplus
extern "C" {
#endif

// Initialise and clear.
void arpra_init (struct arpra_range *x);
void arpra_inits (struct arpra_range *x, ...);
void arpra_init2 (struct arpra_range *x, const arpra_precision prec);
void arpra_inits2 (const arpra_precision prec, struct arpra_range *x, ...);
void arpra_clear (struct arpra_range *x);
void arpra_clears (struct arpra_range *x, ...);
void arpra_clear_terms (struct arpra_range *x);

// Get from an Arpra range.
void arpra_get_bounds (mpfr_ptr lo, mpfr_ptr hi, const struct arpra_range *x);

// Set to an Arpra range.
void arpra_set (struct arpra_range *z, const struct arpra_range *x);
void arpra_set_d (struct arpra_range *z, const double centre);
void arpra_set_d_rad (struct arpra_range *z, const double centre, const double radius);
void arpra_set_str (struct arpra_range *z, const char *centre, const arpra_int base);
void arpra_set_str_rad (struct arpra_range *z, const char *centre, const char *radius, const arpra_int base);
void arpra_set_mpfr (struct arpra_range *z, mpfr_srcptr centre);
void arpra_set_mpfr_rad (struct arpra_range *z, mpfr_srcptr centre, mpfr_srcptr radius);

// Set special values.
void arpra_set_nan (struct arpra_range *z);
void arpra_set_inf (struct arpra_range *z);
void arpra_set_zero (struct arpra_range *z);

// Affine operations.
void arpra_affine_1 (struct arpra_range *z, const struct arpra_range *x,
                     mpfr_srcptr alpha, mpfr_srcptr gamma, mpfr_srcptr delta);
void arpra_affine_2 (struct arpra_range *z, const struct arpra_range *x, const struct arpra_range *y,
                     mpfr_srcptr alpha, mpfr_srcptr beta, mpfr_srcptr gamma, mpfr_srcptr delta);
void arpra_add (struct arpra_range *z, const struct arpra_range *x, const struct arpra_range *y);
void arpra_sub (struct arpra_range *z, const struct arpra_range *x, const struct arpra_range *y);
void arpra_neg (struct arpra_range *z, const struct arpra_range *x);
void arpra_sum (struct arpra_range *z, struct arpra_range **x, const arpra_uint n);

// Non-affine operations.
void arpra_mul (struct arpra_range *z, const struct arpra_range *x, const struct arpra_range *y);
void arpra_div (struct arpra_range *z, const struct arpra_range *x, const struct arpra_range *y);
void arpra_sqrt (struct arpra_range *z, const struct arpra_range *x);
void arpra_exp (struct arpra_range *z, const struct arpra_range *x);
void arpra_log (struct arpra_range *z, const struct arpra_range *x);
void arpra_inv (struct arpra_range *z, const struct arpra_range *x);

// Numerical precision.
arpra_precision arpra_get_precision (const struct arpra_range *x);
arpra_precision arpra_get_default_precision ();
arpra_precision arpra_get_internal_precision ();
void arpra_set_precision (struct arpra_range *z, const arpra_precision precision);
void arpra_set_default_precision (const arpra_precision precision);
void arpra_set_internal_precision (const arpra_precision precision);

// Deviation term reduction.
void arpra_reduce_last_n (struct arpra_range *z, arpra_uint n);
void arpra_reduce_small (struct arpra_range *z, double min_fraction);

// Predicates on Arpra ranges.
int arpra_nan_p (const struct arpra_range *x);
int arpra_inf_p (const struct arpra_range *x);
int arpra_bounded_p (const struct arpra_range *x);
int arpra_zero_p (const struct arpra_range *x);
int arpra_has_zero_p (const struct arpra_range *x);
int arpra_has_pos_p (const struct arpra_range *x);
int arpra_has_neg_p (const struct arpra_range *x);

// Get new deviation symbols.
arpra_uint arpra_next_symbol ();

#ifdef __cplusplus
}
#endif

#endif // ARPRA_H
