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

// Arpra typedefs.
typedef struct arpra_range_struct arpra_range;
typedef __mpfr_struct arpra_mpfr;
typedef mpfr_prec_t arpra_precision;
typedef long int arpra_int;
typedef unsigned long int arpra_uint;

// The Arpra range type.
struct arpra_range_struct
{
    arpra_mpfr centre;
    arpra_mpfr radius;
    arpra_mpfr *deviations;
    arpra_uint *symbols;
    arpra_uint nTerms;
};

#ifdef __cplusplus
extern "C" {
#endif

// Initialise and clear.
void arpra_init (arpra_range *x);
void arpra_inits (arpra_range *x, ...);
void arpra_init2 (arpra_range *x, const arpra_precision prec);
void arpra_inits2 (const arpra_precision prec, arpra_range *x, ...);
void arpra_clear (arpra_range *x);
void arpra_clears (arpra_range *x, ...);
void arpra_clear_terms (arpra_range *x);

// Get from an Arpra range.
void arpra_get_bounds (arpra_mpfr *lo, arpra_mpfr *hi, const arpra_range *x);

// Set to an Arpra range.
void arpra_set (arpra_range *z, const arpra_range *x);
void arpra_set_d (arpra_range *z, const double centre);
void arpra_set_d_rad (arpra_range *z, const double centre, const double radius);
void arpra_set_str (arpra_range *z, const char *centre, const arpra_int base);
void arpra_set_str_rad (arpra_range *z, const char *centre, const char *radius, const arpra_int base);
void arpra_set_mpfr (arpra_range *z, const arpra_mpfr *centre);
void arpra_set_mpfr_rad (arpra_range *z, const arpra_mpfr *centre, const arpra_mpfr *radius);

// Set special values.
void arpra_set_nan (arpra_range *z);
void arpra_set_inf (arpra_range *z);
void arpra_set_zero (arpra_range *z);

// Affine operations.
void arpra_affine_1 (arpra_range *z, const arpra_range *x,
                     const arpra_mpfr *alpha, const arpra_mpfr *gamma,
                     const arpra_mpfr *delta);
void arpra_affine_2 (arpra_range *z, const arpra_range *x, const arpra_range *y,
                     const arpra_mpfr *alpha, const arpra_mpfr *beta,
                     const arpra_mpfr *gamma, const arpra_mpfr *delta);
void arpra_add (arpra_range *z, const arpra_range *x, const arpra_range *y);
void arpra_sub (arpra_range *z, const arpra_range *x, const arpra_range *y);
void arpra_neg (arpra_range *z, const arpra_range *x);
void arpra_sum (arpra_range *z, arpra_range *x, const arpra_uint n);

// Non-affine operations.
void arpra_mul (arpra_range *z, const arpra_range *x, const arpra_range *y);
void arpra_div (arpra_range *z, const arpra_range *x, const arpra_range *y);
void arpra_sqrt (arpra_range *z, const arpra_range *x);
void arpra_exp (arpra_range *z, const arpra_range *x);
void arpra_log (arpra_range *z, const arpra_range *x);
void arpra_inv (arpra_range *z, const arpra_range *x);

// Numerical precision.
arpra_precision arpra_get_precision (const arpra_range *x);
arpra_precision arpra_get_default_precision ();
arpra_precision arpra_get_internal_precision ();
void arpra_set_precision (arpra_range *z, const arpra_precision prec);
void arpra_set_default_precision (const arpra_precision prec);
void arpra_set_internal_precision (const arpra_precision prec);

// Deviation term reduction.
void arpra_reduce_last_n (arpra_range *z, arpra_uint n);
void arpra_reduce_small (arpra_range *z, double min_fraction);

// Predicates on Arpra ranges.
int arpra_nan_p (const arpra_range *x);
int arpra_inf_p (const arpra_range *x);
int arpra_bounded_p (const arpra_range *x);
int arpra_zero_p (const arpra_range *x);
int arpra_has_zero_p (const arpra_range *x);
int arpra_has_pos_p (const arpra_range *x);
int arpra_has_neg_p (const arpra_range *x);

// Get new deviation symbols.
arpra_uint arpra_next_symbol ();

#ifdef __cplusplus
}
#endif

#endif // ARPRA_H
