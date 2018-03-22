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
struct __arpra_struct
{
    __mpfr_struct centre;
    __mpfr_struct radius;
    __mpfr_struct *deviations;
    arpra_uint *symbols;
    arpra_uint nTerms;
};




typedef struct __arpra_struct arpra_t[1];
typedef struct __arpra_struct *arpra_ptr;
typedef const struct __arpra_struct *arpra_srcptr;




#ifdef __cplusplus
extern "C" {
#endif

// Initialise and clear.
void arpra_init (arpra_ptr x);
void arpra_inits (arpra_ptr x, ...);
void arpra_init2 (arpra_ptr x, arpra_precision prec);
void arpra_inits2 (arpra_precision prec, arpra_ptr x, ...);
void arpra_clear (arpra_ptr x);
void arpra_clears (arpra_ptr x, ...);
void arpra_clear_terms (arpra_ptr x);

// Get from an arpra_t.
void arpra_get_bounds (mpfr_ptr lo, mpfr_ptr hi, arpra_srcptr x);

// Set to an arpra_t.
void arpra_set (arpra_ptr z, arpra_srcptr x);
void arpra_set_d (arpra_ptr z, const double centre);
void arpra_set_d_rad (arpra_ptr z, const double centre, const double radius);
void arpra_set_str (arpra_ptr z, const char *centre, arpra_int base);
void arpra_set_str_rad (arpra_ptr z, const char *centre, const char *radius, arpra_int base);
void arpra_set_mpfr (arpra_ptr z, mpfr_srcptr centre);
void arpra_set_mpfr_rad (arpra_ptr z, mpfr_srcptr centre, mpfr_srcptr radius);

// Set special values.
void arpra_set_nan (arpra_ptr z);
void arpra_set_inf (arpra_ptr z);
void arpra_set_zero (arpra_ptr z);

// Affine operations.
void arpra_affine_1 (arpra_ptr z, arpra_srcptr x,
                     mpfr_srcptr alpha, mpfr_srcptr gamma, mpfr_srcptr delta);
void arpra_affine_2 (arpra_ptr z, arpra_srcptr x, arpra_srcptr y,
                     mpfr_srcptr alpha, mpfr_srcptr beta, mpfr_srcptr gamma, mpfr_srcptr delta);
void arpra_add (arpra_ptr z, arpra_srcptr x, arpra_srcptr y);
void arpra_sub (arpra_ptr z, arpra_srcptr x, arpra_srcptr y);
void arpra_neg (arpra_ptr z, arpra_srcptr x);
void arpra_sum (arpra_ptr z, const arpra_ptr *x, arpra_uint n);

// Non-affine operations.
void arpra_mul(arpra_ptr z, arpra_srcptr x, arpra_srcptr y);
void arpra_div (arpra_ptr z, arpra_srcptr x, arpra_srcptr y);
void arpra_sqrt (arpra_ptr z, arpra_srcptr x);
void arpra_exp (arpra_ptr z, arpra_srcptr x);
void arpra_log (arpra_ptr z, arpra_srcptr x);
void arpra_inv (arpra_ptr z, arpra_srcptr x);

// Numerical precision.
arpra_precision arpra_get_precision (arpra_srcptr x);
arpra_precision arpra_get_default_precision ();
arpra_precision arpra_get_internal_precision ();
void arpra_set_precision (arpra_ptr z, arpra_precision precision);
void arpra_set_default_precision (arpra_precision precision);
void arpra_set_internal_precision (arpra_precision precision);

// Deviation term reduction.
void arpra_reduce_last_n (arpra_ptr z, arpra_uint n);
void arpra_reduce_small (arpra_ptr z, double fraction);

// Predicates on arpra_t.
int arpra_nan_p (arpra_srcptr x);
int arpra_inf_p (arpra_srcptr x);
int arpra_bounded_p (arpra_srcptr x);
int arpra_zero_p (arpra_srcptr x);
int arpra_has_zero_p (arpra_srcptr x);
int arpra_has_pos_p (arpra_srcptr x);
int arpra_has_neg_p (arpra_srcptr x);

// Get new deviation symbols.
arpra_uint arpra_next_sym ();

#ifdef __cplusplus
}
#endif

#endif // ARPRA_H
