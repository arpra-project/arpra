/*
 * arpra.h -- Arpra main public header.
 *
 * Copyright 2016-2020 James Paul Turner.
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
#include <mpfi.h>

// Arpra type definitions.
typedef long int arpra_int;
typedef unsigned long int arpra_uint;
typedef mpfr_prec_t arpra_prec;
typedef __mpfr_struct arpra_mpfr;
typedef __mpfi_struct arpra_mpfi;
typedef struct arpra_range_struct arpra_range;

// The Arpra range struct.
struct arpra_range_struct
{
    arpra_prec precision;
    __mpfr_struct centre;
    __mpfr_struct radius;
    __mpfi_struct true_range;
    __mpfr_struct *deviations;
    arpra_uint *symbols;
    arpra_uint nTerms;
};

#ifdef __cplusplus
extern "C" {
#endif

// Initialise and clear.
void arpra_init (arpra_range *z);
void arpra_init2 (arpra_range *z, const arpra_prec prec);
void arpra_clear (arpra_range *z);

// Get from an Arpra range.
void arpra_get_bounds (arpra_mpfr *lo, arpra_mpfr *hi, const arpra_range *x);
void arpra_get_mpfi (arpra_mpfi *z, const arpra_range *x);

// Set to an Arpra range.
void arpra_set_range (arpra_range *z, const arpra_range *x);
void arpra_set_d (arpra_range *z, const double centre);
void arpra_set_d_rad (arpra_range *z, const double centre, const double radius);
void arpra_set_str (arpra_range *z, const char *centre, const arpra_int base);
void arpra_set_str_rad (arpra_range *z, const char *centre, const char *radius, const arpra_int base);
void arpra_set_mpfr (arpra_range *z, const arpra_mpfr *centre);
void arpra_set_mpfr_rad (arpra_range *z, const arpra_mpfr *centre, const arpra_mpfr *radius);
void arpra_set_mpfi (arpra_range *z, const arpra_mpfi *x);

// Set special values.
void arpra_set_nan (arpra_range *z);
void arpra_set_inf (arpra_range *z);
void arpra_set_zero (arpra_range *z);

// Affine operations.
void arpra_affine_1 (arpra_range *y, const arpra_range *x1,
                     mpfr_srcptr alpha, mpfr_srcptr gamma, mpfr_srcptr delta);
void arpra_affine_2 (arpra_range *y, const arpra_range *x1, const arpra_range *x2,
                     mpfr_srcptr alpha, mpfr_srcptr beta, mpfr_srcptr gamma, mpfr_srcptr delta);
void arpra_add (arpra_range *z, const arpra_range *x, const arpra_range *y);
void arpra_sub (arpra_range *z, const arpra_range *x, const arpra_range *y);
void arpra_neg (arpra_range *z, const arpra_range *x);
void arpra_sum_exact (arpra_range *z, arpra_range *x, const arpra_uint n);
void arpra_sum_recursive (arpra_range *z, arpra_range *x, const arpra_uint n);

// Non-affine operations.
void arpra_mul (arpra_range *z, const arpra_range *x, const arpra_range *y);
void arpra_div (arpra_range *z, const arpra_range *x, const arpra_range *y);
void arpra_sqrt (arpra_range *z, const arpra_range *x);
void arpra_exp (arpra_range *z, const arpra_range *x);
void arpra_log (arpra_range *z, const arpra_range *x);
void arpra_inv (arpra_range *z, const arpra_range *x);

// Numerical precision.
arpra_prec arpra_get_precision (const arpra_range *x);
arpra_prec arpra_get_default_precision ();
arpra_prec arpra_get_internal_precision ();
void arpra_set_precision (arpra_range *z, const arpra_prec prec);
void arpra_set_default_precision (const arpra_prec prec);
void arpra_set_internal_precision (const arpra_prec prec);

// Deviation term reduction.
void arpra_reduce_last_n (arpra_range *z, arpra_uint n);
void arpra_reduce_small (arpra_range *z, double min_fraction);
void arpra_reduce_small_abs (arpra_range *z, const arpra_mpfr *abs_threshold);
void arpra_reduce_small_rel (arpra_range *z, const arpra_mpfr *rel_threshold);

// Predicates on Arpra ranges.
int arpra_nan_p (const arpra_range *x);
int arpra_inf_p (const arpra_range *x);
int arpra_bounded_p (const arpra_range *x);
int arpra_zero_p (const arpra_range *x);
int arpra_has_zero_p (const arpra_range *x);
int arpra_has_pos_p (const arpra_range *x);
int arpra_has_neg_p (const arpra_range *x);

// Clear temporary data.
void arpra_clear_buffers ();

// MPFR wrapper functions.
void arpra_mpfr_fn1 (int (*fn) (mpfr_ptr y, mpfr_srcptr x1, mpfr_rnd_t),
                     arpra_range *y, mpfr_srcptr x1);
void arpra_mpfr_fn1_ui (int (*fn) (mpfr_ptr y, unsigned long int x1, mpfr_rnd_t),
                        arpra_range *y, unsigned long int x1);
void arpra_mpfr_fn1_si (int (*fn) (mpfr_ptr y, long int x1, mpfr_rnd_t),
                        arpra_range *y, long int x1);
void arpra_mpfr_fn1_d (int (*fn) (mpfr_ptr y, double x1, mpfr_rnd_t),
                       arpra_range *y, double x1);
void arpra_mpfr_fn2 (int (*fn) (mpfr_ptr y, mpfr_srcptr x1, mpfr_srcptr x2, mpfr_rnd_t),
                     arpra_range *y, mpfr_srcptr x1, mpfr_srcptr x2);
void arpra_mpfr_ui_fn2 (int (*fn) (mpfr_ptr y, unsigned long int x1, mpfr_srcptr x2, mpfr_rnd_t),
                        arpra_range *y, unsigned long int x1, mpfr_srcptr x2);
void arpra_mpfr_fn2_ui (int (*fn) (mpfr_ptr y, mpfr_srcptr x1, unsigned long int x2, mpfr_rnd_t),
                        arpra_range *y, mpfr_srcptr x1, unsigned long int x2);
void arpra_mpfr_si_fn2 (int (*fn) (mpfr_ptr y, long int x1, mpfr_srcptr x2, mpfr_rnd_t),
                        arpra_range *y, long int x1, mpfr_srcptr x2);
void arpra_mpfr_fn2_si (int (*fn) (mpfr_ptr y, mpfr_srcptr x1, long int x2, mpfr_rnd_t),
                        arpra_range *y, mpfr_srcptr x1, long int x2);
void arpra_mpfr_d_fn2 (int (*fn) (mpfr_ptr y, double x1, mpfr_srcptr x2, mpfr_rnd_t),
                       arpra_range *y, double x1, mpfr_srcptr x2);
void arpra_mpfr_fn2_d (int (*fn) (mpfr_ptr y, mpfr_srcptr x1, double x2, mpfr_rnd_t),
                       arpra_range *y, mpfr_srcptr x1, double x2);
void arpra_mpfr_set_str (arpra_range *y, const char *x1, int base);

#ifdef __cplusplus
}
#endif

// MPFR setter wrappers.
#define arpra_set_mpfr(y, x1) arpra_mpfr_fn1(mpfr_set, y, x1)
#define arpra_set_ui(y, x1) arpra_mpfr_fn1_ui(mpfr_set_ui, y, x1)
#define arpra_set_si(y, x1) arpra_mpfr_fn1_si(mpfr_set_si, y, x1)
#define arpra_set_d(y, x1) arpra_mpfr_fn1_d(mpfr_set_d, y, x1)

#endif // ARPRA_H
