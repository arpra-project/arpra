/*
 * mpfa.h -- MPFA public header file.
 *
 * Copyright 2016-2017 James Paul Turner.
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

#ifndef MPFA_H
#define MPFA_H

#include <mpfr.h>

typedef int mpfa_int_t;
typedef unsigned int mpfa_uint_t;

typedef mpfr_prec_t mpfa_prec_t;
typedef mpfr_exp_t mpfa_exp_t;

struct __mpfa_struct {
    __mpfr_struct centre;
    __mpfr_struct radius;
    __mpfr_struct *deviations;
    mpfa_uint_t *symbols;
    mpfa_uint_t nTerms;
};

typedef struct __mpfa_struct mpfa_t[1];
typedef struct __mpfa_struct *mpfa_ptr;
typedef __gmp_const struct __mpfa_struct *mpfa_srcptr;

#ifdef __cplusplus
extern "C" {
#endif

// Initialise and clear
void mpfa_init (mpfa_ptr x);
void mpfa_inits (mpfa_ptr x, ...);
void mpfa_init2 (mpfa_ptr x, mpfa_prec_t prec);
void mpfa_inits2 (mpfa_prec_t prec, mpfa_ptr x, ...);
void mpfa_clear (mpfa_ptr x);
void mpfa_clears (mpfa_ptr x, ...);

// Set affine form
void mpfa_set (mpfa_ptr z, mpfa_srcptr x);
void mpfa_set_d (mpfa_ptr z, const double centre);
void mpfa_set_d_rad (mpfa_ptr z, const double centre, const double radius);
void mpfa_set_str (mpfa_ptr z, const char *centre, mpfa_int_t base);
void mpfa_set_str_rad (mpfa_ptr z, const char *centre, const char *radius, mpfa_int_t base);
void mpfa_set_mpfr (mpfa_ptr z, mpfr_srcptr centre);
void mpfa_set_mpfr_rad (mpfa_ptr z, mpfr_srcptr centre, mpfr_srcptr radius);

// Affine operations
void mpfa_affine_1 (mpfa_ptr z, mpfa_srcptr x, mpfr_srcptr alpha, mpfr_srcptr gamma, mpfr_srcptr delta);
void mpfa_affine_2 (mpfa_ptr z, mpfa_srcptr x, mpfa_srcptr y, mpfr_srcptr alpha, mpfr_srcptr beta, mpfr_srcptr gamma, mpfr_srcptr delta);
void mpfa_add (mpfa_ptr z, mpfa_srcptr x, mpfa_srcptr y);
void mpfa_sub (mpfa_ptr z, mpfa_srcptr x, mpfa_srcptr y);
void mpfa_neg (mpfa_ptr z, mpfa_srcptr x);
void mpfa_sum (mpfa_ptr z, const mpfa_ptr *x, mpfa_uint_t n);

// Non-affine operations
void mpfa_mul(mpfa_ptr z, mpfa_srcptr x, mpfa_srcptr y);
void mpfa_div (mpfa_ptr z, mpfa_srcptr x, mpfa_srcptr y);
void mpfa_sqrt (mpfa_ptr z, mpfa_srcptr x);
void mpfa_exp (mpfa_ptr z, mpfa_srcptr x);
void mpfa_log (mpfa_ptr z, mpfa_srcptr x);
void mpfa_inv (mpfa_ptr z, mpfa_srcptr x);

// Condense functions.
void mpfa_condense_last_n (mpfa_ptr z, mpfa_uint_t n);
void mpfa_condense_small (mpfa_ptr z, double fraction);

// Get and set precision
mpfa_prec_t mpfa_get_prec (mpfa_srcptr x);
void mpfa_set_prec (mpfa_ptr x, mpfa_prec_t prec);
mpfa_prec_t mpfa_get_default_prec ();
void mpfa_set_default_prec (mpfa_prec_t prec);
mpfa_prec_t mpfa_get_internal_prec ();
void mpfa_set_internal_prec (mpfa_prec_t prec);

// Helper functions
mpfa_uint_t mpfa_next_sym();
mpfa_int_t mpfa_term (mpfr_ptr z, mpfr_srcptr x, mpfr_srcptr y, mpfr_srcptr alpha, mpfr_srcptr beta, mpfr_srcptr gamma);
void mpfa_error (mpfr_ptr error, mpfr_srcptr x);

#ifdef __cplusplus
}
#endif

#endif // MPFA_H
