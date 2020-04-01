/*
 * arpra-impl.h -- Private header file, used during compilation.
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

#ifndef ARPRA_IMPL_H
#define ARPRA_IMPL_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif // HAVE_CONFIG_H

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>

#include <arpra.h>
#include <arpra_ode.h>

// Default range analysis method.
#define ARPRA_DEFAULT_RANGE_METHOD ARPRA_MIXED_TRIMMED_IAAA

// Default multiplication method.
#define ARPRA_DEFAULT_MUL_METHOD ARPRA_MUL_RUMP_KASHIWAGI

// Default precisions.
#define ARPRA_DEFAULT_PRECISION 53
#define ARPRA_DEFAULT_INTERNAL_PRECISION 256

// Min-Range approximation.
//#define ARPRA_MIN_RANGE 1

// Temp buffers.
#define ARPRA_BUFFER_RESIZE_FACTOR 256

// Internal auxiliary functions.
void arpra_helper_affine_1 (arpra_range *y, const arpra_range *x1,
                            mpfr_srcptr alpha, mpfr_srcptr gamma, mpfr_srcptr delta);
void arpra_helper_affine_2 (arpra_range *y, const arpra_range *x1, const arpra_range *x2,
                             mpfr_srcptr alpha, mpfr_srcptr beta, mpfr_srcptr gamma, mpfr_srcptr delta);
void arpra_helper_mpfr_rnderr (mpfr_ptr err, mpfr_rnd_t rnd, mpfr_srcptr y);
void arpra_helper_compute_range (arpra_range *y);
void arpra_helper_mix_trim (arpra_range *y, mpfi_srcptr ia_range);
void arpra_helper_check_result (arpra_range *y);
void arpra_helper_set_symbol_count (arpra_uint n);
arpra_uint arpra_helper_get_symbol_count ();
arpra_uint arpra_helper_next_symbol ();
mpfr_ptr *arpra_helper_buffer_mpfr_ptr (arpra_uint n);
mpfr_ptr arpra_helper_buffer_mpfr (arpra_uint n);
void arpra_helper_clear_terms (arpra_range *y);

// Arpra extensions to the MPFR library.
int arpra_ext_mpfr_fmma (mpfr_ptr y, mpfr_srcptr x1, mpfr_srcptr x2,
                         mpfr_srcptr x3, mpfr_srcptr x4, mpfr_rnd_t rnd);
int arpra_ext_mpfr_fmmaa (mpfr_ptr y, mpfr_srcptr x1, mpfr_srcptr x2,
                          mpfr_srcptr x3, mpfr_srcptr x4, mpfr_srcptr x5, mpfr_rnd_t rnd);
int arpra_ext_mpfr_sum (mpfr_ptr y, mpfr_ptr x,
                        const arpra_uint n, const mpfr_rnd_t rnd);
int arpra_ext_mpfr_sumabs (mpfr_ptr y, mpfr_ptr x,
                           const arpra_uint n, const mpfr_rnd_t rnd);

// arpra_helper_mpfr_rnderr function wrapper macros.
#define ARPRA_MPFR_RNDERR(err, rnd, fn, y, ...)                         \
    if (fn(y, __VA_ARGS__, rnd)) arpra_helper_mpfr_rnderr(err, rnd, y)

#define ARPRA_MPFR_RNDERR_SET(err, rnd, y, x1)                          \
    if (mpfr_set(y, x1, rnd)) arpra_helper_mpfr_rnderr(err, rnd, y)

#define ARPRA_MPFR_RNDERR_ADD(err, rnd, y, x1, x2)                      \
    if (mpfr_add(y, x1, x2, rnd)) arpra_helper_mpfr_rnderr(err, rnd, y)

#define ARPRA_MPFR_RNDERR_SUB(err, rnd, y, x1, x2)                      \
    if (mpfr_sub(y, x1, x2, rnd)) arpra_helper_mpfr_rnderr(err, rnd, y)

#define ARPRA_MPFR_RNDERR_MUL(err, rnd, y, x1, x2)                      \
    if (mpfr_mul(y, x1, x2, rnd)) arpra_helper_mpfr_rnderr(err, rnd, y)

#define ARPRA_MPFR_RNDERR_FMA(err, rnd, y, x1, x2, x3)                  \
    if (mpfr_fma(y, x1, x2, x3, rnd)) arpra_helper_mpfr_rnderr(err, rnd, y)

#define ARPRA_MPFR_RNDERR_FMMA(err, rnd, y, x1, x2, x3, x4)             \
    if (arpra_ext_mpfr_fmma(y, x1, x2, x3, x4, rnd)) arpra_helper_mpfr_rnderr(err, rnd, y)

#define ARPRA_MPFR_RNDERR_FMMAA(err, rnd, y, x1, x2, x3, x4, x5)        \
    if (arpra_ext_mpfr_fmmaa(y, x1, x2, x3, x4, x5, rnd)) arpra_helper_mpfr_rnderr(err, rnd, y)

#define ARPRA_MPFR_RNDERR_SUM(err, rnd, y, x, n)                        \
    if (mpfr_sum(y, x, n, rnd)) arpra_helper_mpfr_rnderr(err, rnd, y)

#endif // ARPRA_IMPL_H
