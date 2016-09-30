/*
 * mpfa.h
 *
 *  Created on: 21 Jul 2016
 *      Author: jt273
 */

#ifndef MPFA_H
#define MPFA_H

#include <gmp.h>
#include <mpfr.h>

typedef struct {
	__mpfr_struct centre;
	__mpfr_struct radius;
} __mpfa_struct;

typedef __mpfa_struct mpfa_t[1];
typedef __mpfa_struct *mpfa_ptr;
typedef __gmp_const __mpfa_struct *mpfa_srcptr;

#ifdef __cplusplus
extern "C" {
#endif

// Initialise
void mpfa_init (mpfa_ptr x);
void mpfa_init2 (mpfa_ptr x, mp_prec_t p);

// Clear
void mpfa_clear (mpfa_ptr x);

// Get and set precision
mp_prec_t mpfa_get_prec (mpfa_srcptr x);
void mpfa_set_prec (mpfa_ptr x, mp_prec_t p);

// Set number
int mpfa_set (mpfa_ptr x, mpfa_srcptr newVal);
int mpfa_set_si (mpfa_ptr x, const long newVal);
int mpfa_set_ui (mpfa_ptr x, const long unsigned newVal);
int mpfa_set_d (mpfa_ptr x, const double newVal);
int mpfa_set_z (mpfa_ptr x, mpz_srcptr newVal);
int mpfa_set_q (mpfa_ptr x, mpq_srcptr newVal);
int mpfa_set_fr (mpfa_ptr x, mpfr_srcptr newVal);
int mpfa_set_str (mpfa_ptr x, const char *newVal, int base);

// Combined initialise and set
int mpfa_init_set (mpfa_ptr x, mpfa_srcptr newVal);
int mpfa_init_set_si (mpfa_ptr x, const long newVal);
int mpfa_init_set_ui (mpfa_ptr x, const long unsigned newVal);
int mpfa_init_set_d (mpfa_ptr x, const double newVal);
int mpfa_init_set_z (mpfa_ptr x, mpz_srcptr newVal);
int mpfa_init_set_q (mpfa_ptr x, mpq_srcptr newVal);
int mpfa_init_set_fr (mpfa_ptr x, mpfr_srcptr newVal);
int mpfa_init_set_str (mpfa_ptr x, const char *newVal, int base);

#ifdef __cplusplus
}
#endif

#endif /* MPFA_H */
