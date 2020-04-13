/*
 * morris_lecar_mpfi.c -- Test Morris-Lecar model using MPFI.
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

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <mpfi.h>

/*
 * Global variables
 * ----------------
 * h      Step size (msec)
 * t      Current time (msec)
 *
 * Neuron variables
 * ----------------
 * N      Fraction of open K+ channels
 * V      Membrane potential (mV)
 *
 * Neuron parameters
 * -----------------
 * GL     Maximum leak conductance (mS/cm^2)
 * VL     Equilibrium potential of leak conductance (mV)
 * GCa    Maximum Ca++ conductance (mS/cm^2)
 * VCa    Equilibrium potential of Ca++ conductance (mV)
 * GK     Maximum K+ conductance (mS/cm^2)
 * VK     Equilibrium potential of K+ conductance (mV)
 * V1     Potential at which M_ss(V) = 0.5 (mV)
 * V2     Reciprocal of voltage dependence slope of M_ss(V) (mV)
 * V3     Potential at which N_ss(V) = 0.5 (mV)
 * V4     Reciprocal of voltage dependence slope of N_ss(V) (mV)
 * phi    (s^-1)
 * C      Membrane capacitance (uF/cm^2)
 *
 * Synapse variables
 * -----------------
 * R      Neurotransmitter release
 * S      Neurotransmitter binding
 *
 * Synapse parameters
 * ------------------
 * GSyn   Maximum synapse conductance (mS/cm^2)
 * VSyn   Synapse conductance equilibrium potential (mV)
 * thr    Neuronal spike threshold (mV)
 * a      Transmitter release/bind rise factor
 * b      Transmitter release/bind decay factor
 * k      Steepness of activation function
 */

// General parameters
#define p_h 0.5
#define p_t0 0.0
#define p_prec 53
#define p_error_prec 256
#define p_sim_steps 1000
#define p_report_step 20

// RNG parameters
// Seeds are random if not #defined
#define p_rand_prec 53
//#define p_rng_uf_seed 707135875931353ul
//#define p_rng_nf_seed 503108552855933ul

// Poisson input parameters
#define p_in_size 50
#define p_in_freq 42.0
#define p_in_V_lo -60.0
#define p_in_V_hi 20.0

// Neuron parameters (class 1)
#define p_nrn_size 1
#define p_nrn_N0 0.0
#define p_nrn_V0 -60.0
#define p_nrn_GL 2.0
#define p_nrn_GCa 4.0
#define p_nrn_GK 8.0
#define p_nrn_VL -60.0
#define p_nrn_VCa 120.0
#define p_nrn_VK -80.0
#define p_nrn_V1 -1.2
#define p_nrn_V2 18.0
#define p_nrn_V3 12.0
#define p_nrn_V4 17.4
#define p_nrn_phi 1.0 / 15.0
#define p_nrn_C 20.0

/* // Neuron parameters (class 2) */
/* #define p_nrn_size 1 */
/* #define p_nrn_N0 0.0 */
/* #define p_nrn_V0 -60.0 */
/* #define p_nrn_GL 2.0 */
/* #define p_nrn_GCa 4.4 */
/* #define p_nrn_GK 8.0 */
/* #define p_nrn_VL -60.0 */
/* #define p_nrn_VCa 120.0 */
/* #define p_nrn_VK -80.0 */
/* #define p_nrn_V1 -1.2 */
/* #define p_nrn_V2 18.0 */
/* #define p_nrn_V3 2.0 */
/* #define p_nrn_V4 30.0 */
/* #define p_nrn_phi 1.0 / 25.0 */
/* #define p_nrn_C 20.0 */

// Synapse parameters (excitatory)
#define p_syn_size p_in_size * p_nrn_size
#define p_syn_R0 0.0
#define p_syn_S0 0.0
#define p_syn_GSyn_std 0.05
#define p_syn_GSyn_mean 25.0 / p_in_size
#define p_syn_VSyn 0.0
#define p_syn_thr -50.0
#define p_syn_a 0.25 // in [1/10, 1/2]
#define p_syn_b 0.15 // in [1/20, 1/4]
#define p_syn_k 1.0E6

/* // Synapse parameters (inhibitory) */
/* #define p_syn_size p_in_size * p_nrn_size */
/* #define p_syn_R0 0.0 */
/* #define p_syn_S0 0.0 */
/* #define p_syn_GSyn_std 0.5 */
/* #define p_syn_GSyn_mean 3.0 */
/* #define p_syn_VSyn -80.0 */
/* #define p_syn_thr -50.0 */
/* #define p_syn_a 0.075 // in [1/20, 1/10] */
/* #define p_syn_b 0.035 // in [1/50, 1/20] */
/* #define p_syn_k 1.0E6 */

// ===================== end of model parameters ======================


int *in;
mpfr_t rand_uf, rand_nf, temp_error, in_p0;
mpfr_ptr temp_sum_error, *temp_sum_error_ptr;
mpfi_t nrn_GL, nrn_VL, nrn_GCa, nrn_VCa, nrn_GK, nrn_VK, nrn_V1, nrn_V2, nrn_V3, nrn_V4,
    nrn_phi, nrn_C, syn_VSyn, syn_thr, syn_a, syn_b, syn_k, one, two, neg_two, temp_sum,
    temp1, temp2, M_ss, N_ss, in_V_lo, in_V_hi;
mpfi_ptr syn_GSyn, I;
gmp_randstate_t rng_uf, rng_nf;

// System state variables
mpfi_ptr nrn_N, nrn_V, syn_R, syn_S;
mpfi_ptr d_nrn_N, d_nrn_V, d_syn_R, d_syn_S;

// DEBUG: print MPFI numbers to stderr
void debug_i (mpfi_srcptr x) {
    fputs("[", stderr);
    mpfr_out_str(stderr, 10, 40, &(x->left), MPFR_RNDN);
    fputs(", ", stderr);
    mpfr_out_str(stderr, 10, 40, &(x->right), MPFR_RNDN);
    fputs("]\n", stderr);
}

// DEBUG: print MPFR numbers to stderr
void debug_r (mpfr_srcptr x) {
    mpfr_out_str(stderr, 10, 60, x, MPFR_RNDN);
    fputs("\n", stderr);
}

void file_init (char *grp, unsigned long grp_size, FILE **f)
{
    char fname[20];
    unsigned long i;

    for (i = 0; i < grp_size; i++) {
        sprintf(fname, "%s_%03lu.dat", grp, i);
        f[i] = fopen(fname, "w");
    }
}

void file_clear (unsigned long grp_size, FILE **f)
{
    unsigned long i;

    for (i = 0; i < grp_size; i++) {
        fclose(f[i]);
    }
}

void file_write (mpfi_srcptr A, unsigned long grp_size, FILE **f)
{
    unsigned long i;

    for (i = 0; i < grp_size; i++) {
        mpfr_out_str(f[i], 10, 40, &(A[i].left), MPFR_RNDN);
        fputs(" ", f[i]);
        mpfr_out_str(f[i], 10, 40, &(A[i].right), MPFR_RNDN);
        fputs("\n", f[i]);
    }
}

void dNdt (const unsigned long idx)
{
    mpfi_ptr d_N = &(d_nrn_N[idx]);
    mpfi_srcptr N = &(nrn_N[idx]);
    mpfi_srcptr V = &(nrn_V[idx]);
    mpfi_srcptr V3 = nrn_V3;
    mpfi_srcptr V4 = nrn_V4;
    mpfi_srcptr phi = nrn_phi;

    // K+ channel activation steady-state
    // N_ss = 1 / (1 + exp(-2 (V - V3) / V4))
    mpfi_sub(temp1, V, V3);
    mpfi_mul(N_ss, neg_two, temp1);
    mpfi_div(N_ss, N_ss, V4);
    mpfi_exp(N_ss, N_ss);
    mpfi_add(N_ss, one, N_ss);
    mpfi_ui_div(N_ss, 1, N_ss);

    // tau of K+ channel activation
    // tau = 1 / (phi ((p + q) / 2))
    // p = exp(-(V - V3) / (2 V4))
    // q = exp( (V - V3) / (2 V4))
    mpfi_mul(temp2, two, V4);
    mpfi_div(temp2, temp1, temp2);
    mpfi_neg(temp1, temp2);
    mpfi_exp(temp1, temp1);
    mpfi_exp(temp2, temp2);
    mpfi_add(temp1, temp1, temp2);
    mpfi_div(temp1, temp1, two);
    mpfi_mul(temp1, phi, temp1);
    mpfi_ui_div(temp1, 1, temp1);

    // delta of K+ channel activation
    // dN/dt = (N_ss - N) / tau
    mpfi_sub(d_N, N_ss, N);
    mpfi_div(d_N, d_N, temp1);
}

void dVdt (const unsigned long idx)
{
    unsigned long i;
    unsigned long pre_size = p_in_size;
    mpfi_ptr d_V = &(d_nrn_V[idx]);
    mpfi_srcptr N = &(nrn_N[idx]);
    mpfi_srcptr V = &(nrn_V[idx]);
    mpfi_srcptr S = &(syn_S[idx * pre_size]);
    mpfi_srcptr GSyn = &(syn_GSyn[idx * pre_size]);
    mpfi_srcptr VSyn = syn_VSyn;
    mpfi_srcptr VL = nrn_VL;
    mpfi_srcptr GL = nrn_GL;
    mpfi_srcptr VCa = nrn_VCa;
    mpfi_srcptr GCa = nrn_GCa;
    mpfi_srcptr VK = nrn_VK;
    mpfi_srcptr GK = nrn_GK;
    mpfi_srcptr V1 = nrn_V1;
    mpfi_srcptr V2 = nrn_V2;
    mpfi_srcptr C = nrn_C;

    // Ca++ channel activation steady-state
    // M_ss = 1 / (1 + exp(-2 (V - V1) / V2))
    mpfi_sub(M_ss, V, V1);
    mpfi_mul(M_ss, neg_two, M_ss);
    mpfi_div(M_ss, M_ss, V2);
    mpfi_exp(M_ss, M_ss);
    mpfi_add(M_ss, one, M_ss);
    mpfi_ui_div(M_ss, 1, M_ss);

    // Synapse current
    mpfi_sub(temp1, VSyn, V);
    mpfi_set_ui(temp_sum, 0);
    for (i = 0; i < pre_size; i++) {
        mpfi_mul(&(I[i]), temp1, &(GSyn[i]));
        mpfi_mul(&(I[i]), &(I[i]), &(S[i]));
        mpfi_add(temp_sum, temp_sum, &(I[i]));
        if (mpfr_cmpabs(&(I[i].left), &(I[i].right)) <= 0) {
            mpfr_set(&(temp_sum_error[i]), &(I[i].right), MPFR_RNDU);
        }
        else {
            mpfr_abs(&(temp_sum_error[i]), &(I[i].left), MPFR_RNDU);
        }
    }

    // Compute error(sum(I)) = (n - 1) u sum(|I|)
    mpfr_set_si_2exp(temp_error, 1, -p_prec, MPFR_RNDU);
    mpfr_mul_ui(temp_error, temp_error, (pre_size - 1), MPFR_RNDU);
    mpfr_sum(temp_sum_error_ptr[0], temp_sum_error_ptr, pre_size, MPFR_RNDU);
    mpfr_mul(temp_error, temp_error, temp_sum_error_ptr[0], MPFR_RNDU);
    mpfr_sub(&(temp_sum->left), &(temp_sum->left), temp_error, MPFR_RNDD);
    mpfr_add(&(temp_sum->right), &(temp_sum->right), temp_error, MPFR_RNDU);
    mpfi_set(d_V, temp_sum);



    // ======== TEMP DEBUG ==========
    //mpfi_set_d(d_V, 80.0); // bifurcation at sum(I) = 80.0
    if (idx == 0) {
        //fprintf(stderr, "sum(I): "); debug_i(d_V);
        //fprintf(stderr, "I[0]: "); debug_i(I);
        for (i = 0; i < pre_size; i++) {
            //fprintf(stderr, "I[%lu]: ", i); debug_i(&(I[i]));
        }
    }



    // Leak current
    mpfi_sub(temp1, V, VL);
    mpfi_mul(temp1, temp1, GL);
    mpfi_sub(d_V, d_V, temp1);

    // Ca++ current
    mpfi_sub(temp1, V, VCa);
    mpfi_mul(temp1, temp1, GCa);
    mpfi_mul(temp1, temp1, M_ss);
    mpfi_sub(d_V, d_V, temp1);

    // K+ current
    mpfi_sub(temp1, V, VK);
    mpfi_mul(temp1, temp1, GK);
    mpfi_mul(temp1, temp1, N);
    mpfi_sub(d_V, d_V, temp1);

    // delta of membrane potential
    // dV/dt = (I + GL (VL - V) + GCa M (VCa - V) + GK N (VK - V)) / C
    mpfi_div(d_V, d_V, C);
}

void dRdt (const unsigned long idx)
{
    mpfi_ptr d_R = &(d_syn_R[idx]);
    mpfi_srcptr R = &(syn_R[idx]);
    mpfi_srcptr a = syn_a;
    mpfi_srcptr b = syn_b;
    mpfi_srcptr k = syn_k;
    mpfi_srcptr VPre = in[idx % p_in_size] ? in_V_hi : in_V_lo;
    mpfi_srcptr threshold = syn_thr;

    // Sigmoid of threshold difference
    mpfi_sub(temp1, VPre, threshold);
    mpfi_mul(temp1, temp1, k);
    mpfi_exp(temp1, temp1);
    mpfi_add(temp1, temp1, one);
    mpfi_ui_div(temp1, 1, temp1);

    // Presynaptic transmitter release rise
    mpfi_mul(temp1, a, temp1);

    // Presynaptic transmitter release decay
    mpfi_mul(temp2, b, R);

    // delta of presynaptic transmitter release
    // dR/dt = a Q - b R
    // Q = 1 / (1 + e^(k(V - threshold)))
    mpfi_sub(d_R, temp1, temp2);
}

void dSdt (const unsigned long idx)
{
    mpfi_ptr d_S = &(d_syn_S[idx]);
    mpfi_srcptr R = &(syn_R[idx]);
    mpfi_srcptr S = &(syn_S[idx]);
    mpfi_srcptr a = syn_a;
    mpfi_srcptr b = syn_b;

    // Postsynaptic transmitter binding rise
    mpfi_mul(temp1, a, R);

    // Postsynaptic transmitter binding decay
    mpfi_mul(temp2, b, S);

    // delta of postsynaptic transmitter binding
    // dS/dt = a R - b S
    mpfi_sub(d_S, temp1, temp2);
}

int main (int argc, char *argv[])
{
    mpfi_t h, t;
    unsigned long i, j;

    // Allocate system state
    nrn_N = malloc(p_nrn_size * sizeof(mpfi_t));
    nrn_V = malloc(p_nrn_size * sizeof(mpfi_t));
    syn_R = malloc(p_syn_size * sizeof(mpfi_t));
    syn_S = malloc(p_syn_size * sizeof(mpfi_t));

    // Allocate other arrays
    d_nrn_N = malloc(p_nrn_size * sizeof(mpfi_t));
    d_nrn_V = malloc(p_nrn_size * sizeof(mpfi_t));
    d_syn_R = malloc(p_syn_size * sizeof(mpfi_t));
    d_syn_S = malloc(p_syn_size * sizeof(mpfi_t));
    syn_GSyn = malloc(p_syn_size * sizeof(mpfi_t));
    temp_sum_error = malloc(p_in_size * sizeof(mpfr_t));
    temp_sum_error_ptr = malloc(p_in_size * sizeof(mpfr_ptr));
    I = malloc(p_in_size * sizeof(mpfi_t));
    in = malloc(p_in_size * sizeof(int));

    struct timespec clock_time;

    // Initialise uniform float RNG
    gmp_randinit_default(rng_uf);
    clock_gettime(CLOCK_REALTIME, &clock_time);
#ifdef p_rng_uf_seed
    unsigned long rng_uf_seed = p_rng_uf_seed;
#else
    unsigned long rng_uf_seed = clock_time.tv_sec + clock_time.tv_nsec;
#endif
    gmp_randseed_ui(rng_uf, rng_uf_seed);
    printf("GMP rand uniform float seed: %lu\n", rng_uf_seed);

    // Initialise normal float RNG
    gmp_randinit_default(rng_nf);
    clock_gettime(CLOCK_REALTIME, &clock_time);
#ifdef p_rng_nf_seed
    unsigned long rng_nf_seed = p_rng_nf_seed;
#else
    unsigned long rng_nf_seed = clock_time.tv_sec + clock_time.tv_nsec;
#endif
    gmp_randseed_ui(rng_nf, rng_nf_seed);
    printf("GMP rand normal float seed: %lu\n", rng_nf_seed);

    // Initialise system state
    mpfi_init2(h, p_prec);
    mpfi_init2(t, p_prec);
    for (i = 0; i < p_nrn_size; i++) {
        mpfi_init2(&(nrn_N[i]), p_prec);
        mpfi_init2(&(nrn_V[i]), p_prec);
        mpfi_init2(&(d_nrn_N[i]), p_prec);
        mpfi_init2(&(d_nrn_V[i]), p_prec);
    }
    for (i = 0; i < p_syn_size; i++) {
        mpfi_init2(&(syn_R[i]), p_prec);
        mpfi_init2(&(syn_S[i]), p_prec);
        mpfi_init2(&(d_syn_R[i]), p_prec);
        mpfi_init2(&(d_syn_S[i]), p_prec);
    }

    // Initialise Poisson input parameters
    mpfr_init2(in_p0, p_prec);
    mpfi_init2(in_V_lo, p_prec);
    mpfi_init2(in_V_hi, p_prec);

    // Initialise neuron parameters
    mpfi_init2(nrn_GL, p_prec);
    mpfi_init2(nrn_GCa, p_prec);
    mpfi_init2(nrn_GK, p_prec);
    mpfi_init2(nrn_VL, p_prec);
    mpfi_init2(nrn_VCa, p_prec);
    mpfi_init2(nrn_VK, p_prec);
    mpfi_init2(nrn_V1, p_prec);
    mpfi_init2(nrn_V2, p_prec);
    mpfi_init2(nrn_V3, p_prec);
    mpfi_init2(nrn_V4, p_prec);
    mpfi_init2(nrn_phi, p_prec);
    mpfi_init2(nrn_C, p_prec);

    // Initialise synapse parameters
    for (i = 0; i < p_syn_size; i++) {
        mpfi_init2(&(syn_GSyn[i]), p_prec);
    }
    mpfi_init2(syn_VSyn, p_prec);
    mpfi_init2(syn_thr, p_prec);
    mpfi_init2(syn_a, p_prec);
    mpfi_init2(syn_b, p_prec);
    mpfi_init2(syn_k, p_prec);

    // Initialise constants
    mpfi_init2(one, p_prec);
    mpfi_init2(two, p_prec);
    mpfi_init2(neg_two, p_prec);

    // Initialise scratch space
    mpfr_init2(rand_uf, p_rand_prec);
    mpfr_init2(rand_nf, p_rand_prec);
    mpfr_init2(temp_error, p_error_prec);
    mpfi_init2(temp_sum, p_error_prec);
    mpfi_init2(temp1, p_prec);
    mpfi_init2(temp2, p_prec);
    mpfi_init2(M_ss, p_prec);
    mpfi_init2(N_ss, p_prec);

    // Set system state
    mpfi_set_d(h, p_h);
    mpfi_set_d(t, p_t0);
    for (i = 0; i < p_nrn_size; i++) {
        mpfi_set_d(&(nrn_N[i]), p_nrn_N0);
        mpfi_set_d(&(nrn_V[i]), p_nrn_V0);
    }
    for (i = 0; i < p_syn_size; i++) {
        mpfi_set_d(&(syn_R[i]), p_syn_R0);
        mpfi_set_d(&(syn_S[i]), p_syn_S0);
    }

    // Set Poisson input parameters
    mpfr_set_d(in_p0, -(p_in_freq / 1000) * p_h, MPFR_RNDN);
    mpfr_exp(in_p0, in_p0, MPFR_RNDN);
    mpfi_set_d(in_V_lo, p_in_V_lo);
    mpfi_set_d(in_V_hi, p_in_V_hi);
    for (i = 0; i < p_in_size; i++) {
        mpfr_init2(&(temp_sum_error[i]), p_error_prec);
        temp_sum_error_ptr[i] = &(temp_sum_error[i]);
        mpfi_init2(&(I[i]), p_prec);
    }

    // Set neuron parameters
    mpfi_set_d(nrn_GL, p_nrn_GL);
    mpfi_set_d(nrn_GCa, p_nrn_GCa);
    mpfi_set_d(nrn_GK, p_nrn_GK);
    mpfi_set_d(nrn_VL, p_nrn_VL);
    mpfi_set_d(nrn_VCa, p_nrn_VCa);
    mpfi_set_d(nrn_VK, p_nrn_VK);
    mpfi_set_d(nrn_V1, p_nrn_V1);
    mpfi_set_d(nrn_V2, p_nrn_V2);
    mpfi_set_d(nrn_V3, p_nrn_V3);
    mpfi_set_d(nrn_V4, p_nrn_V4);
    mpfi_set_d(nrn_phi, p_nrn_phi);
    mpfi_set_d(nrn_C, p_nrn_C);

    // Set synapse parameters
    for (i = 0; i < p_syn_size; i++) {
        //mpfr_nrandom(rand_nf, rng_nf, MPFR_RNDN);
        mpfr_grandom(rand_nf, NULL, rng_nf, MPFR_RNDN);
        mpfr_mul_d(rand_nf, rand_nf, p_syn_GSyn_std, MPFR_RNDN);
        mpfr_add_d(rand_nf, rand_nf, p_syn_GSyn_mean, MPFR_RNDN);
        mpfi_set_fr(&(syn_GSyn[i]), rand_nf);
    }
    mpfi_set_d(syn_VSyn, p_syn_VSyn);
    mpfi_set_d(syn_thr, p_syn_thr);
    mpfi_set_d(syn_a, p_syn_a);
    mpfi_set_d(syn_b, p_syn_b);
    mpfi_set_d(syn_k, p_syn_k);
    mpfi_neg(syn_k, syn_k);

    // Set constants
    mpfi_set_d(one, 1.0);
    mpfi_set_d(two, 2.0);
    mpfi_set_d(neg_two, -2.0);

    // Initialise report files
    FILE **f_time = malloc(sizeof(FILE *));
    file_init("time", 1, f_time);

    FILE **f_nrn_N = malloc(p_nrn_size * sizeof(FILE *));
    file_init("nrn_N", p_nrn_size, f_nrn_N);
    FILE **f_nrn_V = malloc(p_nrn_size * sizeof(FILE *));
    file_init("nrn_V", p_nrn_size, f_nrn_V);

    /* FILE **f_syn_R = malloc(p_syn_size * sizeof(FILE *)); */
    /* file_init("syn_R", p_syn_size, f_syn_R); */
    /* FILE **f_syn_S = malloc(p_syn_size * sizeof(FILE *)); */
    /* file_init("syn_S", p_syn_size, f_syn_S); */


    // Begin simulation loop
    // =====================

    clock_t run_time = clock();

    for (i = 0; i < p_sim_steps; i++) {
        if (i % p_report_step == 0) printf("%lu\n", i);

        // Event(s) occur if urandom >= e^-rate
        for (j = 0; j < p_in_size; j++) {
            mpfr_urandom(rand_uf, rng_uf, MPFR_RNDN);
            in[j] = mpfr_greaterequal_p(rand_uf, in_p0);
            fprintf(stderr, "%s", (in[j] ? "\x1B[31m\xE2\x96\xA3\x1B[0m" : "\xE2\x96\xA3"));
        }
        fprintf(stderr, "\n");

        // Compute derivatives
        for (j = 0; j < p_nrn_size; j++) {
            dNdt(j);
            dVdt(j);
        }
        for (j = 0; j < p_syn_size; j++) {
            dRdt(j);
            dSdt(j);
        }

        // Step system
        mpfi_add(t, t, h);
        for (j = 0; j < p_nrn_size; j++) {
            mpfi_mul(&(d_nrn_N[j]), &(d_nrn_N[j]), h);
            mpfi_add(&(nrn_N[j]), &(nrn_N[j]), &(d_nrn_N[j]));
            mpfi_mul(&(d_nrn_V[j]), &(d_nrn_V[j]), h);
            mpfi_add(&(nrn_V[j]), &(nrn_V[j]), &(d_nrn_V[j]));
        }
        for (j = 0; j < p_syn_size; j++) {
            mpfi_mul(&(d_syn_R[j]), &(d_syn_R[j]), h);
            mpfi_add(&(syn_R[j]), &(syn_R[j]), &(d_syn_R[j]));
            mpfi_mul(&(d_syn_S[j]), &(d_syn_S[j]), h);
            mpfi_add(&(syn_S[j]), &(syn_S[j]), &(d_syn_S[j]));
        }

        file_write(t, 1, f_time);

        file_write(nrn_N, p_nrn_size, f_nrn_N);
        file_write(nrn_V, p_nrn_size, f_nrn_V);

        /* file_write(syn_R, p_syn_size, f_syn_R); */
        /* file_write(syn_S, p_syn_size, f_syn_S); */
    }

    run_time = clock() - run_time;
    printf("Finished in %f seconds.\n", ((float) run_time) / CLOCKS_PER_SEC);

    // End simulation loop
    // ===================


    // Clear system state
    mpfi_clear(h);
    mpfi_clear(t);
    for (i = 0; i < p_nrn_size; i++) {
        mpfi_clear(&(nrn_N[i]));
        mpfi_clear(&(nrn_V[i]));
        mpfi_clear(&(d_nrn_N[i]));
        mpfi_clear(&(d_nrn_V[i]));
    }
    for (i = 0; i < p_syn_size; i++) {
        mpfi_clear(&(syn_R[i]));
        mpfi_clear(&(syn_S[i]));
        mpfi_clear(&(d_syn_R[i]));
        mpfi_clear(&(d_syn_S[i]));
    }

    // Clear Poisson input parameters
    mpfr_clear(in_p0);
    mpfi_clear(in_V_lo);
    mpfi_clear(in_V_hi);
    for (i = 0; i < p_in_size; i++) {
        mpfr_clear(&(temp_sum_error[i]));
        mpfi_clear(&(I[i]));
    }

    // Clear neuron parameters
    mpfi_clear(nrn_GL);
    mpfi_clear(nrn_GCa);
    mpfi_clear(nrn_GK);
    mpfi_clear(nrn_VL);
    mpfi_clear(nrn_VCa);
    mpfi_clear(nrn_VK);
    mpfi_clear(nrn_V1);
    mpfi_clear(nrn_V2);
    mpfi_clear(nrn_V3);
    mpfi_clear(nrn_V4);
    mpfi_clear(nrn_phi);
    mpfi_clear(nrn_C);

    // Clear synapse parameters
    for (i = 0; i < p_syn_size; i++) {
        mpfi_clear(&(syn_GSyn[i]));
    }
    mpfi_clear(syn_VSyn);
    mpfi_clear(syn_thr);
    mpfi_clear(syn_a);
    mpfi_clear(syn_b);
    mpfi_clear(syn_k);

    // Clear constants
    mpfi_clear(one);
    mpfi_clear(two);
    mpfi_clear(neg_two);

    // Clear scratch space
    mpfr_clear(rand_uf);
    mpfr_clear(rand_nf);
    mpfr_clear(temp_error);
    mpfi_clear(temp_sum);
    mpfi_clear(temp1);
    mpfi_clear(temp2);
    mpfi_clear(M_ss);
    mpfi_clear(N_ss);

    // Free system state
    free(nrn_N);
    free(nrn_V);
    free(syn_R);
    free(syn_S);

    // Free other arrays
    free(d_nrn_N);
    free(d_nrn_V);
    free(d_syn_R);
    free(d_syn_S);
    free(syn_GSyn);
    free(temp_sum_error);
    free(temp_sum_error_ptr);
    free(I);
    free(in);

    // Clear report files
    file_clear(1, f_time);
    free(f_time);

    file_clear(p_nrn_size, f_nrn_N);
    free(f_nrn_N);
    file_clear(p_nrn_size, f_nrn_V);
    free(f_nrn_V);

    /* file_clear(p_syn_size, f_syn_R); */
    /* free(f_syn_R); */
    /* file_clear(p_syn_size, f_syn_S); */
    /* free(f_syn_S); */

    gmp_randclear(rng_uf);
    gmp_randclear(rng_nf);
    mpfr_free_cache();

    return 0;
}
