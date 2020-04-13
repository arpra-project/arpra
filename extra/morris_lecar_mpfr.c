/*
 * morris_lecar_mpfr.c -- Test Morris-Lecar model using MPFR.
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
#include <mpfr.h>

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
#define p_sim_steps 1000
#define p_report_step 20

// RNG parameters
// Seeds are random if not #defined
#define p_rand_prec 53
//#define p_rng_uf_seed 707135875931353ul
//#define p_rng_nf_seed 503108552855933ul
//#define p_rng_uz_seed 2071328946103ul

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
mpz_t rand_uz;
mpfr_t nrn_GL, nrn_VL, nrn_GCa, nrn_VCa, nrn_GK, nrn_VK, nrn_V1, nrn_V2, nrn_V3, nrn_V4,
    nrn_phi, nrn_C, syn_VSyn, syn_thr, syn_a, syn_b, syn_k, one, two, neg_two,
    temp1, temp2, M_ss, N_ss, in_V_lo, in_V_hi, in_p0, rand_uf, rand_nf;
mpfr_ptr syn_GSyn, I;
gmp_randstate_t rng_uf, rng_nf, rng_uz;

// System state variables
mpfr_ptr nrn_N, nrn_V, syn_R, syn_S;
mpfr_ptr d_nrn_N, d_nrn_V, d_syn_R, d_syn_S;

// DEBUG: print MPFR numbers to stderr
void debug (mpfr_srcptr x) {
    mpfr_out_str(stderr, 10, 80, x, MPFR_RNDN);
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

void file_write (mpfr_srcptr A, unsigned long grp_size, FILE **f)
{
    unsigned long i;

    for (i = 0; i < grp_size; i++) {
        mpfr_out_str(f[i], 10, 80, &(A[i]), MPFR_RNDN);
        fputc('\n', f[i]);
    }
}

void dNdt (const unsigned long idx)
{
    mpfr_ptr d_N = &(d_nrn_N[idx]);
    mpfr_srcptr N = &(nrn_N[idx]);
    mpfr_srcptr V = &(nrn_V[idx]);
    mpfr_srcptr V3 = nrn_V3;
    mpfr_srcptr V4 = nrn_V4;
    mpfr_srcptr phi = nrn_phi;

    // K+ channel activation steady-state
    // N_ss = 1 / (1 + exp(-2 (V - V3) / V4))
    mpfr_sub(temp1, V, V3, MPFR_RNDN);
    mpfr_mul(N_ss, neg_two, temp1, MPFR_RNDN);
    mpfr_div(N_ss, N_ss, V4, MPFR_RNDN);
    mpfr_exp(N_ss, N_ss, MPFR_RNDN);
    mpfr_add(N_ss, one, N_ss, MPFR_RNDN);
    mpfr_ui_div(N_ss, 1, N_ss, MPFR_RNDN);

    // tau of K+ channel activation
    // tau = 1 / (phi ((p + q) / 2))
    // p = exp(-(V - V3) / (2 V4))
    // q = exp( (V - V3) / (2 V4))
    mpfr_mul(temp2, two, V4, MPFR_RNDN);
    mpfr_div(temp2, temp1, temp2, MPFR_RNDN);
    mpfr_neg(temp1, temp2, MPFR_RNDN);
    mpfr_exp(temp1, temp1, MPFR_RNDN);
    mpfr_exp(temp2, temp2, MPFR_RNDN);
    mpfr_add(temp1, temp1, temp2, MPFR_RNDN);
    mpfr_div(temp1, temp1, two, MPFR_RNDN);
    mpfr_mul(temp1, phi, temp1, MPFR_RNDN);
    mpfr_ui_div(temp1, 1, temp1, MPFR_RNDN);

    // delta of K+ channel activation
    // dN/dt = (N_ss - N) / tau
    mpfr_sub(d_N, N_ss, N, MPFR_RNDN);
    mpfr_div(d_N, d_N, temp1, MPFR_RNDN);
}

void dVdt (const unsigned long idx)
{
    unsigned long i, j;
    unsigned long pre_size = p_in_size;
    mpfr_ptr d_V = &(d_nrn_V[idx]);
    mpfr_srcptr N = &(nrn_N[idx]);
    mpfr_srcptr V = &(nrn_V[idx]);
    mpfr_srcptr S = &(syn_S[idx * pre_size]);
    mpfr_srcptr GSyn = &(syn_GSyn[idx * pre_size]);
    mpfr_srcptr VSyn = syn_VSyn;
    mpfr_srcptr VL = nrn_VL;
    mpfr_srcptr GL = nrn_GL;
    mpfr_srcptr VCa = nrn_VCa;
    mpfr_srcptr GCa = nrn_GCa;
    mpfr_srcptr VK = nrn_VK;
    mpfr_srcptr GK = nrn_GK;
    mpfr_srcptr V1 = nrn_V1;
    mpfr_srcptr V2 = nrn_V2;
    mpfr_srcptr C = nrn_C;

    // Ca++ channel activation steady-state
    // M_ss = 1 / (1 + exp(-2 (V - V1) / V2))
    mpfr_sub(M_ss, V, V1, MPFR_RNDN);
    mpfr_mul(M_ss, neg_two, M_ss, MPFR_RNDN);
    mpfr_div(M_ss, M_ss, V2, MPFR_RNDN);
    mpfr_exp(M_ss, M_ss, MPFR_RNDN);
    mpfr_add(M_ss, one, M_ss, MPFR_RNDN);
    mpfr_ui_div(M_ss, 1, M_ss, MPFR_RNDN);

    // Synapse current
    mpfr_sub(temp1, VSyn, V, MPFR_RNDN);
    mpfr_set_ui(d_V, 0, MPFR_RNDN);
    for (i = 0; i < pre_size; i++) {
        mpfr_mul(&(I[i]), temp1, &(GSyn[i]), MPFR_RNDN);
        mpfr_mul(&(I[i]), &(I[i]), &(S[i]), MPFR_RNDN);
    }

    // Shuffle input currents with Fisher-Yates, then sum.
    for (i = 0; i < pre_size; i++) {
        mpz_set_ui(rand_uz, pre_size - i);
        mpz_urandomm(rand_uz, rng_uz, rand_uz);
        j = i + mpz_get_ui(rand_uz);
        mpfr_swap(&(I[i]), &(I[j]));
        mpfr_add(d_V, d_V, &(I[i]), MPFR_RNDN);
    }




    // ======== TEMP DEBUG ==========
    //mpfr_set_d(d_V, 80.0, MPFR_RNDN); // bifurcation at sum(I) = 80.0
    if (idx == 0) {
        //fprintf(stderr, "sum(I): "); debug(d_V);
        //fprintf(stderr, "I[0]: "); debug(I);
        for (i = 0; i < pre_size; i++) {
            //fprintf(stderr, "I[%lu]: ", i); debug(&(I[i]));
        }
    }



    // Leak current
    mpfr_sub(temp1, V, VL, MPFR_RNDN);
    mpfr_mul(temp1, temp1, GL, MPFR_RNDN);
    mpfr_sub(d_V, d_V, temp1, MPFR_RNDN);

    // Ca++ current
    mpfr_sub(temp1, V, VCa, MPFR_RNDN);
    mpfr_mul(temp1, temp1, GCa, MPFR_RNDN);
    mpfr_mul(temp1, temp1, M_ss, MPFR_RNDN);
    mpfr_sub(d_V, d_V, temp1, MPFR_RNDN);

    // K+ current
    mpfr_sub(temp1, V, VK, MPFR_RNDN);
    mpfr_mul(temp1, temp1, GK, MPFR_RNDN);
    mpfr_mul(temp1, temp1, N, MPFR_RNDN);
    mpfr_sub(d_V, d_V, temp1, MPFR_RNDN);

    // delta of membrane potential
    // dV/dt = (I + GL (VL - V) + GCa M (VCa - V) + GK N (VK - V)) / C
    mpfr_div(d_V, d_V, C, MPFR_RNDN);
}

void dRdt (const unsigned long idx)
{
    mpfr_ptr d_R = &(d_syn_R[idx]);
    mpfr_srcptr R = &(syn_R[idx]);
    mpfr_srcptr a = syn_a;
    mpfr_srcptr b = syn_b;
    mpfr_srcptr k = syn_k;
    mpfr_srcptr VPre = in[idx % p_in_size] ? in_V_hi : in_V_lo;
    mpfr_srcptr threshold = syn_thr;

    // Sigmoid of threshold difference
    mpfr_sub(temp1, VPre, threshold, MPFR_RNDN);
    mpfr_mul(temp1, temp1, k, MPFR_RNDN);
    mpfr_exp(temp1, temp1, MPFR_RNDN);
    mpfr_add(temp1, temp1, one, MPFR_RNDN);
    mpfr_ui_div(temp1, 1, temp1, MPFR_RNDN);

    // Presynaptic transmitter release rise
    mpfr_mul(temp1, a, temp1, MPFR_RNDN);

    // Presynaptic transmitter release decay
    mpfr_mul(temp2, b, R, MPFR_RNDN);

    // delta of presynaptic transmitter release
    // dR/dt = a Q - b R
    // Q = 1 / (1 + e^(k(V - threshold)))
    mpfr_sub(d_R, temp1, temp2, MPFR_RNDN);
}

void dSdt (const unsigned long idx)
{
    mpfr_ptr d_S = &(d_syn_S[idx]);
    mpfr_srcptr R = &(syn_R[idx]);
    mpfr_srcptr S = &(syn_S[idx]);
    mpfr_srcptr a = syn_a;
    mpfr_srcptr b = syn_b;

    // Postsynaptic transmitter binding rise
    mpfr_mul(temp1, a, R, MPFR_RNDN);

    // Postsynaptic transmitter binding decay
    mpfr_mul(temp2, b, S, MPFR_RNDN);

    // delta of postsynaptic transmitter binding
    // dS/dt = a R - b S
    mpfr_sub(d_S, temp1, temp2, MPFR_RNDN);
}

int main (int argc, char *argv[])
{
    mpfr_t h, t;
    unsigned long i, j;

    // Allocate system state
    nrn_N = malloc(p_nrn_size * sizeof(mpfr_t));
    nrn_V = malloc(p_nrn_size * sizeof(mpfr_t));
    syn_R = malloc(p_syn_size * sizeof(mpfr_t));
    syn_S = malloc(p_syn_size * sizeof(mpfr_t));

    // Allocate other arrays
    d_nrn_N = malloc(p_nrn_size * sizeof(mpfr_t));
    d_nrn_V = malloc(p_nrn_size * sizeof(mpfr_t));
    d_syn_R = malloc(p_syn_size * sizeof(mpfr_t));
    d_syn_S = malloc(p_syn_size * sizeof(mpfr_t));
    syn_GSyn = malloc(p_syn_size * sizeof(mpfr_t));
    I = malloc(p_in_size * sizeof(mpfr_t));
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

    // Initialise uniform integer RNG
    gmp_randinit_default(rng_uz);
    clock_gettime(CLOCK_REALTIME, &clock_time);
#ifdef p_rng_uz_seed
    unsigned long rng_uz_seed = p_rng_uz_seed;
#else
    unsigned long rng_uz_seed = clock_time.tv_sec + clock_time.tv_nsec;
#endif
    gmp_randseed_ui(rng_uz, rng_uz_seed);
    printf("GMP rand uniform integer seed: %lu\n", rng_uz_seed);

    // Initialise system state
    mpfr_init2(h, p_prec);
    mpfr_init2(t, p_prec);
    for (i = 0; i < p_nrn_size; i++) {
        mpfr_init2(&(nrn_N[i]), p_prec);
        mpfr_init2(&(nrn_V[i]), p_prec);
        mpfr_init2(&(d_nrn_N[i]), p_prec);
        mpfr_init2(&(d_nrn_V[i]), p_prec);
    }
    for (i = 0; i < p_syn_size; i++) {
        mpfr_init2(&(syn_R[i]), p_prec);
        mpfr_init2(&(syn_S[i]), p_prec);
        mpfr_init2(&(d_syn_R[i]), p_prec);
        mpfr_init2(&(d_syn_S[i]), p_prec);
    }

    // Initialise Poisson input parameters
    mpfr_init2(in_p0, p_prec);
    mpfr_init2(in_V_lo, p_prec);
    mpfr_init2(in_V_hi, p_prec);

    // Initialise neuron parameters
    mpfr_init2(nrn_GL, p_prec);
    mpfr_init2(nrn_GCa, p_prec);
    mpfr_init2(nrn_GK, p_prec);
    mpfr_init2(nrn_VL, p_prec);
    mpfr_init2(nrn_VCa, p_prec);
    mpfr_init2(nrn_VK, p_prec);
    mpfr_init2(nrn_V1, p_prec);
    mpfr_init2(nrn_V2, p_prec);
    mpfr_init2(nrn_V3, p_prec);
    mpfr_init2(nrn_V4, p_prec);
    mpfr_init2(nrn_phi, p_prec);
    mpfr_init2(nrn_C, p_prec);

    // Initialise synapse parameters
    for (i = 0; i < p_syn_size; i++) {
        mpfr_init2(&(syn_GSyn[i]), p_prec);
    }
    mpfr_init2(syn_VSyn, p_prec);
    mpfr_init2(syn_thr, p_prec);
    mpfr_init2(syn_a, p_prec);
    mpfr_init2(syn_b, p_prec);
    mpfr_init2(syn_k, p_prec);

    // Initialise constants
    mpfr_init2(one, p_prec);
    mpfr_init2(two, p_prec);
    mpfr_init2(neg_two, p_prec);

    // Initialise scratch space
    mpfr_init2(rand_uf, p_rand_prec);
    mpfr_init2(rand_nf, p_rand_prec);
    mpz_init(rand_uz);
    mpfr_init2(temp1, p_prec);
    mpfr_init2(temp2, p_prec);
    mpfr_init2(M_ss, p_prec);
    mpfr_init2(N_ss, p_prec);
    for (i = 0; i < p_in_size; i++) {
        mpfr_init2(&(I[i]), p_prec);
    }

    // Set system state
    mpfr_set_d(h, p_h, MPFR_RNDN);
    mpfr_set_d(t, p_t0, MPFR_RNDN);
    for (i = 0; i < p_nrn_size; i++) {
        mpfr_set_d(&(nrn_N[i]), p_nrn_N0, MPFR_RNDN);
        mpfr_set_d(&(nrn_V[i]), p_nrn_V0, MPFR_RNDN);
    }
    for (i = 0; i < p_syn_size; i++) {
        mpfr_set_d(&(syn_R[i]), p_syn_R0, MPFR_RNDN);
        mpfr_set_d(&(syn_S[i]), p_syn_S0, MPFR_RNDN);
    }

    // Set Poisson input parameters
    mpfr_set_d(in_p0, -(p_in_freq / 1000) * p_h, MPFR_RNDN);
    mpfr_exp(in_p0, in_p0, MPFR_RNDN);
    mpfr_set_d(in_V_lo, p_in_V_lo, MPFR_RNDN);
    mpfr_set_d(in_V_hi, p_in_V_hi, MPFR_RNDN);

    // Set neuron parameters
    mpfr_set_d(nrn_GL, p_nrn_GL, MPFR_RNDN);
    mpfr_set_d(nrn_GCa, p_nrn_GCa, MPFR_RNDN);
    mpfr_set_d(nrn_GK, p_nrn_GK, MPFR_RNDN);
    mpfr_set_d(nrn_VL, p_nrn_VL, MPFR_RNDN);
    mpfr_set_d(nrn_VCa, p_nrn_VCa, MPFR_RNDN);
    mpfr_set_d(nrn_VK, p_nrn_VK, MPFR_RNDN);
    mpfr_set_d(nrn_V1, p_nrn_V1, MPFR_RNDN);
    mpfr_set_d(nrn_V2, p_nrn_V2, MPFR_RNDN);
    mpfr_set_d(nrn_V3, p_nrn_V3, MPFR_RNDN);
    mpfr_set_d(nrn_V4, p_nrn_V4, MPFR_RNDN);
    mpfr_set_d(nrn_phi, p_nrn_phi, MPFR_RNDN);
    mpfr_set_d(nrn_C, p_nrn_C, MPFR_RNDN);

    // Set synapse parameters
    for (i = 0; i < p_syn_size; i++) {
        //mpfr_nrandom(rand_nf, rng_nf, MPFR_RNDN);
        mpfr_grandom(rand_nf, NULL, rng_nf, MPFR_RNDN);
        mpfr_mul_d(rand_nf, rand_nf, p_syn_GSyn_std, MPFR_RNDN);
        mpfr_add_d(&(syn_GSyn[i]), rand_nf, p_syn_GSyn_mean, MPFR_RNDN);
    }
    mpfr_set_d(syn_VSyn, p_syn_VSyn, MPFR_RNDN);
    mpfr_set_d(syn_thr, p_syn_thr, MPFR_RNDN);
    mpfr_set_d(syn_a, p_syn_a, MPFR_RNDN);
    mpfr_set_d(syn_b, p_syn_b, MPFR_RNDN);
    mpfr_set_d(syn_k, p_syn_k, MPFR_RNDN);
    mpfr_neg(syn_k, syn_k, MPFR_RNDN);

    // Set constants
    mpfr_set_d(one, 1.0, MPFR_RNDN);
    mpfr_set_d(two, 2.0, MPFR_RNDN);
    mpfr_set_d(neg_two, -2.0, MPFR_RNDN);

    // Initialise report files
    FILE **f_time = malloc(sizeof(FILE *));;
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
        mpfr_add(t, t, h, MPFR_RNDN);
        for (j = 0; j < p_nrn_size; j++) {
            mpfr_mul(&(d_nrn_N[j]), &(d_nrn_N[j]), h, MPFR_RNDN);
            mpfr_add(&(nrn_N[j]), &(nrn_N[j]), &(d_nrn_N[j]), MPFR_RNDN);
            mpfr_mul(&(d_nrn_V[j]), &(d_nrn_V[j]), h, MPFR_RNDN);
            mpfr_add(&(nrn_V[j]), &(nrn_V[j]), &(d_nrn_V[j]), MPFR_RNDN);
        }
        for (j = 0; j < p_syn_size; j++) {
            mpfr_mul(&(d_syn_R[j]), &(d_syn_R[j]), h, MPFR_RNDN);
            mpfr_add(&(syn_R[j]), &(syn_R[j]), &(d_syn_R[j]), MPFR_RNDN);
            mpfr_mul(&(d_syn_S[j]), &(d_syn_S[j]), h, MPFR_RNDN);
            mpfr_add(&(syn_S[j]), &(syn_S[j]), &(d_syn_S[j]), MPFR_RNDN);
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
    mpfr_clear(h);
    mpfr_clear(t);
    for (i = 0; i < p_nrn_size; i++) {
        mpfr_clear(&(nrn_N[i]));
        mpfr_clear(&(nrn_V[i]));
        mpfr_clear(&(d_nrn_N[i]));
        mpfr_clear(&(d_nrn_V[i]));
    }
    for (i = 0; i < p_syn_size; i++) {
        mpfr_clear(&(syn_R[i]));
        mpfr_clear(&(syn_S[i]));
        mpfr_clear(&(d_syn_R[i]));
        mpfr_clear(&(d_syn_S[i]));
    }

    // Clear Poisson input parameters
    mpfr_clear(in_p0);
    mpfr_clear(in_V_lo);
    mpfr_clear(in_V_hi);

    // Clear neuron parameters
    mpfr_clear(nrn_GL);
    mpfr_clear(nrn_GCa);
    mpfr_clear(nrn_GK);
    mpfr_clear(nrn_VL);
    mpfr_clear(nrn_VCa);
    mpfr_clear(nrn_VK);
    mpfr_clear(nrn_V1);
    mpfr_clear(nrn_V2);
    mpfr_clear(nrn_V3);
    mpfr_clear(nrn_V4);
    mpfr_clear(nrn_phi);
    mpfr_clear(nrn_C);

    // Clear synapse parameters
    for (i = 0; i < p_syn_size; i++) {
        mpfr_clear(&(syn_GSyn[i]));
    }
    mpfr_clear(syn_VSyn);
    mpfr_clear(syn_thr);
    mpfr_clear(syn_a);
    mpfr_clear(syn_b);
    mpfr_clear(syn_k);

    // Clear constants
    mpfr_clear(one);
    mpfr_clear(two);
    mpfr_clear(neg_two);

    // Clear scratch space
    mpfr_clear(rand_uf);
    mpfr_clear(rand_nf);
    mpz_clear(rand_uz);
    mpfr_clear(temp1);
    mpfr_clear(temp2);
    mpfr_clear(M_ss);
    mpfr_clear(N_ss);
    for (i = 0; i < p_in_size; i++) {
        mpfr_clear(&(I[i]));
    }

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
    gmp_randclear(rng_uz);
    mpfr_free_cache();

    return 0;
}
