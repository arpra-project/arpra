/*
 * morris_lecar_mpfr.c -- Test Morris-Lecar model using MPFR.
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

// Poisson input parameters (group 1)
#define p_in1_size 50
#define p_in1_freq 5.0
#define p_in1_V_lo -60.0
#define p_in1_V_hi 20.0

// Poisson input parameters (group 2)
#define p_in2_size 0
#define p_in2_freq 5.0
#define p_in2_V_lo -60.0
#define p_in2_V_hi 20.0

// Neuron parameters (group 1)
#define p_nrn1_size 1
#define p_nrn1_N0 0.0
#define p_nrn1_V0 -60.0
#define p_nrn1_class 1

// Neuron parameters (group 2)
#define p_nrn2_size 0
#define p_nrn2_N0 0.0
#define p_nrn2_V0 -60.0
#define p_nrn2_class 1

// Neuron parameters (common)
#define p_GL 2.0
#define p_GCa 4.0 // Class 1
//#define p_GCa 4.4 // Class 2
#define p_GK 8.0
#define p_VL -60.0
#define p_VCa 120.0
#define p_VK -80.0
#define p_V1 -1.2
#define p_V2 18.0
#define p_V3 12.0 // Class 1
//#define p_V3 2.0 // Class 2
#define p_V4 17.4 // Class 1
//#define p_V4 30.0 // Class 2
#define p_phi 1.0 / 15.0 // Class 1
//#define p_phi 1.0 / 25.0 // Class 2
#define p_C 20.0

// Synapse parameters (excitatory)
#define p_syn_exc_size p_in1_size * p_nrn1_size
#define p_syn_exc_R0 0.0
#define p_syn_exc_S0 0.0
#define p_syn_exc_GSyn_std 0.5
#define p_syn_exc_GSyn_mean 3.0
#define p_syn_exc_VSyn 0.0
#define p_syn_exc_thr -50.0
#define p_syn_exc_a 0.25 // in [1/10, 1/2]
#define p_syn_exc_b 0.15 // in [1/20, 1/4]
#define p_syn_exc_k 1.0E6

// Synapse parameters (inhibitory)
#define p_syn_inh_size 0
#define p_syn_inh_R0 0.0
#define p_syn_inh_S0 0.0
#define p_syn_inh_GSyn_std 0.5
#define p_syn_inh_GSyn_mean 3.0
#define p_syn_inh_VSyn -80.0
#define p_syn_inh_thr -50.0
#define p_syn_inh_a 0.075 // in [1/20, 1/10]
#define p_syn_inh_b 0.035 // in [1/50, 1/20]
#define p_syn_inh_k 1.0E6

// ===================== end of model parameters ======================


int *in1, *in2;
mpfr_t GL, VL, GCa, VCa, GK, VK, V1, V2, V3, V4, phi, C, syn_exc_VSyn, syn_exc_thr,
       syn_exc_a, syn_exc_b, syn_exc_k, syn_inh_VSyn, syn_inh_thr, syn_inh_a,
       syn_inh_b, syn_inh_k, one, two, neg_two, temp1, temp2, M_ss, N_ss, in1_V_lo,
       in1_V_hi, in2_V_lo, in2_V_hi, in1_p0, in2_p0, rand_uf, rand_nf;
mpz_t rand_uz;
mpfr_ptr syn_exc_GSyn, syn_inh_GSyn, I1, I2;
gmp_randstate_t rng_uf, rng_nf, rng_uz;
unsigned long rng_uf_seed, rng_nf_seed, rng_uz_seed;

// System state variables
mpfr_ptr nrn1_N, nrn2_N, nrn1_V, nrn2_V, syn_exc_R, syn_inh_R, syn_exc_S, syn_inh_S;
mpfr_ptr d_nrn1_N, d_nrn2_N, d_nrn1_V, d_nrn2_V, d_syn_exc_R, d_syn_inh_R, d_syn_exc_S, d_syn_inh_S;

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
        sprintf(fname, "%s_%03u.dat", grp, (unsigned) i);
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
    unsigned long i, j;

    for (i = 0; i < grp_size; i++) {
        mpfr_out_str(f[i], 10, 80, &(A[i]), MPFR_RNDN);
        fputc('\n', f[i]);
    }
}

void dNdt (const unsigned long idx, int grp)
{
    mpfr_ptr d_N;
    mpfr_srcptr N, V;

    if (grp == 1) {
        d_N = d_nrn1_N + idx;
        N = nrn1_N + idx;
        V = nrn1_V + idx;
    }
    //else if (grp == 2) {
    //   d_N = d_nrn2_N + idx;
    //   N = nrn2_N + idx;
    //   V = nrn2_V + idx;
    //}

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

void dVdt (const unsigned long idx, int grp)
{
    mpfr_ptr d_V;
    mpfr_srcptr N, V;
    mpfr_srcptr S, GSyn, VSyn;
    mpfr_ptr I;
    unsigned long i, j, pre_size;

    if (grp == 1) {
        d_V = d_nrn1_V + idx;
        N = nrn1_N + idx;
        V = nrn1_V + idx;
        pre_size = p_in1_size;
        S = syn_exc_S + (idx * pre_size);
        GSyn = syn_exc_GSyn + (idx * pre_size);
        VSyn = syn_exc_VSyn;
        I = I1;
    }
    //else if (grp == 2) {
    //   d_V = d_nrn2_V + idx;
    //   N = nrn2_N + idx;
    //   V = nrn2_V + idx;
    //   pre_size = p_in2_size;
    //   S = syn_exc_S + (idx * pre_size);
    //   GSyn = syn_exc_GSyn + (idx * pre_size);
    //   VSyn = syn_exc_VSyn;
    //   I = I2;
    //}

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

void dRdt (const unsigned long idx, int grp)
{
    mpfr_ptr d_R;
    mpfr_srcptr R;
    mpfr_srcptr a, b, k;
    mpfr_srcptr VPre, threshold;

    if (grp == 1) {
        d_R = d_syn_exc_R + idx;
        R = syn_exc_R + idx;
        a = syn_exc_a;
        b = syn_exc_b;
        k = syn_exc_k;
        VPre = in1[idx % p_in1_size] ? in1_V_hi : in1_V_lo;
        threshold = syn_exc_thr;
    }
    //else if (grp == 2) {
    //   d_R = d_syn_inh_R + idx;
    //   R = syn_inh_R + idx;
    //   a = syn_inh_a;
    //   b = syn_inh_b;
    //   k = syn_inh_k;
    //   VPre = in2[idx % p_in2_size] ? in2_V_hi : in2_V_lo;
    //   threshold = syn_inh_thr;
    //}

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

void dSdt (const unsigned long idx, int grp)
{
    mpfr_ptr d_S;
    mpfr_srcptr R, S;
    mpfr_srcptr a, b;

    if (grp == 1) {
        d_S = d_syn_exc_S + idx;
        R = syn_exc_R + idx;
        S = syn_exc_S + idx;
        a = syn_exc_a;
        b = syn_exc_b;
    }
    //else if (grp == 2) {
    //   d_S = d_syn_inh_S + idx;
    //   R = syn_inh_R + idx;
    //   S = syn_inh_S + idx;
    //   a = syn_inh_a;
    //   b = syn_inh_b;
    //}

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
    struct timespec sys_t;
    clock_t run_time;
    unsigned long i, j;

    // Allocate dynamic arrays
    nrn1_N = malloc(p_nrn1_size * sizeof(mpfr_t));
    nrn2_N = malloc(p_nrn2_size * sizeof(mpfr_t));
    nrn1_V = malloc(p_nrn1_size * sizeof(mpfr_t));
    nrn2_V = malloc(p_nrn2_size * sizeof(mpfr_t));
    syn_exc_R = malloc(p_syn_exc_size * sizeof(mpfr_t));
    syn_inh_R = malloc(p_syn_inh_size * sizeof(mpfr_t));
    syn_exc_S = malloc(p_syn_exc_size * sizeof(mpfr_t));
    syn_inh_S = malloc(p_syn_inh_size * sizeof(mpfr_t));
    d_nrn1_N = malloc(p_nrn1_size * sizeof(mpfr_t));
    d_nrn2_N = malloc(p_nrn2_size * sizeof(mpfr_t));
    d_nrn1_V = malloc(p_nrn1_size * sizeof(mpfr_t));
    d_nrn2_V = malloc(p_nrn2_size * sizeof(mpfr_t));
    d_syn_exc_R = malloc(p_syn_exc_size * sizeof(mpfr_t));
    d_syn_inh_R = malloc(p_syn_inh_size * sizeof(mpfr_t));
    d_syn_exc_S = malloc(p_syn_exc_size * sizeof(mpfr_t));
    d_syn_inh_S = malloc(p_syn_inh_size * sizeof(mpfr_t));
    syn_exc_GSyn = malloc(p_syn_exc_size * sizeof(mpfr_t));
    syn_inh_GSyn = malloc(p_syn_inh_size * sizeof(mpfr_t));
    I1 = malloc(p_in1_size * sizeof(mpfr_t));
    I2 = malloc(p_in2_size * sizeof(mpfr_t));
    in1 = malloc(p_in1_size * sizeof(int));
    in2 = malloc(p_in2_size * sizeof(int));

    // Initialise uniform float RNG
    gmp_randinit_default(rng_uf);
    clock_gettime(CLOCK_REALTIME, &sys_t);
    //rng_uf_seed = 707135875931353ul;
    rng_uf_seed = sys_t.tv_sec + sys_t.tv_nsec;
    gmp_randseed_ui(rng_uf, rng_uf_seed);
    printf("GMP rand uniform float seed: %lu\n", rng_uf_seed);

    // Initialise normal float RNG
    gmp_randinit_default(rng_nf);
    clock_gettime(CLOCK_REALTIME, &sys_t);
    //rng_nf_seed = 503108552855933ul;
    rng_nf_seed = sys_t.tv_sec + sys_t.tv_nsec;
    gmp_randseed_ui(rng_nf, rng_nf_seed);
    printf("GMP rand normal float seed: %lu\n", rng_nf_seed);

    // Initialise uniform integer RNG
    gmp_randinit_default(rng_uz);
    clock_gettime(CLOCK_REALTIME, &sys_t);
    //rng_uz_seed = 2071328946103ul;
    rng_uz_seed = sys_t.tv_sec + sys_t.tv_nsec;
    gmp_randseed_ui(rng_uz, rng_uz_seed);
    printf("GMP rand uniform integer seed: %lu\n", rng_uz_seed);

    // Initialise system state
    mpfr_init2(h, p_prec);
    mpfr_init2(t, p_prec);
    for (i = 0; i < p_nrn1_size; i++) {
        mpfr_init2(&(nrn1_N[i]), p_prec);
        mpfr_init2(&(nrn1_V[i]), p_prec);
        mpfr_init2(&(d_nrn1_N[i]), p_prec);
        mpfr_init2(&(d_nrn1_V[i]), p_prec);
    }
    for (i = 0; i < p_nrn2_size; i++) {
        mpfr_init2(&(nrn2_N[i]), p_prec);
        mpfr_init2(&(nrn2_V[i]), p_prec);
        mpfr_init2(&(d_nrn2_N[i]), p_prec);
        mpfr_init2(&(d_nrn2_V[i]), p_prec);
    }
    for (i = 0; i < p_syn_exc_size; i++) {
        mpfr_init2(&(syn_exc_R[i]), p_prec);
        mpfr_init2(&(syn_exc_S[i]), p_prec);
        mpfr_init2(&(d_syn_exc_R[i]), p_prec);
        mpfr_init2(&(d_syn_exc_S[i]), p_prec);
    }
    for (i = 0; i < p_syn_inh_size; i++) {
        mpfr_init2(&(syn_inh_R[i]), p_prec);
        mpfr_init2(&(syn_inh_S[i]), p_prec);
        mpfr_init2(&(d_syn_inh_R[i]), p_prec);
        mpfr_init2(&(d_syn_inh_S[i]), p_prec);
    }

    // Initialise Poisson input parameters (group 1)
    mpfr_init2(in1_p0, p_prec);
    mpfr_init2(in1_V_lo, p_prec);
    mpfr_init2(in1_V_hi, p_prec);

    // Initialise Poisson input parameters (group 2)
    mpfr_init2(in2_p0, p_prec);
    mpfr_init2(in2_V_lo, p_prec);
    mpfr_init2(in2_V_hi, p_prec);

    // Initialise neuron parameters
    mpfr_init2(GL, p_prec);
    mpfr_init2(GCa, p_prec);
    mpfr_init2(GK, p_prec);
    mpfr_init2(VL, p_prec);
    mpfr_init2(VCa, p_prec);
    mpfr_init2(VK, p_prec);
    mpfr_init2(V1, p_prec);
    mpfr_init2(V2, p_prec);
    mpfr_init2(V3, p_prec);
    mpfr_init2(V4, p_prec);
    mpfr_init2(phi, p_prec);
    mpfr_init2(C, p_prec);

    // Initialise excitatory synapse parameters
    for (i = 0; i < p_syn_exc_size; i++) {
        mpfr_init2(&(syn_exc_GSyn[i]), p_prec);
    }
    mpfr_init2(syn_exc_VSyn, p_prec);
    mpfr_init2(syn_exc_thr, p_prec);
    mpfr_init2(syn_exc_a, p_prec);
    mpfr_init2(syn_exc_b, p_prec);
    mpfr_init2(syn_exc_k, p_prec);

    // Initialise inhibitory synapse parameters
    for (i = 0; i < p_syn_inh_size; i++) {
        mpfr_init2(&(syn_inh_GSyn[i]), p_prec);
    }
    mpfr_init2(syn_inh_VSyn, p_prec);
    mpfr_init2(syn_inh_thr, p_prec);
    mpfr_init2(syn_inh_a, p_prec);
    mpfr_init2(syn_inh_b, p_prec);
    mpfr_init2(syn_inh_k, p_prec);

    // Initialise constants
    mpfr_init2(one, p_prec);
    mpfr_init2(two, p_prec);
    mpfr_init2(neg_two, p_prec);

    // Initialise scratch space
    mpfr_init2(rand_uf, p_prec);
    mpfr_init2(rand_nf, p_prec);
    mpz_init(rand_uz);
    mpfr_init2(temp1, p_prec);
    mpfr_init2(temp2, p_prec);
    mpfr_init2(M_ss, p_prec);
    mpfr_init2(N_ss, p_prec);
    for (i = 0; i < p_in1_size; i++) {
        mpfr_init2(&(I1[i]), p_prec);
    }
    for (i = 0; i < p_in2_size; i++) {
        mpfr_init2(&(I2[i]), p_prec);
    }

    // Set system state
    mpfr_set_d(h, p_h, MPFR_RNDN);
    mpfr_set_d(t, p_t0, MPFR_RNDN);
    for (i = 0; i < p_nrn1_size; i++) {
        mpfr_set_d(&(nrn1_N[i]), p_nrn1_N0, MPFR_RNDN);
        mpfr_set_d(&(nrn1_V[i]), p_nrn1_V0, MPFR_RNDN);
    }
    for (i = 0; i < p_nrn2_size; i++) {
        mpfr_set_d(&(nrn2_N[i]), p_nrn2_N0, MPFR_RNDN);
        mpfr_set_d(&(nrn2_V[i]), p_nrn2_V0, MPFR_RNDN);
    }
    for (i = 0; i < p_syn_exc_size; i++) {
        mpfr_set_d(&(syn_exc_R[i]), p_syn_exc_R0, MPFR_RNDN);
        mpfr_set_d(&(syn_exc_S[i]), p_syn_exc_S0, MPFR_RNDN);
    }
    for (i = 0; i < p_syn_inh_size; i++) {
        mpfr_set_d(&(syn_inh_R[i]), p_syn_inh_R0, MPFR_RNDN);
        mpfr_set_d(&(syn_inh_S[i]), p_syn_inh_S0, MPFR_RNDN);
    }

    // Set Poisson input parameters (group 1)
    mpfr_set_d(in1_p0, -(p_in1_freq / 1000) * p_h, MPFR_RNDN);
    mpfr_exp(in1_p0, in1_p0, MPFR_RNDN);
    mpfr_set_d(in1_V_lo, p_in1_V_lo, MPFR_RNDN);
    mpfr_set_d(in1_V_hi, p_in1_V_hi, MPFR_RNDN);

    // Set Poisson input parameters (group 2)
    mpfr_set_d(in2_p0, -(p_in2_freq / 1000) * p_h, MPFR_RNDN);
    mpfr_exp(in2_p0, in2_p0, MPFR_RNDN);
    mpfr_set_d(in2_V_lo, p_in2_V_lo, MPFR_RNDN);
    mpfr_set_d(in2_V_hi, p_in2_V_hi, MPFR_RNDN);

    // Set neuron parameters
    mpfr_set_d(GL, p_GL, MPFR_RNDN);
    mpfr_set_d(GCa, p_GCa, MPFR_RNDN);
    mpfr_set_d(GK, p_GK, MPFR_RNDN);
    mpfr_set_d(VL, p_VL, MPFR_RNDN);
    mpfr_set_d(VCa, p_VCa, MPFR_RNDN);
    mpfr_set_d(VK, p_VK, MPFR_RNDN);
    mpfr_set_d(V1, p_V1, MPFR_RNDN);
    mpfr_set_d(V2, p_V2, MPFR_RNDN);
    mpfr_set_d(V3, p_V3, MPFR_RNDN);
    mpfr_set_d(V4, p_V4, MPFR_RNDN);
    mpfr_set_d(phi, p_phi, MPFR_RNDN);
    mpfr_set_d(C, p_C, MPFR_RNDN);

    // Set excitatory synapse parameters
    for (i = 0; i < p_syn_exc_size; i++) {
        //mpfr_nrandom(rand_nf, rng_nf, MPFR_RNDN);
        mpfr_grandom(rand_nf, NULL, rng_nf, MPFR_RNDN);
        mpfr_mul_d(rand_nf, rand_nf, p_syn_exc_GSyn_std, MPFR_RNDN);
        mpfr_add_d(&(syn_exc_GSyn[i]), rand_nf, p_syn_exc_GSyn_mean, MPFR_RNDN);
    }
    mpfr_set_d(syn_exc_VSyn, p_syn_exc_VSyn, MPFR_RNDN);
    mpfr_set_d(syn_exc_thr, p_syn_exc_thr, MPFR_RNDN);
    mpfr_set_d(syn_exc_a, p_syn_exc_a, MPFR_RNDN);
    mpfr_set_d(syn_exc_b, p_syn_exc_b, MPFR_RNDN);
    mpfr_set_d(syn_exc_k, p_syn_exc_k, MPFR_RNDN);
    mpfr_neg(syn_exc_k, syn_exc_k, MPFR_RNDN);

    // Set inhibitory synapse parameters
    for (i = 0; i < p_syn_inh_size; i++) {
        //mpfr_nrandom(rand_nf, rng_nf, MPFR_RNDN);
        mpfr_grandom(rand_nf, NULL, rng_nf, MPFR_RNDN);
        mpfr_mul_d(rand_nf, rand_nf, p_syn_inh_GSyn_std, MPFR_RNDN);
        mpfr_add_d(&(syn_inh_GSyn[i]), rand_nf, p_syn_inh_GSyn_mean, MPFR_RNDN);
    }
    mpfr_set_d(syn_inh_VSyn, p_syn_inh_VSyn, MPFR_RNDN);
    mpfr_set_d(syn_inh_thr, p_syn_inh_thr, MPFR_RNDN);
    mpfr_set_d(syn_inh_a, p_syn_inh_a, MPFR_RNDN);
    mpfr_set_d(syn_inh_b, p_syn_inh_b, MPFR_RNDN);
    mpfr_set_d(syn_inh_k, p_syn_inh_k, MPFR_RNDN);
    mpfr_neg(syn_inh_k, syn_inh_k, MPFR_RNDN);

    // Set constants
    mpfr_set_d(one, 1.0, MPFR_RNDN);
    mpfr_set_d(two, 2.0, MPFR_RNDN);
    mpfr_set_d(neg_two, -2.0, MPFR_RNDN);

    // Initialise report files
    FILE **f_time = malloc(sizeof(FILE *));;
    file_init("time", 1, f_time);

    FILE **f_nrn1_N = malloc(p_nrn1_size * sizeof(FILE *));
    file_init("nrn1_N", p_nrn1_size, f_nrn1_N);
    FILE **f_nrn1_V = malloc(p_nrn1_size * sizeof(FILE *));
    file_init("nrn1_V", p_nrn1_size, f_nrn1_V);

    //FILE **f_nrn2_N = malloc(p_nrn2_size * sizeof(FILE *));
    //file_init("nrn2_N", p_nrn2_size, f_nrn2_N);
    //FILE **f_nrn2_V = malloc(p_nrn2_size * sizeof(FILE *));
    //file_init("nrn2_V", p_nrn2_size, f_nrn2_V);

    FILE **f_syn_exc_R = malloc(p_syn_exc_size * sizeof(FILE *));
    file_init("syn_exc_R", p_syn_exc_size, f_syn_exc_R);
    FILE **f_syn_exc_S = malloc(p_syn_exc_size * sizeof(FILE *));
    file_init("syn_exc_S", p_syn_exc_size, f_syn_exc_S);

    //FILE **f_syn_inh_R = malloc(p_syn_inh_size * sizeof(FILE *));
    //file_init("syn_inh_R", p_syn_inh_size, f_syn_inh_R);
    //FILE **f_syn_inh_S = malloc(p_syn_inh_size * sizeof(FILE *));
    //file_init("syn_inh_S", p_syn_inh_size, f_syn_inh_S);


    // Begin simulation loop
    // =====================

    run_time = clock();

    for (i = 0; i < p_sim_steps; i++) {
        if (i % p_report_step == 0) printf("%lu\n", i);

        // Event(s) occur if urandom >= e^-rate
        for (j = 0; j < p_in1_size; j++) {
            mpfr_urandom(rand_uf, rng_uf, MPFR_RNDN);
            in1[j] = mpfr_greaterequal_p(rand_uf, in1_p0);
            fprintf(stderr, "%s", (in1[j] ? "\x1B[31m\xE2\x96\xA3\x1B[0m" : "\xE2\x96\xA3"));
        }
        fprintf(stderr, "  ");
        for (j = 0; j < p_in2_size; j++) {
            mpfr_urandom(rand_uf, rng_uf, MPFR_RNDN);
            in2[j] = mpfr_greaterequal_p(rand_uf, in2_p0);
            fprintf(stderr, "%s", (in2[j] ? "\x1B[31m\xE2\x96\xA3\x1B[0m" : "\xE2\x96\xA3"));
        }
        fprintf(stderr, "\n");

        // Compute derivatives
        for (j = 0; j < p_nrn1_size; j++) {
            dNdt(j, 1);
            dVdt(j, 1);
        }
        for (j = 0; j < p_nrn2_size; j++) {
            dNdt(j, 2);
            dVdt(j, 2);
        }
        for (j = 0; j < p_syn_exc_size; j++) {
            dRdt(j, 1);
            dSdt(j, 1);
        }
        for (j = 0; j < p_syn_inh_size; j++) {
            dRdt(j, 2);
            dSdt(j, 2);
        }

        // Step system
        mpfr_add(t, t, h, MPFR_RNDN);
        for (j = 0; j < p_nrn1_size; j++) {
            mpfr_mul(&(d_nrn1_N[j]), &(d_nrn1_N[j]), h, MPFR_RNDN);
            mpfr_add(&(nrn1_N[j]), &(nrn1_N[j]), &(d_nrn1_N[j]), MPFR_RNDN);
            mpfr_mul(&(d_nrn1_V[j]), &(d_nrn1_V[j]), h, MPFR_RNDN);
            mpfr_add(&(nrn1_V[j]), &(nrn1_V[j]), &(d_nrn1_V[j]), MPFR_RNDN);
        }
        for (j = 0; j < p_nrn2_size; j++) {
            mpfr_mul(&(d_nrn2_N[j]), &(d_nrn2_N[j]), h, MPFR_RNDN);
            mpfr_add(&(nrn2_N[j]), &(nrn2_N[j]), &(d_nrn2_N[j]), MPFR_RNDN);
            mpfr_mul(&(d_nrn2_V[j]), &(d_nrn2_V[j]), h, MPFR_RNDN);
            mpfr_add(&(nrn2_V[j]), &(nrn2_V[j]), &(d_nrn2_V[j]), MPFR_RNDN);
        }
        for (j = 0; j < p_syn_exc_size; j++) {
            mpfr_mul(&(d_syn_exc_R[j]), &(d_syn_exc_R[j]), h, MPFR_RNDN);
            mpfr_add(&(syn_exc_R[j]), &(syn_exc_R[j]), &(d_syn_exc_R[j]), MPFR_RNDN);
            mpfr_mul(&(d_syn_exc_S[j]), &(d_syn_exc_S[j]), h, MPFR_RNDN);
            mpfr_add(&(syn_exc_S[j]), &(syn_exc_S[j]), &(d_syn_exc_S[j]), MPFR_RNDN);
        }
        for (j = 0; j < p_syn_inh_size; j++) {
            mpfr_mul(&(d_syn_inh_R[j]), &(d_syn_inh_R[j]), h, MPFR_RNDN);
            mpfr_add(&(syn_inh_R[j]), &(syn_inh_R[j]), &(d_syn_inh_R[j]), MPFR_RNDN);
            mpfr_mul(&(d_syn_inh_S[j]), &(d_syn_inh_S[j]), h, MPFR_RNDN);
            mpfr_add(&(syn_inh_S[j]), &(syn_inh_S[j]), &(d_syn_inh_S[j]), MPFR_RNDN);
        }

        file_write(t, 1, f_time);

        file_write(nrn1_N, p_nrn1_size, f_nrn1_N);
        file_write(nrn1_V, p_nrn1_size, f_nrn1_V);

        //file_write(nrn2_N, p_nrn2_size, f_nrn2_N);
        //file_write(nrn2_V, p_nrn2_size, f_nrn2_V);

        file_write(syn_exc_R, p_syn_exc_size, f_syn_exc_R);
        file_write(syn_exc_S, p_syn_exc_size, f_syn_exc_S);

        //file_write(syn_inh_R, p_syn_inh_size, f_syn_inh_R);
        //file_write(syn_inh_S, p_syn_inh_size, f_syn_inh_S);
    }

    run_time = clock() - run_time;
    printf("Finished in %f seconds.\n", ((float) run_time) / CLOCKS_PER_SEC);

    // End simulation loop
    // ===================


    // Clear system state
    mpfr_clear(h);
    mpfr_clear(t);
    for (i = 0; i < p_nrn1_size; i++) {
        mpfr_clear(&(nrn1_N[i]));
        mpfr_clear(&(nrn1_V[i]));
        mpfr_clear(&(d_nrn1_N[i]));
        mpfr_clear(&(d_nrn1_V[i]));
    }
    for (i = 0; i < p_nrn2_size; i++) {
        mpfr_clear(&(nrn2_N[i]));
        mpfr_clear(&(nrn2_V[i]));
        mpfr_clear(&(d_nrn2_N[i]));
        mpfr_clear(&(d_nrn2_V[i]));
    }
    for (i = 0; i < p_syn_exc_size; i++) {
        mpfr_clear(&(syn_exc_R[i]));
        mpfr_clear(&(syn_exc_S[i]));
        mpfr_clear(&(d_syn_exc_R[i]));
        mpfr_clear(&(d_syn_exc_S[i]));
    }
    for (i = 0; i < p_syn_inh_size; i++) {
        mpfr_clear(&(syn_inh_R[i]));
        mpfr_clear(&(syn_inh_S[i]));
        mpfr_clear(&(d_syn_inh_R[i]));
        mpfr_clear(&(d_syn_inh_S[i]));
    }

    // Clear Poisson input parameters (group 1)
    mpfr_clear(in1_p0);
    mpfr_clear(in1_V_lo);
    mpfr_clear(in1_V_hi);

    // Clear Poisson input parameters (group 2)
    mpfr_clear(in2_p0);
    mpfr_clear(in2_V_lo);
    mpfr_clear(in2_V_hi);

    // Clear neuron parameters
    mpfr_clear(GL);
    mpfr_clear(GCa);
    mpfr_clear(GK);
    mpfr_clear(VL);
    mpfr_clear(VCa);
    mpfr_clear(VK);
    mpfr_clear(V1);
    mpfr_clear(V2);
    mpfr_clear(V3);
    mpfr_clear(V4);
    mpfr_clear(phi);
    mpfr_clear(C);

    // Clear excitatory synapse parameters
    for (i = 0; i < p_syn_exc_size; i++) {
        mpfr_clear(&(syn_exc_GSyn[i]));
    }
    mpfr_clear(syn_exc_VSyn);
    mpfr_clear(syn_exc_thr);
    mpfr_clear(syn_exc_a);
    mpfr_clear(syn_exc_b);
    mpfr_clear(syn_exc_k);

    // Clear inhibitory synapse parameters
    for (i = 0; i < p_syn_inh_size; i++) {
        mpfr_clear(&(syn_inh_GSyn[i]));
    }
    mpfr_clear(syn_inh_VSyn);
    mpfr_clear(syn_inh_thr);
    mpfr_clear(syn_inh_a);
    mpfr_clear(syn_inh_b);
    mpfr_clear(syn_inh_k);

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
    for (i = 0; i < p_in1_size; i++) {
        mpfr_clear(&(I1[i]));
    }
    for (i = 0; i < p_in2_size; i++) {
        mpfr_clear(&(I2[i]));
    }

    // Free dynamic arrays
    free(nrn1_N);
    free(nrn2_N);
    free(nrn1_V);
    free(nrn2_V);
    free(syn_exc_R);
    free(syn_inh_R);
    free(syn_exc_S);
    free(syn_inh_S);
    free(d_nrn1_N);
    free(d_nrn2_N);
    free(d_nrn1_V);
    free(d_nrn2_V);
    free(d_syn_exc_R);
    free(d_syn_inh_R);
    free(d_syn_exc_S);
    free(d_syn_inh_S);
    free(syn_exc_GSyn);
    free(syn_inh_GSyn);
    free(I1);
    free(I2);
    free(in1);
    free(in2);

    // Clear report files
    file_clear(1, f_time);
    free(f_time);

    file_clear(p_nrn1_size, f_nrn1_N);
    free(f_nrn1_N);
    file_clear(p_nrn1_size, f_nrn1_V);
    free(f_nrn1_V);

    //file_clear(p_nrn2_size, f_nrn2_N);
    //free(f_nrn2_N);
    //file_clear(p_nrn2_size, f_nrn2_V);
    //free(f_nrn2_V);

    file_clear(p_syn_exc_size, f_syn_exc_R);
    free(f_syn_exc_R);
    file_clear(p_syn_exc_size, f_syn_exc_S);
    free(f_syn_exc_S);

    //file_clear(p_syn_inh_size, f_syn_inh_R);
    //free(f_syn_inh_R);
    //file_clear(p_syn_inh_size, f_syn_inh_S);
    //free(f_syn_inh_S);

    gmp_randclear(rng_uf);
    gmp_randclear(rng_nf);
    gmp_randclear(rng_uz);
    mpfr_free_cache();

    return 0;
}
