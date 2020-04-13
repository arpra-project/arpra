/*
 * morris_lecar.c -- Test Morris-Lecar model.
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
#include <arpra_ode.h>

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
#define p_prec_internal 256
#define p_sim_steps 1000
#define p_report_step 20
#define p_reduce_step 100
#define p_reduce_rel 0.3

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


// DEBUG: print MPFR numbers to stderr
void debug (mpfr_srcptr x) {
    mpfr_out_str(stderr, 10, 80, x, MPFR_RNDN);
    fputs("\n", stderr);
}

void file_init (char *grp, arpra_uint grp_size,
                FILE **c, FILE **r, FILE **n, FILE **s, FILE **d)
{
    char fname[20];
    arpra_uint i;

    for (i = 0; i < grp_size; i++) {
        sprintf(fname, "%s_%03lu_c.dat", grp, i);
        c[i] = fopen(fname, "w");
        sprintf(fname, "%s_%03lu_r.dat", grp, i);
        r[i] = fopen(fname, "w");
        sprintf(fname, "%s_%03lu_n.dat", grp, i);
        n[i] = fopen(fname, "w");
        sprintf(fname, "%s_%03lu_s.dat", grp, i);
        s[i] = fopen(fname, "w");
        sprintf(fname, "%s_%03lu_d.dat", grp, i);
        d[i] = fopen(fname, "w");
    }
}

void file_clear (arpra_uint grp_size, FILE **c, FILE **r, FILE **n, FILE **s, FILE **d)
{
    arpra_uint i;

    for (i = 0; i < grp_size; i++) {
        fclose(c[i]);
        fclose(r[i]);
        fclose(n[i]);
        fclose(s[i]);
        fclose(d[i]);
    }
}

void file_write (const arpra_range *A, arpra_uint grp_size,
                 FILE **c, FILE **r, FILE **n, FILE **s, FILE **d)
{
    arpra_uint i, j;

    for (i = 0; i < grp_size; i++) {
        mpfr_out_str(c[i], 10, 80, &(A[i].centre), MPFR_RNDN);
        fputc('\n', c[i]);
        mpfr_out_str(r[i], 10, 80, &(A[i].radius), MPFR_RNDN);
        fputc('\n', r[i]);
        fprintf(n[i], "%lu\n", A[i].nTerms);
        for (j = 0; j < A[i].nTerms; j++) {
            fprintf(s[i], "%lu ", A[i].symbols[j]);
            mpfr_out_str(d[i], 10, 80, &(A[i].deviations[j]), MPFR_RNDN);
            fputc(' ', d[i]);
        }
        fputc('\n', s[i]);
        fputc('\n', d[i]);
    }
}

struct dNdt_params
{
    arpra_uint grp_V;
    arpra_range *V3;
    arpra_range *V4;
    arpra_range *phi;
    arpra_range *one;
    arpra_range *two;
    arpra_range *neg_two;
    arpra_range *N_ss;
    arpra_range *temp1;
    arpra_range *temp2;
};

void dNdt (arpra_range *y, const void *params,
           const arpra_range *t, const arpra_range **x,
           const arpra_uint x_grp, const arpra_uint x_dim)
{
    const struct dNdt_params *p = (struct dNdt_params *) params;
    const arpra_range *N = &(x[x_grp][x_dim]);
    const arpra_range *V = &(x[p->grp_V][x_dim]);
    const arpra_range *V3 = p->V3;
    const arpra_range *V4 = p->V4;
    const arpra_range *phi = p->phi;
    const arpra_range *one = p->one;
    const arpra_range *two = p->two;
    const arpra_range *neg_two = p->neg_two;
    arpra_range *N_ss = p->N_ss;
    arpra_range *temp1 = p->temp1;
    arpra_range *temp2 = p->temp2;

    // K+ channel activation steady-state
    // N_ss = 1 / (1 + exp(-2 (V - V3) / V4))
    arpra_sub(temp1, V, V3);
    arpra_mul(N_ss, neg_two, temp1);
    arpra_div(N_ss, N_ss, V4);
    arpra_exp(N_ss, N_ss);
    arpra_add(N_ss, one, N_ss);
    arpra_inv(N_ss, N_ss);

    // tau of K+ channel activation
    // tau = 1 / (phi ((p + q) / 2))
    // p = exp(-(V - V3) / (2 V4))
    // q = exp( (V - V3) / (2 V4))
    arpra_mul(temp2, two, V4);
    arpra_div(temp2, temp1, temp2);
    arpra_neg(temp1, temp2);
    arpra_exp(temp1, temp1);
    arpra_exp(temp2, temp2);
    arpra_add(temp1, temp1, temp2);
    arpra_div(temp1, temp1, two);
    arpra_mul(temp1, phi, temp1);
    arpra_inv(temp1, temp1);

    // delta of K+ channel activation
    // dN/dt = (N_ss - N) / tau
    arpra_sub(y, N_ss, N);
    arpra_div(y, y, temp1);
}

struct dVdt_params
{
    arpra_uint grp_N;
    arpra_uint grp_S;
    arpra_range *GSyn;
    arpra_range *VSyn;
    arpra_range *GL;
    arpra_range *VL;
    arpra_range *GCa;
    arpra_range *VCa;
    arpra_range *GK;
    arpra_range *VK;
    arpra_range *V1;
    arpra_range *V2;
    arpra_range *C;
    arpra_uint pre_syn_size;
    arpra_range *one;
    arpra_range *neg_two;
    arpra_range *I;
    arpra_range *M_ss;
    arpra_range *temp1;
};

void dVdt (arpra_range *y, const void *params,
           const arpra_range *t, const arpra_range **x,
           const arpra_uint x_grp, const arpra_uint x_dim)
{
    const struct dVdt_params *p = (struct dVdt_params *) params;
    const arpra_range *V = &(x[x_grp][x_dim]);
    const arpra_range *N = &(x[p->grp_N][x_dim]);
    const arpra_range *S = &(x[p->grp_S][x_dim * p->pre_syn_size]);
    const arpra_range *GSyn = &(p->GSyn[x_dim * p->pre_syn_size]);
    const arpra_range *VSyn = p->VSyn;
    const arpra_range *GL = p->GL;
    const arpra_range *VL = p->VL;
    const arpra_range *GCa = p->GCa;
    const arpra_range *VCa = p->VCa;
    const arpra_range *GK = p->GK;
    const arpra_range *VK = p->VK;
    const arpra_range *V1 = p->V1;
    const arpra_range *V2 = p->V2;
    const arpra_range *C = p->C;
    const arpra_range *one = p->one;
    const arpra_range *neg_two = p->neg_two;
    arpra_range *I = p->I;
    arpra_range *M_ss = p->M_ss;
    arpra_range *temp1 = p->temp1;
    arpra_uint i;

    // Ca++ channel activation steady-state
    // M_ss = 1 / (1 + exp(-2 (V - V1) / V2))
    arpra_sub(M_ss, V, V1);
    arpra_mul(M_ss, neg_two, M_ss);
    arpra_div(M_ss, M_ss, V2);
    arpra_exp(M_ss, M_ss);
    arpra_add(M_ss, one, M_ss);
    arpra_inv(M_ss, M_ss);

    // Synapse current
    arpra_sub(temp1, VSyn, V);
    for (i = 0; i < p->pre_syn_size; i++) {
        arpra_mul(&(I[i]), temp1, &(GSyn[i]));
        arpra_mul(&(I[i]), &(I[i]), &(S[i]));
    }
    arpra_sum_recursive(y, I, p->pre_syn_size);
    //arpra_sum(y, I, p->pre_syn_size);




    // ======== TEMP DEBUG ==========
    //arpra_set_d(y, 80.0); // bifurcation at sum(I) = 80.0
    if (x_dim == 0) {
        //fprintf(stderr, "sum(I): "); debug(y->centre);
        //fprintf(stderr, "I[0]: "); debug(I[0].centre);
        for (i = 0; i < p->pre_syn_size; i++) {
            //fprintf(stderr, "I[%lu]: ", i); debug(I[i].centre);
        }
    }



    // Leak current
    arpra_sub(temp1, V, VL);
    arpra_mul(temp1, temp1, GL);
    arpra_sub(y, y, temp1);

    // Ca++ current
    arpra_sub(temp1, V, VCa);
    arpra_mul(temp1, temp1, GCa);
    arpra_mul(temp1, temp1, M_ss);
    arpra_sub(y, y, temp1);

    // K+ current
    arpra_sub(temp1, V, VK);
    arpra_mul(temp1, temp1, GK);
    arpra_mul(temp1, temp1, N);
    arpra_sub(y, y, temp1);

    // delta of membrane potential
    // dV/dt = (I + GL (VL - V) + GCa M (VCa - V) + GK N (VK - V)) / C
    arpra_div(y, y, C);
}

struct dRdt_params
{
    arpra_range *a;
    arpra_range *b;
    arpra_range *k;
    arpra_range *VPre_lo;
    arpra_range *VPre_hi;
    arpra_range *threshold;
    int *in;
    arpra_uint pre_syn_size;
    arpra_range *one;
    arpra_range *temp1;
    arpra_range *temp2;
};

void dRdt (arpra_range *y, const void *params,
           const arpra_range *t, const arpra_range **x,
           const arpra_uint x_grp, const arpra_uint x_dim)
{
    const struct dRdt_params *p = (struct dRdt_params *) params;
    const arpra_range *R = &(x[x_grp][x_dim]);
    const arpra_range *a = p->a;
    const arpra_range *b = p->b;
    const arpra_range *k = p->k;
    const arpra_range *VPre = p->in[x_dim % p->pre_syn_size] ? p->VPre_hi : p->VPre_lo;
    const arpra_range *threshold = p->threshold;
    const arpra_range *one = p->one;
    arpra_range *temp1 = p->temp1;
    arpra_range *temp2 = p->temp2;

    // Sigmoid of threshold difference
    arpra_sub(temp1, VPre, threshold);
    arpra_mul(temp1, temp1, k);
    arpra_exp(temp1, temp1);
    arpra_add(temp1, temp1, one);
    arpra_inv(temp1, temp1);

    // Presynaptic transmitter release rise
    arpra_mul(temp1, a, temp1);

    // Presynaptic transmitter release decay
    arpra_mul(temp2, b, R);

    // delta of presynaptic transmitter release
    // dR/dt = a Q - b R
    // Q = 1 / (1 + e^(k(V - threshold)))
    arpra_sub(y, temp1, temp2);
}

struct dSdt_params
{
    arpra_uint grp_R;
    arpra_range *a;
    arpra_range *b;
    arpra_range *temp1;
    arpra_range *temp2;
};

void dSdt (arpra_range *y, const void *params,
           const arpra_range *t, const arpra_range **x,
           const arpra_uint x_grp, const arpra_uint x_dim)
{
    const struct dSdt_params *p = (struct dSdt_params *) params;
    const arpra_range *S = &(x[x_grp][x_dim]);
    const arpra_range *R = &(x[p->grp_R][x_dim]);
    const arpra_range *a = p->a;
    const arpra_range *b = p->b;
    arpra_range *temp1 = p->temp1;
    arpra_range *temp2 = p->temp2;

    // Postsynaptic transmitter binding rise
    arpra_mul(temp1, a, R);

    // Postsynaptic transmitter binding decay
    arpra_mul(temp2, b, S);

    // delta of postsynaptic transmitter binding
    // dS/dt = a R - b S
    arpra_sub(y, temp1, temp2);
}

int main (int argc, char *argv[])
{
    arpra_uint i, j;
    arpra_range h, sys_t;

    enum grps {
        grp_nrn_N, grp_nrn_V,
        grp_syn_R, grp_syn_S,
    };

    arpra_set_range_method(ARPRA_MIXED_TRIMMED_IAAA);
    arpra_set_internal_precision(p_prec_internal);

    // Initialise arpra_reduce_small_rel threshold.
    mpfr_t reduce_rel;
    mpfr_init2(reduce_rel, 53);
    mpfr_set_d(reduce_rel, p_reduce_rel, MPFR_RNDN);

    // Allocate system state
    arpra_range *nrn_N = malloc(p_nrn_size * sizeof(arpra_range));
    arpra_range *nrn_V = malloc(p_nrn_size * sizeof(arpra_range));
    arpra_range *syn_R = malloc(p_syn_size * sizeof(arpra_range));
    arpra_range *syn_S = malloc(p_syn_size * sizeof(arpra_range));

    // Allocate other arrays
    arpra_range *syn_GSyn = malloc(p_syn_size * sizeof(arpra_range));
    arpra_range *I = malloc(p_in_size * sizeof(arpra_range));
    int *in = malloc(p_in_size * sizeof(int));
    arpra_uint *nrn_N_reduce_epoch = malloc(p_nrn_size * sizeof(arpra_uint));
    arpra_uint *nrn_V_reduce_epoch = malloc(p_nrn_size * sizeof(arpra_uint));
    arpra_uint *syn_R_reduce_epoch = malloc(p_syn_size * sizeof(arpra_uint));
    arpra_uint *syn_S_reduce_epoch = malloc(p_syn_size * sizeof(arpra_uint));

    mpfr_t in_p0, rand_uf, rand_nf;
    arpra_range nrn_GL, nrn_VL, nrn_GCa, nrn_VCa, nrn_GK, nrn_VK, nrn_V1, nrn_V2,
        nrn_V3, nrn_V4, nrn_phi, nrn_C, syn_VSyn, syn_thr, syn_a, syn_b, syn_k,
        one, two, neg_two, temp1, temp2, M_ss, N_ss, in_V_lo, in_V_hi;

    struct timespec clock_time;

    // Initialise uniform float RNG
    gmp_randstate_t rng_uf;
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
    gmp_randstate_t rng_nf;
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
    arpra_init2(&h, p_prec);
    arpra_init2(&sys_t, p_prec);
    for (i = 0; i < p_nrn_size; i++) {
        arpra_init2(&(nrn_N[i]), p_prec);
        arpra_init2(&(nrn_V[i]), p_prec);
    }
    for (i = 0; i < p_syn_size; i++) {
        arpra_init2(&(syn_R[i]), p_prec);
        arpra_init2(&(syn_S[i]), p_prec);
    }

    // Initialise Poisson input parameters
    mpfr_init2(in_p0, p_prec);
    arpra_init2(&in_V_lo, p_prec);
    arpra_init2(&in_V_hi, p_prec);

    // Initialise neuron parameters
    arpra_init2(&nrn_GL, p_prec);
    arpra_init2(&nrn_GCa, p_prec);
    arpra_init2(&nrn_GK, p_prec);
    arpra_init2(&nrn_VL, p_prec);
    arpra_init2(&nrn_VCa, p_prec);
    arpra_init2(&nrn_VK, p_prec);
    arpra_init2(&nrn_V1, p_prec);
    arpra_init2(&nrn_V2, p_prec);
    arpra_init2(&nrn_V3, p_prec);
    arpra_init2(&nrn_V4, p_prec);
    arpra_init2(&nrn_phi, p_prec);
    arpra_init2(&nrn_C, p_prec);

    // Initialise synapse parameters
    for (i = 0; i < p_syn_size; i++) {
        arpra_init2(&(syn_GSyn[i]), p_prec);
    }
    arpra_init2(&syn_VSyn, p_prec);
    arpra_init2(&syn_thr, p_prec);
    arpra_init2(&syn_a, p_prec);
    arpra_init2(&syn_b, p_prec);
    arpra_init2(&syn_k, p_prec);

    // Initialise constants
    arpra_init2(&one, p_prec);
    arpra_init2(&two, p_prec);
    arpra_init2(&neg_two, p_prec);

    // Initialise scratch space
    mpfr_init2(rand_uf, p_rand_prec);
    mpfr_init2(rand_nf, p_rand_prec);
    arpra_init2(&temp1, p_prec);
    arpra_init2(&temp2, p_prec);
    arpra_init2(&M_ss, p_prec);
    arpra_init2(&N_ss, p_prec);
    for (i = 0; i < p_in_size; i++) {
        arpra_init2(&(I[i]), p_prec);
    }

    // Set system state
    arpra_set_d(&h, p_h);
    arpra_set_d(&sys_t, p_t0);
    for (i = 0; i < p_nrn_size; i++) {
        arpra_set_d(&(nrn_N[i]), p_nrn_N0);
        arpra_set_d(&(nrn_V[i]), p_nrn_V0);
    }
    for (i = 0; i < p_syn_size; i++) {
        arpra_set_d(&(syn_R[i]), p_syn_R0);
        arpra_set_d(&(syn_S[i]), p_syn_S0);
    }

    // Set Poisson input parameters
    mpfr_set_d(in_p0, -(p_in_freq / 1000) * p_h, MPFR_RNDN);
    mpfr_exp(in_p0, in_p0, MPFR_RNDN);
    arpra_set_d(&in_V_lo, p_in_V_lo);
    arpra_set_d(&in_V_hi, p_in_V_hi);

    // Set neuron parameters
    arpra_set_d(&nrn_GL, p_nrn_GL);
    arpra_set_d(&nrn_GCa, p_nrn_GCa);
    arpra_set_d(&nrn_GK, p_nrn_GK);
    arpra_set_d(&nrn_VL, p_nrn_VL);
    arpra_set_d(&nrn_VCa, p_nrn_VCa);
    arpra_set_d(&nrn_VK, p_nrn_VK);
    arpra_set_d(&nrn_V1, p_nrn_V1);
    arpra_set_d(&nrn_V2, p_nrn_V2);
    arpra_set_d(&nrn_V3, p_nrn_V3);
    arpra_set_d(&nrn_V4, p_nrn_V4);
    arpra_set_d(&nrn_phi, p_nrn_phi);
    arpra_set_d(&nrn_C, p_nrn_C);

    // Set synapse parameters
    for (i = 0; i < p_syn_size; i++) {
        //mpfr_nrandom(rand_nf, rng_nf, MPFR_RNDN);
        mpfr_grandom(rand_nf, NULL, rng_nf, MPFR_RNDN);
        mpfr_mul_d(rand_nf, rand_nf, p_syn_GSyn_std, MPFR_RNDN);
        mpfr_add_d(rand_nf, rand_nf, p_syn_GSyn_mean, MPFR_RNDN);
        arpra_set_mpfr(&(syn_GSyn[i]), rand_nf);
    }
    arpra_set_d(&syn_VSyn, p_syn_VSyn);
    arpra_set_d(&syn_thr, p_syn_thr);
    arpra_set_d(&syn_a, p_syn_a);
    arpra_set_d(&syn_b, p_syn_b);
    arpra_set_d(&syn_k, p_syn_k);
    arpra_neg(&syn_k, &syn_k);

    // Set constants
    arpra_set_d(&one, 1.0);
    arpra_set_d(&two, 2.0);
    arpra_set_d(&neg_two, -2.0);

    // Initialise report files
    FILE **f_time_c = malloc(sizeof(FILE *));
    FILE **f_time_r = malloc(sizeof(FILE *));
    FILE **f_time_n = malloc(sizeof(FILE *));
    FILE **f_time_s = malloc(sizeof(FILE *));
    FILE **f_time_d = malloc(sizeof(FILE *));
    file_init("time", 1, f_time_c, f_time_r, f_time_n, f_time_s, f_time_d);

    FILE **f_nrn_N_c = malloc(p_nrn_size * sizeof(FILE *));
    FILE **f_nrn_N_r = malloc(p_nrn_size * sizeof(FILE *));
    FILE **f_nrn_N_n = malloc(p_nrn_size * sizeof(FILE *));
    FILE **f_nrn_N_s = malloc(p_nrn_size * sizeof(FILE *));
    FILE **f_nrn_N_d = malloc(p_nrn_size * sizeof(FILE *));
    file_init("nrn_N", p_nrn_size, f_nrn_N_c, f_nrn_N_r, f_nrn_N_n, f_nrn_N_s, f_nrn_N_d);
    FILE **f_nrn_V_c = malloc(p_nrn_size * sizeof(FILE *));
    FILE **f_nrn_V_r = malloc(p_nrn_size * sizeof(FILE *));
    FILE **f_nrn_V_n = malloc(p_nrn_size * sizeof(FILE *));
    FILE **f_nrn_V_s = malloc(p_nrn_size * sizeof(FILE *));
    FILE **f_nrn_V_d = malloc(p_nrn_size * sizeof(FILE *));
    file_init("nrn_V", p_nrn_size, f_nrn_V_c, f_nrn_V_r, f_nrn_V_n, f_nrn_V_s, f_nrn_V_d);

    /* FILE **f_syn_R_c = malloc(p_syn_size * sizeof(FILE *)); */
    /* FILE **f_syn_R_r = malloc(p_syn_size * sizeof(FILE *)); */
    /* FILE **f_syn_R_n = malloc(p_syn_size * sizeof(FILE *)); */
    /* FILE **f_syn_R_s = malloc(p_syn_size * sizeof(FILE *)); */
    /* FILE **f_syn_R_d = malloc(p_syn_size * sizeof(FILE *)); */
    /* file_init("syn_R", p_syn_size, f_syn_R_c, f_syn_R_r, f_syn_R_n, f_syn_R_s, f_syn_R_d); */
    /* FILE **f_syn_S_c = malloc(p_syn_size * sizeof(FILE *)); */
    /* FILE **f_syn_S_r = malloc(p_syn_size * sizeof(FILE *)); */
    /* FILE **f_syn_S_n = malloc(p_syn_size * sizeof(FILE *)); */
    /* FILE **f_syn_S_s = malloc(p_syn_size * sizeof(FILE *)); */
    /* FILE **f_syn_S_d = malloc(p_syn_size * sizeof(FILE *)); */
    /* file_init("syn_S", p_syn_size, f_syn_S_c, f_syn_S_r, f_syn_S_n, f_syn_S_s, f_syn_S_d); */

    // Set parameter structs
    struct dNdt_params params_nrn_N = {
        .grp_V = grp_nrn_V,
        .V3 = &nrn_V3,
        .V4 = &nrn_V4,
        .phi = &nrn_phi,
        .one = &one,
        .two = &two,
        .neg_two = &neg_two,
        .N_ss = &N_ss,
        .temp1 = &temp1,
        .temp2 = &temp2,
    };

    struct dVdt_params params_nrn_V = {
        .grp_N = grp_nrn_N,
        .grp_S = grp_syn_S,
        .GSyn = syn_GSyn,
        .VSyn = &syn_VSyn,
        .GL = &nrn_GL,
        .VL = &nrn_VL,
        .GCa = &nrn_GCa,
        .VCa = &nrn_VCa,
        .GK = &nrn_GK,
        .VK = &nrn_VK,
        .V1 = &nrn_V1,
        .V2 = &nrn_V2,
        .C = &nrn_C,
        .pre_syn_size = p_in_size,
        .one = &one,
        .neg_two = &neg_two,
        .I = I,
        .M_ss = &M_ss,
        .temp1 = &temp1,
    };

    struct dRdt_params params_syn_R = {
        .a = &syn_a,
        .b = &syn_b,
        .k = &syn_k,
        .VPre_lo = &in_V_lo,
        .VPre_hi = &in_V_hi,
        .threshold = &syn_thr,
        .in = in,
        .pre_syn_size = p_in_size,
        .one = &one,
        .temp1 = &temp1,
        .temp2 = &temp2,
    };

    struct dSdt_params params_syn_S = {
        .grp_R = grp_syn_R,
        .a = &syn_a,
        .b = &syn_b,
        .temp1 = &temp1,
        .temp2 = &temp2,
    };

    // ODE system
    arpra_uint sys_grps = 4;
    arpra_uint sys_dims[4] = {
        p_nrn_size, p_nrn_size,
        p_syn_size, p_syn_size,
    };
    arpra_ode_f sys_f[4] = {
        dNdt, dVdt, dRdt, dSdt,
    };
    void *sys_params[4] = {
        &params_nrn_N, &params_nrn_V,
        &params_syn_R, &params_syn_S,
    };
    arpra_range *sys_x[4] = {
        nrn_N, nrn_V, syn_R, syn_S,
    };
    arpra_ode_system ode_system = {
        .f = sys_f,
        .params = sys_params,
        .t = &sys_t,
        .x = sys_x,
        .grps = sys_grps,
        .dims = sys_dims,
    };

    // ODE stepper
    arpra_ode_stepper ode_stepper;
    arpra_ode_stepper_init(&ode_stepper, &ode_system, arpra_ode_euler);
    //arpra_ode_stepper_init(&ode_stepper, &ode_system, arpra_ode_trapezoidal);
    //arpra_ode_stepper_init(&ode_stepper, &ode_system, arpra_ode_bogsham32);
    //arpra_ode_stepper_init(&ode_stepper, &ode_system, arpra_ode_dopri54);
    //arpra_ode_stepper_init(&ode_stepper, &ode_system, arpra_ode_dopri87);


    // Begin simulation loop
    // =====================

    clock_t run_time = clock();

    for (i = 0; i < p_sim_steps; i++) {
        if (i % p_report_step == 0) printf("%lu\n", i);

        for (j = 0; j < p_nrn_size; j++) {
            nrn_N_reduce_epoch[j] = ode_system.x[grp_nrn_N][j].nTerms;
            nrn_V_reduce_epoch[j] = ode_system.x[grp_nrn_V][j].nTerms;
        }
        for (j = 0; j < p_syn_size; j++) {
            syn_R_reduce_epoch[j] = ode_system.x[grp_syn_R][j].nTerms;
            syn_S_reduce_epoch[j] = ode_system.x[grp_syn_S][j].nTerms;
        }

        // Event(s) occur if urandom >= e^-rate
        for (j = 0; j < p_in_size; j++) {
            mpfr_urandom(rand_uf, rng_uf, MPFR_RNDN);
            in[j] = mpfr_greaterequal_p(rand_uf, in_p0);
            fprintf(stderr, "%s", (in[j] ? "\x1B[31m\xE2\x96\xA3\x1B[0m" : "\xE2\x96\xA3"));
        }
        fprintf(stderr, "\n");

        // Step system
        arpra_ode_stepper_step(&ode_stepper, &h);

        arpra_uint reduce_n;
        for (j = 0; j < p_nrn_size; j++) {
            reduce_n = ode_system.x[grp_nrn_N][j].nTerms - nrn_N_reduce_epoch[j];
            arpra_reduce_last_n(&(ode_system.x[grp_nrn_N][j]), &(ode_system.x[grp_nrn_N][j]), reduce_n);
            reduce_n = ode_system.x[grp_nrn_V][j].nTerms - nrn_V_reduce_epoch[j];
            arpra_reduce_last_n(&(ode_system.x[grp_nrn_V][j]), &(ode_system.x[grp_nrn_V][j]), reduce_n);
        }
        for (j = 0; j < p_syn_size; j++) {
            reduce_n = ode_system.x[grp_syn_R][j].nTerms - syn_R_reduce_epoch[j];
            arpra_reduce_last_n(&(ode_system.x[grp_syn_R][j]), &(ode_system.x[grp_syn_R][j]), reduce_n);
            reduce_n = ode_system.x[grp_syn_S][j].nTerms - syn_S_reduce_epoch[j];
            arpra_reduce_last_n(&(ode_system.x[grp_syn_S][j]), &(ode_system.x[grp_syn_S][j]), reduce_n);
        }

        if (i % p_reduce_step == 0) {
            for (j = 0; j < p_nrn_size; j++) {
                arpra_reduce_small_rel(&(ode_system.x[grp_nrn_N][j]), &(ode_system.x[grp_nrn_N][j]), reduce_rel);
                arpra_reduce_small_rel(&(ode_system.x[grp_nrn_V][j]), &(ode_system.x[grp_nrn_V][j]), reduce_rel);
            }
            for (j = 0; j < p_syn_size; j++) {
                arpra_reduce_small_rel(&(ode_system.x[grp_syn_R][j]), &(ode_system.x[grp_syn_R][j]), reduce_rel);
                arpra_reduce_small_rel(&(ode_system.x[grp_syn_S][j]), &(ode_system.x[grp_syn_S][j]), reduce_rel);
            }
        }

        file_write(&sys_t, 1, f_time_c, f_time_r, f_time_n, f_time_s, f_time_d);

        file_write(nrn_N, p_nrn_size, f_nrn_N_c, f_nrn_N_r, f_nrn_N_n, f_nrn_N_s, f_nrn_N_d);
        file_write(nrn_V, p_nrn_size, f_nrn_V_c, f_nrn_V_r, f_nrn_V_n, f_nrn_V_s, f_nrn_V_d);

        /* file_write(syn_R, p_syn_size, f_syn_R_c, f_syn_R_r, f_syn_R_n, f_syn_R_s, f_syn_R_d); */
        /* file_write(syn_S, p_syn_size, f_syn_S_c, f_syn_S_r, f_syn_S_n, f_syn_S_s, f_syn_S_d); */
    }

    run_time = clock() - run_time;
    printf("Finished in %f seconds.\n", ((float) run_time) / CLOCKS_PER_SEC);

    // End simulation loop
    // ===================


    // Clear arpra_reduce_small_rel threshold.
    mpfr_clear(reduce_rel);

    // Clear system state
    arpra_clear(&h);
    arpra_clear(&sys_t);
    for (i = 0; i < p_nrn_size; i++) {
        arpra_clear(&(nrn_N[i]));
        arpra_clear(&(nrn_V[i]));
    }
    for (i = 0; i < p_syn_size; i++) {
        arpra_clear(&(syn_R[i]));
        arpra_clear(&(syn_S[i]));
    }

    // Clear Poisson input parameters
    mpfr_clear(in_p0);
    arpra_clear(&in_V_lo);
    arpra_clear(&in_V_hi);

    // Clear neuron parameters
    arpra_clear(&nrn_GL);
    arpra_clear(&nrn_GCa);
    arpra_clear(&nrn_GK);
    arpra_clear(&nrn_VL);
    arpra_clear(&nrn_VCa);
    arpra_clear(&nrn_VK);
    arpra_clear(&nrn_V1);
    arpra_clear(&nrn_V2);
    arpra_clear(&nrn_V3);
    arpra_clear(&nrn_V4);
    arpra_clear(&nrn_phi);
    arpra_clear(&nrn_C);

    // Clear synapse parameters
    for (i = 0; i < p_syn_size; i++) {
        arpra_clear(&(syn_GSyn[i]));
    }
    arpra_clear(&syn_VSyn);
    arpra_clear(&syn_thr);
    arpra_clear(&syn_a);
    arpra_clear(&syn_b);
    arpra_clear(&syn_k);

    // Clear constants
    arpra_clear(&one);
    arpra_clear(&two);
    arpra_clear(&neg_two);

    // Clear scratch space
    mpfr_clear(rand_uf);
    mpfr_clear(rand_nf);
    arpra_clear(&temp1);
    arpra_clear(&temp2);
    arpra_clear(&M_ss);
    arpra_clear(&N_ss);
    for (i = 0; i < p_in_size; i++) {
        arpra_clear(&(I[i]));
    }

    // Free system state
    free(nrn_N);
    free(nrn_V);
    free(syn_R);
    free(syn_S);

    // Free other arrays
    free(syn_GSyn);
    free(I);
    free(in);
    free(nrn_N_reduce_epoch);
    free(nrn_V_reduce_epoch);
    free(syn_R_reduce_epoch);
    free(syn_S_reduce_epoch);

    // Clear report files
    file_clear(1, f_time_c, f_time_r, f_time_n, f_time_s, f_time_d);
    free(f_time_c);
    free(f_time_r);
    free(f_time_n);
    free(f_time_s);
    free(f_time_d);

    file_clear(p_nrn_size, f_nrn_N_c, f_nrn_N_r, f_nrn_N_n, f_nrn_N_s, f_nrn_N_d);
    free(f_nrn_N_c);
    free(f_nrn_N_r);
    free(f_nrn_N_n);
    free(f_nrn_N_s);
    free(f_nrn_N_d);
    file_clear(p_nrn_size, f_nrn_V_c, f_nrn_V_r, f_nrn_V_n, f_nrn_V_s, f_nrn_V_d);
    free(f_nrn_V_c);
    free(f_nrn_V_r);
    free(f_nrn_V_n);
    free(f_nrn_V_s);
    free(f_nrn_V_d);

    /* file_clear(p_syn_size, f_syn_R_c, f_syn_R_r, f_syn_R_n, f_syn_R_s, f_syn_R_d); */
    /* free(f_syn_R_c); */
    /* free(f_syn_R_r); */
    /* free(f_syn_R_n); */
    /* free(f_syn_R_s); */
    /* free(f_syn_R_d); */
    /* file_clear(p_syn_size, f_syn_S_c, f_syn_S_r, f_syn_S_n, f_syn_S_s, f_syn_S_d); */
    /* free(f_syn_S_c); */
    /* free(f_syn_S_r); */
    /* free(f_syn_S_n); */
    /* free(f_syn_S_s); */
    /* free(f_syn_S_d); */

    arpra_ode_stepper_clear(&ode_stepper);
    arpra_clear_buffers();
    gmp_randclear(rng_uf);
    gmp_randclear(rng_nf);
    mpfr_free_cache();

    return 0;
}
