/*
 * morris_lecar.c -- Test Morris-Lecar model.
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
#define p_reduce_ratio 0.3
#define p_prec 53
#define p_prec_internal 2048
#define p_sim_steps 1000
#define p_report_step 20
#define p_reduce_step 50

// RNG parameters
// Seeds are random if not #defined
#define p_rand_prec 53
//#define p_rng_uf_seed 707135875931353ul
//#define p_rng_nf_seed 503108552855933ul

// Poisson input parameters (group 1)
#define p_in1_size 50
#define p_in1_freq 10.0
#define p_in1_V_lo -60.0
#define p_in1_V_hi 20.0

// Poisson input parameters (group 2)
#define p_in2_size 0
#define p_in2_freq 10.0
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


// DEBUG: print MPFR numbers to stderr
void debug (const arpra_mpfr x) {
    mpfr_out_str(stderr, 10, 80, &x, MPFR_RNDN);
    fputs("\n", stderr);
}

void file_init (char *grp, arpra_uint grp_size,
                FILE **c, FILE **r, FILE **n, FILE **s, FILE **d)
{
    char fname[20];
    arpra_uint i;

    for (i = 0; i < grp_size; i++) {
        sprintf(fname, "%s_%03u_c.dat", grp, (unsigned) i);
        c[i] = fopen(fname, "w");
        sprintf(fname, "%s_%03u_r.dat", grp, (unsigned) i);
        r[i] = fopen(fname, "w");
        sprintf(fname, "%s_%03u_n.dat", grp, (unsigned) i);
        n[i] = fopen(fname, "w");
        sprintf(fname, "%s_%03u_s.dat", grp, (unsigned) i);
        s[i] = fopen(fname, "w");
        sprintf(fname, "%s_%03u_d.dat", grp, (unsigned) i);
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

void dNdt (arpra_range *out, const void *params,
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
    arpra_sub(out, N_ss, N);
    arpra_div(out, out, temp1);
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

void dVdt (arpra_range *out, const void *params,
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
    arpra_sum_recursive(out, I, p->pre_syn_size);
    //arpra_sum_exact(out, I, p->pre_syn_size);




    // ======== TEMP DEBUG ==========
    //arpra_set_d(out, 80.0); // bifurcation at sum(I) = 80.0
    if (x_dim == 0) {
        //fprintf(stderr, "sum(I): "); debug(out->centre);
        //fprintf(stderr, "I[0]: "); debug(I[0].centre);
        for (i = 0; i < p->pre_syn_size; i++) {
            //fprintf(stderr, "I[%lu]: ", i); debug(I[i].centre);
        }
    }



    // Leak current
    arpra_sub(temp1, V, VL);
    arpra_mul(temp1, temp1, GL);
    arpra_sub(out, out, temp1);

    // Ca++ current
    arpra_sub(temp1, V, VCa);
    arpra_mul(temp1, temp1, GCa);
    arpra_mul(temp1, temp1, M_ss);
    arpra_sub(out, out, temp1);

    // K+ current
    arpra_sub(temp1, V, VK);
    arpra_mul(temp1, temp1, GK);
    arpra_mul(temp1, temp1, N);
    arpra_sub(out, out, temp1);

    // delta of membrane potential
    // dV/dt = (I + GL (VL - V) + GCa M (VCa - V) + GK N (VK - V)) / C
    arpra_div(out, out, C);
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

void dRdt (arpra_range *out, const void *params,
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
    arpra_sub(out, temp1, temp2);
}

struct dSdt_params
{
    arpra_uint grp_R;
    arpra_range *a;
    arpra_range *b;
    arpra_range *temp1;
    arpra_range *temp2;
};

void dSdt (arpra_range *out, const void *params,
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
    arpra_sub(out, temp1, temp2);
}

int main (int argc, char *argv[])
{
    arpra_uint i, j;
    arpra_range h, sys_t;

    enum grps {
        grp_nrn1_N, grp_nrn1_V, grp_nrn2_N, grp_nrn2_V,
        grp_syn_exc_R, grp_syn_exc_S, grp_syn_inh_R, grp_syn_inh_S
    };

    arpra_set_internal_precision(p_prec_internal);

    // Allocate system state
    arpra_range *nrn1_N = malloc(p_nrn1_size * sizeof(arpra_range));
    arpra_range *nrn1_V = malloc(p_nrn1_size * sizeof(arpra_range));
    arpra_range *nrn2_N = malloc(p_nrn2_size * sizeof(arpra_range));
    arpra_range *nrn2_V = malloc(p_nrn2_size * sizeof(arpra_range));
    arpra_range *syn_exc_R = malloc(p_syn_exc_size * sizeof(arpra_range));
    arpra_range *syn_exc_S = malloc(p_syn_exc_size * sizeof(arpra_range));
    arpra_range *syn_inh_R = malloc(p_syn_inh_size * sizeof(arpra_range));
    arpra_range *syn_inh_S = malloc(p_syn_inh_size * sizeof(arpra_range));

    // Allocate other arrays
    arpra_range *syn_exc_GSyn = malloc(p_syn_exc_size * sizeof(arpra_range));
    arpra_range *syn_inh_GSyn = malloc(p_syn_inh_size * sizeof(arpra_range));
    arpra_range *I1 = malloc(p_in1_size * sizeof(arpra_range));
    arpra_range *I2 = malloc(p_in2_size * sizeof(arpra_range));
    int *in1 = malloc(p_in1_size * sizeof(int));
    int *in2 = malloc(p_in2_size * sizeof(int));
    arpra_uint *nrn1_N_reduce_epoch = malloc(p_nrn1_size * sizeof(arpra_uint));
    arpra_uint *nrn1_V_reduce_epoch = malloc(p_nrn1_size * sizeof(arpra_uint));
    arpra_uint *nrn2_N_reduce_epoch = malloc(p_nrn2_size * sizeof(arpra_uint));
    arpra_uint *nrn2_V_reduce_epoch = malloc(p_nrn2_size * sizeof(arpra_uint));
    arpra_uint *syn_exc_R_reduce_epoch = malloc(p_syn_exc_size * sizeof(arpra_uint));
    arpra_uint *syn_exc_S_reduce_epoch = malloc(p_syn_exc_size * sizeof(arpra_uint));
    arpra_uint *syn_inh_R_reduce_epoch = malloc(p_syn_inh_size * sizeof(arpra_uint));
    arpra_uint *syn_inh_S_reduce_epoch = malloc(p_syn_inh_size * sizeof(arpra_uint));

    arpra_mpfr in1_p0, in2_p0, rand_uf, rand_nf;
    arpra_range GL, VL, GCa, VCa, GK, VK, V1, V2, V3, V4, phi, C,
                syn_exc_VSyn, syn_exc_thr, syn_exc_a, syn_exc_b, syn_exc_k,
                syn_inh_VSyn, syn_inh_thr, syn_inh_a, syn_inh_b, syn_inh_k, one, two, neg_two,
                temp1, temp2, M_ss, N_ss, in1_V_lo, in1_V_hi, in2_V_lo, in2_V_hi;

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
    for (i = 0; i < p_nrn1_size; i++) {
        arpra_init2(&(nrn1_N[i]), p_prec);
        arpra_init2(&(nrn1_V[i]), p_prec);
    }
    for (i = 0; i < p_nrn2_size; i++) {
        arpra_init2(&(nrn2_N[i]), p_prec);
        arpra_init2(&(nrn2_V[i]), p_prec);
    }
    for (i = 0; i < p_syn_exc_size; i++) {
        arpra_init2(&(syn_exc_R[i]), p_prec);
        arpra_init2(&(syn_exc_S[i]), p_prec);
    }
    for (i = 0; i < p_syn_inh_size; i++) {
        arpra_init2(&(syn_inh_R[i]), p_prec);
        arpra_init2(&(syn_inh_S[i]), p_prec);
    }

    // Initialise Poisson input parameters (group 1)
    mpfr_init2(&in1_p0, p_prec);
    arpra_init2(&in1_V_lo, p_prec);
    arpra_init2(&in1_V_hi, p_prec);

    // Initialise Poisson input parameters (group 2)
    mpfr_init2(&in2_p0, p_prec);
    arpra_init2(&in2_V_lo, p_prec);
    arpra_init2(&in2_V_hi, p_prec);

    // Initialise neuron parameters
    arpra_init2(&GL, p_prec);
    arpra_init2(&GCa, p_prec);
    arpra_init2(&GK, p_prec);
    arpra_init2(&VL, p_prec);
    arpra_init2(&VCa, p_prec);
    arpra_init2(&VK, p_prec);
    arpra_init2(&V1, p_prec);
    arpra_init2(&V2, p_prec);
    arpra_init2(&V3, p_prec);
    arpra_init2(&V4, p_prec);
    arpra_init2(&phi, p_prec);
    arpra_init2(&C, p_prec);

    // Initialise excitatory synapse parameters
    for (i = 0; i < p_syn_exc_size; i++) {
        arpra_init2(&(syn_exc_GSyn[i]), p_prec);
    }
    arpra_init2(&syn_exc_VSyn, p_prec);
    arpra_init2(&syn_exc_thr, p_prec);
    arpra_init2(&syn_exc_a, p_prec);
    arpra_init2(&syn_exc_b, p_prec);
    arpra_init2(&syn_exc_k, p_prec);

    // Initialise inhibitory synapse parameters
    for (i = 0; i < p_syn_inh_size; i++) {
        arpra_init2(&(syn_inh_GSyn[i]), p_prec);
    }
    arpra_init2(&syn_inh_VSyn, p_prec);
    arpra_init2(&syn_inh_thr, p_prec);
    arpra_init2(&syn_inh_a, p_prec);
    arpra_init2(&syn_inh_b, p_prec);
    arpra_init2(&syn_inh_k, p_prec);

    // Initialise constants
    arpra_init2(&one, p_prec);
    arpra_init2(&two, p_prec);
    arpra_init2(&neg_two, p_prec);

    // Initialise scratch space
    mpfr_init2(&rand_uf, p_rand_prec);
    mpfr_init2(&rand_nf, p_rand_prec);
    arpra_init2(&temp1, p_prec);
    arpra_init2(&temp2, p_prec);
    arpra_init2(&M_ss, p_prec);
    arpra_init2(&N_ss, p_prec);
    for (i = 0; i < p_in1_size; i++) {
        arpra_init2(&(I1[i]), p_prec);
    }
    for (i = 0; i < p_in2_size; i++) {
        arpra_init2(&(I2[i]), p_prec);
    }

    // Set system state
    arpra_set_d(&h, p_h);
    arpra_set_d(&sys_t, p_t0);
    for (i = 0; i < p_nrn1_size; i++) {
        arpra_set_d(&(nrn1_N[i]), p_nrn1_N0);
        arpra_set_d(&(nrn1_V[i]), p_nrn1_V0);
    }
    for (i = 0; i < p_nrn2_size; i++) {
        arpra_set_d(&(nrn2_N[i]), p_nrn2_N0);
        arpra_set_d(&(nrn2_V[i]), p_nrn2_V0);
    }
    for (i = 0; i < p_syn_exc_size; i++) {
        arpra_set_d(&(syn_exc_R[i]), p_syn_exc_R0);
        arpra_set_d(&(syn_exc_S[i]), p_syn_exc_S0);
    }
    for (i = 0; i < p_syn_inh_size; i++) {
        arpra_set_d(&(syn_inh_R[i]), p_syn_inh_R0);
        arpra_set_d(&(syn_inh_S[i]), p_syn_inh_S0);
    }

    // Set Poisson input parameters (group 1)
    mpfr_set_d(&in1_p0, -(p_in1_freq / 1000) * p_h, MPFR_RNDN);
    mpfr_exp(&in1_p0, &in1_p0, MPFR_RNDN);
    arpra_set_d(&in1_V_lo, p_in1_V_lo);
    arpra_set_d(&in1_V_hi, p_in1_V_hi);

    // Set Poisson input parameters (group 2)
    mpfr_set_d(&in2_p0, -(p_in2_freq / 1000) * p_h, MPFR_RNDN);
    mpfr_exp(&in2_p0, &in2_p0, MPFR_RNDN);
    arpra_set_d(&in2_V_lo, p_in2_V_lo);
    arpra_set_d(&in2_V_hi, p_in2_V_hi);

    // Set neuron parameters
    arpra_set_d(&GL, p_GL);
    arpra_set_d(&GCa, p_GCa);
    arpra_set_d(&GK, p_GK);
    arpra_set_d(&VL, p_VL);
    arpra_set_d(&VCa, p_VCa);
    arpra_set_d(&VK, p_VK);
    arpra_set_d(&V1, p_V1);
    arpra_set_d(&V2, p_V2);
    arpra_set_d(&V3, p_V3);
    arpra_set_d(&V4, p_V4);
    arpra_set_d(&phi, p_phi);
    arpra_set_d(&C, p_C);

    // Set excitatory synapse parameters
    for (i = 0; i < p_syn_exc_size; i++) {
        //mpfr_nrandom(&rand_nf, rng_nf, MPFR_RNDN);
        mpfr_grandom(&rand_nf, NULL, rng_nf, MPFR_RNDN);
        mpfr_mul_d(&rand_nf, &rand_nf, p_syn_exc_GSyn_std, MPFR_RNDN);
        mpfr_add_d(&rand_nf, &rand_nf, p_syn_exc_GSyn_mean, MPFR_RNDN);
        arpra_set_mpfr(&(syn_exc_GSyn[i]), &rand_nf);
    }
    arpra_set_d(&syn_exc_VSyn, p_syn_exc_VSyn);
    arpra_set_d(&syn_exc_thr, p_syn_exc_thr);
    arpra_set_d(&syn_exc_a, p_syn_exc_a);
    arpra_set_d(&syn_exc_b, p_syn_exc_b);
    arpra_set_d(&syn_exc_k, p_syn_exc_k);
    arpra_neg(&syn_exc_k, &syn_exc_k);

    // Set inhibitory synapse parameters
    for (i = 0; i < p_syn_inh_size; i++) {
        //mpfr_nrandom(&rand_nf, rng_nf, MPFR_RNDN);
        mpfr_grandom(&rand_nf, NULL, rng_nf, MPFR_RNDN);
        mpfr_mul_d(&rand_nf, &rand_nf, p_syn_inh_GSyn_std, MPFR_RNDN);
        mpfr_add_d(&rand_nf, &rand_nf, p_syn_inh_GSyn_mean, MPFR_RNDN);
        arpra_set_mpfr(&(syn_inh_GSyn[i]), &rand_nf);
    }
    arpra_set_d(&syn_inh_VSyn, p_syn_inh_VSyn);
    arpra_set_d(&syn_inh_thr, p_syn_inh_thr);
    arpra_set_d(&syn_inh_a, p_syn_inh_a);
    arpra_set_d(&syn_inh_b, p_syn_inh_b);
    arpra_set_d(&syn_inh_k, p_syn_inh_k);
    arpra_neg(&syn_inh_k, &syn_inh_k);

    // Set constants
    arpra_set_d(&one, 1.0);
    arpra_set_d(&two, 2.0);
    arpra_set_d(&neg_two, -2.0);

    // Initialise report files
    FILE **f_time_c = malloc(sizeof(FILE *));;
    FILE **f_time_r = malloc(sizeof(FILE *));;
    FILE **f_time_n = malloc(sizeof(FILE *));;
    FILE **f_time_s = malloc(sizeof(FILE *));;
    FILE **f_time_d = malloc(sizeof(FILE *));;
    file_init("time", 1, f_time_c, f_time_r, f_time_n, f_time_s, f_time_d);

    FILE **f_nrn1_N_c = malloc(p_nrn1_size * sizeof(FILE *));
    FILE **f_nrn1_N_r = malloc(p_nrn1_size * sizeof(FILE *));
    FILE **f_nrn1_N_n = malloc(p_nrn1_size * sizeof(FILE *));
    FILE **f_nrn1_N_s = malloc(p_nrn1_size * sizeof(FILE *));
    FILE **f_nrn1_N_d = malloc(p_nrn1_size * sizeof(FILE *));
    file_init("nrn1_N", p_nrn1_size, f_nrn1_N_c, f_nrn1_N_r, f_nrn1_N_n, f_nrn1_N_s, f_nrn1_N_d);
    FILE **f_nrn1_V_c = malloc(p_nrn1_size * sizeof(FILE *));
    FILE **f_nrn1_V_r = malloc(p_nrn1_size * sizeof(FILE *));
    FILE **f_nrn1_V_n = malloc(p_nrn1_size * sizeof(FILE *));
    FILE **f_nrn1_V_s = malloc(p_nrn1_size * sizeof(FILE *));
    FILE **f_nrn1_V_d = malloc(p_nrn1_size * sizeof(FILE *));
    file_init("nrn1_V", p_nrn1_size, f_nrn1_V_c, f_nrn1_V_r, f_nrn1_V_n, f_nrn1_V_s, f_nrn1_V_d);

    /* FILE **f_nrn2_N_c = malloc(p_nrn2_size * sizeof(FILE *)); */
    /* FILE **f_nrn2_N_r = malloc(p_nrn2_size * sizeof(FILE *)); */
    /* FILE **f_nrn2_N_n = malloc(p_nrn2_size * sizeof(FILE *)); */
    /* FILE **f_nrn2_N_s = malloc(p_nrn2_size * sizeof(FILE *)); */
    /* FILE **f_nrn2_N_d = malloc(p_nrn2_size * sizeof(FILE *)); */
    /* file_init("nrn2_N", p_nrn2_size, f_nrn2_N_c, f_nrn2_N_r, f_nrn2_N_n, f_nrn2_N_s, f_nrn2_N_d); */
    /* FILE **f_nrn2_V_c = malloc(p_nrn2_size * sizeof(FILE *)); */
    /* FILE **f_nrn2_V_r = malloc(p_nrn2_size * sizeof(FILE *)); */
    /* FILE **f_nrn2_V_n = malloc(p_nrn2_size * sizeof(FILE *)); */
    /* FILE **f_nrn2_V_s = malloc(p_nrn2_size * sizeof(FILE *)); */
    /* FILE **f_nrn2_V_d = malloc(p_nrn2_size * sizeof(FILE *)); */
    /* file_init("nrn2_V", p_nrn2_size, f_nrn2_V_c, f_nrn2_V_r, f_nrn2_V_n, f_nrn2_V_s, f_nrn2_V_d); */

    /* FILE **f_syn_exc_R_c = malloc(p_syn_exc_size * sizeof(FILE *)); */
    /* FILE **f_syn_exc_R_r = malloc(p_syn_exc_size * sizeof(FILE *)); */
    /* FILE **f_syn_exc_R_n = malloc(p_syn_exc_size * sizeof(FILE *)); */
    /* FILE **f_syn_exc_R_s = malloc(p_syn_exc_size * sizeof(FILE *)); */
    /* FILE **f_syn_exc_R_d = malloc(p_syn_exc_size * sizeof(FILE *)); */
    /* file_init("syn_exc_R", p_syn_exc_size, f_syn_exc_R_c, f_syn_exc_R_r, f_syn_exc_R_n, f_syn_exc_R_s, f_syn_exc_R_d); */
    /* FILE **f_syn_exc_S_c = malloc(p_syn_exc_size * sizeof(FILE *)); */
    /* FILE **f_syn_exc_S_r = malloc(p_syn_exc_size * sizeof(FILE *)); */
    /* FILE **f_syn_exc_S_n = malloc(p_syn_exc_size * sizeof(FILE *)); */
    /* FILE **f_syn_exc_S_s = malloc(p_syn_exc_size * sizeof(FILE *)); */
    /* FILE **f_syn_exc_S_d = malloc(p_syn_exc_size * sizeof(FILE *)); */
    /* file_init("syn_exc_S", p_syn_exc_size, f_syn_exc_S_c, f_syn_exc_S_r, f_syn_exc_S_n, f_syn_exc_S_s, f_syn_exc_S_d); */

    /* FILE **f_syn_inh_R_c = malloc(p_syn_inh_size * sizeof(FILE *)); */
    /* FILE **f_syn_inh_R_r = malloc(p_syn_inh_size * sizeof(FILE *)); */
    /* FILE **f_syn_inh_R_n = malloc(p_syn_inh_size * sizeof(FILE *)); */
    /* FILE **f_syn_inh_R_s = malloc(p_syn_inh_size * sizeof(FILE *)); */
    /* FILE **f_syn_inh_R_d = malloc(p_syn_inh_size * sizeof(FILE *)); */
    /* file_init("syn_inh_R", p_syn_inh_size, f_syn_inh_R_c, f_syn_inh_R_r, f_syn_inh_R_n, f_syn_inh_R_s, f_syn_inh_R_d); */
    /* FILE **f_syn_inh_S_c = malloc(p_syn_inh_size * sizeof(FILE *)); */
    /* FILE **f_syn_inh_S_r = malloc(p_syn_inh_size * sizeof(FILE *)); */
    /* FILE **f_syn_inh_S_n = malloc(p_syn_inh_size * sizeof(FILE *)); */
    /* FILE **f_syn_inh_S_s = malloc(p_syn_inh_size * sizeof(FILE *)); */
    /* FILE **f_syn_inh_S_d = malloc(p_syn_inh_size * sizeof(FILE *)); */
    /* file_init("syn_inh_S", p_syn_inh_size, f_syn_inh_S_c, f_syn_inh_S_r, f_syn_inh_S_n, f_syn_inh_S_s, f_syn_inh_S_d); */

    // Set parameter structs
    struct dNdt_params params_nrn1_N = {
        .grp_V = grp_nrn1_V,
        .V3 = &V3,
        .V4 = &V4,
        .phi = &phi,
        .one = &one,
        .two = &two,
        .neg_two = &neg_two,
        .N_ss = &N_ss,
        .temp1 = &temp1,
        .temp2 = &temp2,
    };

    struct dVdt_params params_nrn1_V = {
        .grp_N = grp_nrn1_N,
        .grp_S = grp_syn_exc_S,
        .GSyn = syn_exc_GSyn,
        .VSyn = &syn_exc_VSyn,
        .GL = &GL,
        .VL = &VL,
        .GCa = &GCa,
        .VCa = &VCa,
        .GK = &GK,
        .VK = &VK,
        .V1 = &V1,
        .V2 = &V2,
        .V1 = &V1,
        .V2 = &V2,
        .C = &C,
        .pre_syn_size = p_in1_size,
        .one = &one,
        .neg_two = &neg_two,
        .I = I1,
        .M_ss = &M_ss,
        .temp1 = &temp1,
    };

    struct dNdt_params params_nrn2_N = {
        .grp_V = grp_nrn2_V,
        .V3 = &V3,
        .V4 = &V4,
        .phi = &phi,
        .one = &one,
        .two = &two,
        .neg_two = &neg_two,
        .N_ss = &N_ss,
        .temp1 = &temp1,
        .temp2 = &temp2,
    };

    struct dVdt_params params_nrn2_V = {
        .grp_N = grp_nrn2_N,
        .grp_S = grp_syn_inh_S,
        .GSyn = syn_inh_GSyn,
        .VSyn = &syn_inh_VSyn,
        .GL = &GL,
        .VL = &VL,
        .GCa = &GCa,
        .VCa = &VCa,
        .GK = &GK,
        .VK = &VK,
        .V1 = &V1,
        .V2 = &V2,
        .C = &C,
        .pre_syn_size = p_in2_size,
        .one = &one,
        .neg_two = &neg_two,
        .I = I2,
        .M_ss = &M_ss,
        .temp1 = &temp1,
    };

    struct dRdt_params params_syn_exc_R = {
        .a = &syn_exc_a,
        .b = &syn_exc_b,
        .k = &syn_exc_k,
        .VPre_lo = &in1_V_lo,
        .VPre_hi = &in1_V_hi,
        .threshold = &syn_exc_thr,
        .in = in1,
        .pre_syn_size = p_in1_size,
        .one = &one,
        .temp1 = &temp1,
        .temp2 = &temp2,
    };

    struct dSdt_params params_syn_exc_S = {
        .grp_R = grp_syn_exc_R,
        .a = &syn_exc_a,
        .b = &syn_exc_b,
        .temp1 = &temp1,
        .temp2 = &temp2,
    };

    struct dRdt_params params_syn_inh_R = {
        .a = &syn_inh_a,
        .b = &syn_inh_b,
        .k = &syn_inh_k,
        .VPre_lo = &in2_V_lo,
        .VPre_hi = &in2_V_hi,
        .threshold = &syn_inh_thr,
        .in = in2,
        .pre_syn_size = p_in2_size,
        .one = &one,
        .temp1 = &temp1,
        .temp2 = &temp2,
    };

    struct dSdt_params params_syn_inh_S = {
        .grp_R = grp_syn_inh_R,
        .a = &syn_inh_a,
        .b = &syn_inh_b,
        .temp1 = &temp1,
        .temp2 = &temp2,
    };

    // ODE system
    arpra_uint sys_grps = 8;
    arpra_uint sys_dims[8] = {
        p_nrn1_size, p_nrn1_size, p_nrn2_size, p_nrn2_size,
        p_syn_exc_size, p_syn_exc_size, p_syn_inh_size, p_syn_inh_size
    };
    arpra_ode_f sys_f[8] = {
        dNdt, dVdt, dNdt, dVdt, dRdt, dSdt, dRdt, dSdt
    };
    void *sys_params[8] = {
        &params_nrn1_N, &params_nrn1_V, &params_nrn2_N, &params_nrn2_V,
        &params_syn_exc_R, &params_syn_exc_S, &params_syn_inh_R, &params_syn_inh_S
    };
    arpra_range *sys_x[8] = {
        nrn1_N, nrn1_V, nrn2_N, nrn2_V, syn_exc_R, syn_exc_S, syn_inh_R, syn_inh_S
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

        for (j = 0; j < p_nrn1_size; j++) {
            nrn1_N_reduce_epoch[j] = ode_system.x[grp_nrn1_N][j].nTerms;
            nrn1_V_reduce_epoch[j] = ode_system.x[grp_nrn1_V][j].nTerms;
        }
        for (j = 0; j < p_nrn2_size; j++) {
            nrn2_N_reduce_epoch[j] = ode_system.x[grp_nrn2_N][j].nTerms;
            nrn2_V_reduce_epoch[j] = ode_system.x[grp_nrn2_V][j].nTerms;
        }
        for (j = 0; j < p_syn_exc_size; j++) {
            syn_exc_R_reduce_epoch[j] = ode_system.x[grp_syn_exc_R][j].nTerms;
            syn_exc_S_reduce_epoch[j] = ode_system.x[grp_syn_exc_S][j].nTerms;
        }
        for (j = 0; j < p_syn_inh_size; j++) {
            syn_inh_R_reduce_epoch[j] = ode_system.x[grp_syn_inh_R][j].nTerms;
            syn_inh_S_reduce_epoch[j] = ode_system.x[grp_syn_inh_S][j].nTerms;
        }

        // Event(s) occur if urandom >= e^-rate
        for (j = 0; j < p_in1_size; j++) {
            mpfr_urandom(&rand_uf, rng_uf, MPFR_RNDN);
            in1[j] = mpfr_greaterequal_p(&rand_uf, &in1_p0);
            fprintf(stderr, "%s", (in1[j] ? "\x1B[31m\xE2\x96\xA3\x1B[0m" : "\xE2\x96\xA3"));
        }
        fprintf(stderr, "  ");
        for (j = 0; j < p_in2_size; j++) {
            mpfr_urandom(&rand_uf, rng_uf, MPFR_RNDN);
            in2[j] = mpfr_greaterequal_p(&rand_uf, &in2_p0);
            fprintf(stderr, "%s", (in2[j] ? "\x1B[31m\xE2\x96\xA3\x1B[0m" : "\xE2\x96\xA3"));
        }
        fprintf(stderr, "\n");

        // Step system
        arpra_ode_stepper_step(&ode_stepper, &h);

        arpra_uint reduce_n;
        for (j = 0; j < p_nrn1_size; j++) {
            reduce_n = ode_system.x[grp_nrn1_N][j].nTerms - nrn1_N_reduce_epoch[j];
            arpra_reduce_last_n(&(ode_system.x[grp_nrn1_N][j]), reduce_n);
            reduce_n = ode_system.x[grp_nrn1_V][j].nTerms - nrn1_V_reduce_epoch[j];
            arpra_reduce_last_n(&(ode_system.x[grp_nrn1_V][j]), reduce_n);
        }
        for (j = 0; j < p_nrn2_size; j++) {
            reduce_n = ode_system.x[grp_nrn2_N][j].nTerms - nrn2_N_reduce_epoch[j];
            arpra_reduce_last_n(&(ode_system.x[grp_nrn2_N][j]), reduce_n);
            reduce_n = ode_system.x[grp_nrn2_V][j].nTerms - nrn2_V_reduce_epoch[j];
            arpra_reduce_last_n(&(ode_system.x[grp_nrn2_V][j]), reduce_n);
        }
        for (j = 0; j < p_syn_exc_size; j++) {
            reduce_n = ode_system.x[grp_syn_exc_R][j].nTerms - syn_exc_R_reduce_epoch[j];
            arpra_reduce_last_n(&(ode_system.x[grp_syn_exc_R][j]), reduce_n);
            reduce_n = ode_system.x[grp_syn_exc_S][j].nTerms - syn_exc_S_reduce_epoch[j];
            arpra_reduce_last_n(&(ode_system.x[grp_syn_exc_S][j]), reduce_n);
        }
        for (j = 0; j < p_syn_inh_size; j++) {
            reduce_n = ode_system.x[grp_syn_inh_R][j].nTerms - syn_inh_R_reduce_epoch[j];
            arpra_reduce_last_n(&(ode_system.x[grp_syn_inh_R][j]), reduce_n);
            reduce_n = ode_system.x[grp_syn_inh_S][j].nTerms - syn_inh_S_reduce_epoch[j];
            arpra_reduce_last_n(&(ode_system.x[grp_syn_inh_S][j]), reduce_n);
        }

        if (i % p_reduce_step == 0) {
            for (j = 0; j < p_nrn1_size; j++) {
                arpra_reduce_small(&(ode_system.x[grp_nrn1_N][j]), p_reduce_ratio);
                arpra_reduce_small(&(ode_system.x[grp_nrn1_V][j]), p_reduce_ratio);
            }
            for (j = 0; j < p_nrn2_size; j++) {
                arpra_reduce_small(&(ode_system.x[grp_nrn2_N][j]), p_reduce_ratio);
                arpra_reduce_small(&(ode_system.x[grp_nrn2_V][j]), p_reduce_ratio);
            }
            for (j = 0; j < p_syn_exc_size; j++) {
                arpra_reduce_small(&(ode_system.x[grp_syn_exc_R][j]), p_reduce_ratio);
                arpra_reduce_small(&(ode_system.x[grp_syn_exc_S][j]), p_reduce_ratio);
            }
            for (j = 0; j < p_syn_inh_size; j++) {
                arpra_reduce_small(&(ode_system.x[grp_syn_inh_R][j]), p_reduce_ratio);
                arpra_reduce_small(&(ode_system.x[grp_syn_inh_S][j]), p_reduce_ratio);
            }
        }

        file_write(&sys_t, 1, f_time_c, f_time_r, f_time_n, f_time_s, f_time_d);

        file_write(nrn1_N, p_nrn1_size, f_nrn1_N_c, f_nrn1_N_r, f_nrn1_N_n, f_nrn1_N_s, f_nrn1_N_d);
        file_write(nrn1_V, p_nrn1_size, f_nrn1_V_c, f_nrn1_V_r, f_nrn1_V_n, f_nrn1_V_s, f_nrn1_V_d);

        /* file_write(nrn2_N, p_nrn2_size, f_nrn2_N_c, f_nrn2_N_r, f_nrn2_N_n, f_nrn2_N_s, f_nrn2_N_d); */
        /* file_write(nrn2_V, p_nrn2_size, f_nrn2_V_c, f_nrn2_V_r, f_nrn2_V_n, f_nrn2_V_s, f_nrn2_V_d); */

        /* file_write(syn_exc_R, p_syn_exc_size, f_syn_exc_R_c, f_syn_exc_R_r, f_syn_exc_R_n, f_syn_exc_R_s, f_syn_exc_R_d); */
        /* file_write(syn_exc_S, p_syn_exc_size, f_syn_exc_S_c, f_syn_exc_S_r, f_syn_exc_S_n, f_syn_exc_S_s, f_syn_exc_S_d); */

        /* file_write(syn_inh_R, p_syn_inh_size, f_syn_inh_R_c, f_syn_inh_R_r, f_syn_inh_R_n, f_syn_inh_R_s, f_syn_inh_R_d); */
        /* file_write(syn_inh_S, p_syn_inh_size, f_syn_inh_S_c, f_syn_inh_S_r, f_syn_inh_S_n, f_syn_inh_S_s, f_syn_inh_S_d); */
    }

    run_time = clock() - run_time;
    printf("Finished in %f seconds.\n", ((float) run_time) / CLOCKS_PER_SEC);

    // End simulation loop
    // ===================


    // Clear system state
    arpra_clear(&h);
    arpra_clear(&sys_t);
    for (i = 0; i < p_nrn1_size; i++) {
        arpra_clear(&(nrn1_N[i]));
        arpra_clear(&(nrn1_V[i]));
    }
    for (i = 0; i < p_nrn2_size; i++) {
        arpra_clear(&(nrn2_N[i]));
        arpra_clear(&(nrn2_V[i]));
    }
    for (i = 0; i < p_syn_exc_size; i++) {
        arpra_clear(&(syn_exc_R[i]));
        arpra_clear(&(syn_exc_S[i]));
    }
    for (i = 0; i < p_syn_inh_size; i++) {
        arpra_clear(&(syn_inh_R[i]));
        arpra_clear(&(syn_inh_S[i]));
    }

    // Clear Poisson input parameters (group 1)
    mpfr_clear(&in1_p0);
    arpra_clear(&in1_V_lo);
    arpra_clear(&in1_V_hi);

    // Clear Poisson input parameters (group 2)
    mpfr_clear(&in2_p0);
    arpra_clear(&in2_V_lo);
    arpra_clear(&in2_V_hi);

    // Clear neuron parameters
    arpra_clear(&GL);
    arpra_clear(&GCa);
    arpra_clear(&GK);
    arpra_clear(&VL);
    arpra_clear(&VCa);
    arpra_clear(&VK);
    arpra_clear(&V1);
    arpra_clear(&V2);
    arpra_clear(&V3);
    arpra_clear(&V4);
    arpra_clear(&phi);
    arpra_clear(&C);

    // Clear excitatory synapse parameters
    for (i = 0; i < p_syn_exc_size; i++) {
        arpra_clear(&(syn_exc_GSyn[i]));
    }
    arpra_clear(&syn_exc_VSyn);
    arpra_clear(&syn_exc_thr);
    arpra_clear(&syn_exc_a);
    arpra_clear(&syn_exc_b);
    arpra_clear(&syn_exc_k);

    // Clear inhibitory synapse parameters
    for (i = 0; i < p_syn_inh_size; i++) {
        arpra_clear(&(syn_inh_GSyn[i]));
    }
    arpra_clear(&syn_inh_VSyn);
    arpra_clear(&syn_inh_thr);
    arpra_clear(&syn_inh_a);
    arpra_clear(&syn_inh_b);
    arpra_clear(&syn_inh_k);

    // Clear constants
    arpra_clear(&one);
    arpra_clear(&two);
    arpra_clear(&neg_two);

    // Clear scratch space
    mpfr_clear(&rand_uf);
    mpfr_clear(&rand_nf);
    arpra_clear(&temp1);
    arpra_clear(&temp2);
    arpra_clear(&M_ss);
    arpra_clear(&N_ss);
    for (i = 0; i < p_in1_size; i++) {
        arpra_clear(&(I1[i]));
    }
    for (i = 0; i < p_in2_size; i++) {
        arpra_clear(&(I2[i]));
    }

    // Free system state
    free(nrn1_N);
    free(nrn1_V);
    free(nrn2_N);
    free(nrn2_V);
    free(syn_exc_R);
    free(syn_exc_S);
    free(syn_inh_R);
    free(syn_inh_S);

    // Free other arrays
    free(syn_exc_GSyn);
    free(syn_inh_GSyn);
    free(I1);
    free(I2);
    free(in1);
    free(in2);
    free(nrn1_N_reduce_epoch);
    free(nrn1_V_reduce_epoch);
    free(nrn2_N_reduce_epoch);
    free(nrn2_V_reduce_epoch);
    free(syn_exc_R_reduce_epoch);
    free(syn_exc_S_reduce_epoch);
    free(syn_inh_R_reduce_epoch);
    free(syn_inh_S_reduce_epoch);

    // Clear report files
    file_clear(1, f_time_c, f_time_r, f_time_n, f_time_s, f_time_d);
    free(f_time_c);
    free(f_time_r);
    free(f_time_n);
    free(f_time_s);
    free(f_time_d);

    file_clear(p_nrn1_size, f_nrn1_N_c, f_nrn1_N_r, f_nrn1_N_n, f_nrn1_N_s, f_nrn1_N_d);
    free(f_nrn1_N_c);
    free(f_nrn1_N_r);
    free(f_nrn1_N_n);
    free(f_nrn1_N_s);
    free(f_nrn1_N_d);
    file_clear(p_nrn1_size, f_nrn1_V_c, f_nrn1_V_r, f_nrn1_V_n, f_nrn1_V_s, f_nrn1_V_d);
    free(f_nrn1_V_c);
    free(f_nrn1_V_r);
    free(f_nrn1_V_n);
    free(f_nrn1_V_s);
    free(f_nrn1_V_d);

    /* file_clear(p_nrn2_size, f_nrn2_N_c, f_nrn2_N_r, f_nrn2_N_n, f_nrn2_N_s, f_nrn2_N_d); */
    /* free(f_nrn2_N_c); */
    /* free(f_nrn2_N_r); */
    /* free(f_nrn2_N_n); */
    /* free(f_nrn2_N_s); */
    /* free(f_nrn2_N_d); */
    /* file_clear(p_nrn2_size, f_nrn2_V_c, f_nrn2_V_r, f_nrn2_V_n, f_nrn2_V_s, f_nrn2_V_d); */
    /* free(f_nrn2_V_c); */
    /* free(f_nrn2_V_r); */
    /* free(f_nrn2_V_n); */
    /* free(f_nrn2_V_s); */
    /* free(f_nrn2_V_d); */

    /* file_clear(p_syn_exc_size, f_syn_exc_R_c, f_syn_exc_R_r, f_syn_exc_R_n, f_syn_exc_R_s, f_syn_exc_R_d); */
    /* free(f_syn_exc_R_c); */
    /* free(f_syn_exc_R_r); */
    /* free(f_syn_exc_R_n); */
    /* free(f_syn_exc_R_s); */
    /* free(f_syn_exc_R_d); */
    /* file_clear(p_syn_exc_size, f_syn_exc_S_c, f_syn_exc_S_r, f_syn_exc_S_n, f_syn_exc_S_s, f_syn_exc_S_d); */
    /* free(f_syn_exc_S_c); */
    /* free(f_syn_exc_S_r); */
    /* free(f_syn_exc_S_n); */
    /* free(f_syn_exc_S_s); */
    /* free(f_syn_exc_S_d); */

    /* file_clear(p_syn_inh_size, f_syn_inh_R_c, f_syn_inh_R_r, f_syn_inh_R_n, f_syn_inh_R_s, f_syn_inh_R_d); */
    /* free(f_syn_inh_R_c); */
    /* free(f_syn_inh_R_r); */
    /* free(f_syn_inh_R_n); */
    /* free(f_syn_inh_R_s); */
    /* free(f_syn_inh_R_d); */
    /* file_clear(p_syn_inh_size, f_syn_inh_S_c, f_syn_inh_S_r, f_syn_inh_S_n, f_syn_inh_S_s, f_syn_inh_S_d); */
    /* free(f_syn_inh_S_c); */
    /* free(f_syn_inh_S_r); */
    /* free(f_syn_inh_S_n); */
    /* free(f_syn_inh_S_s); */
    /* free(f_syn_inh_S_d); */

    arpra_ode_stepper_clear(&ode_stepper);
    arpra_clear_buffers();
    gmp_randclear(rng_uf);
    gmp_randclear(rng_nf);
    mpfr_free_cache();

    return 0;
}
