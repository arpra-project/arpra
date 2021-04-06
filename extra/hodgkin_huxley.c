/*
 * hodgkin_huxley.c -- Test Hodgkin-Huxley model.
 *
 * Copyright 2021 James Paul Turner.
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
 * M      Probability of Na+ channel activation
 * H      Probability of Na+ channel not blocking
 * N      Probability of K+ channel activation
 * V      Membrane potential (mV)
 *
 * Neuron parameters
 * -----------------
 * GL     Maximum leak conductance (mS/cm^2)
 * VL     Equilibrium potential of leak conductance (mV)
 * GNa    Maximum Na+ conductance (mS/cm^2)
 * VNa    Equilibrium potential of Na+ conductance (mV)
 * GK     Maximum K+ conductance (mS/cm^2)
 * VK     Equilibrium potential of K+ conductance (mV)
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
//#define p_h 0.5
//#define p_h 0.1
#define p_h 0.05
#define p_t0 0.0
#define p_prec 53
#define p_prec_internal 256
#define p_sim_steps 10000
#define p_report_step 20
#define p_reduce_step 100
#define p_reduce_rel 0.1

// RNG parameters
// Seeds are random if not #defined
#define p_rand_prec 53
//#define p_rng_uf_seed 707135875931ul
//#define p_rng_nf_seed 503108552855ul

// Poisson input parameters
//#define p_in_size 50
#define p_in_size in_size_arg
int in_size_arg;
//#define p_in_freq 10.0
#define p_in_freq in_freq_arg
int in_freq_arg;
#define p_in_V_lo -60.0
#define p_in_V_hi 20.0

// Neuron parameters
#define p_nrn_size 1
#define p_nrn_M0 0.0
#define p_nrn_H0 0.0
#define p_nrn_N0 0.0
#define p_nrn_V0 -60.0
#define p_nrn_GL 0.02672
#define p_nrn_GNa 7.15
#define p_nrn_GK 1.43
#define p_nrn_VL -63.563
#define p_nrn_VNa 50.0
#define p_nrn_VK -95.0
#define p_nrn_C 0.143

// Synapse parameters (excitatory)
#define p_syn_size p_in_size * p_nrn_size
#define p_syn_R0 0.0
#define p_syn_S0 0.0
#define p_syn_GSyn_std 0.1
#define p_syn_GSyn_mean 0.0
#define p_syn_VSyn 0.0
#define p_syn_thr -50.0
#define p_syn_a 0.25 // in [1/10, 1/2]
#define p_syn_b 0.15 // in [1/20, 1/4]
#define p_syn_k 1.0E6

/* // Synapse parameters (inhibitory) */
/* #define p_syn_size p_in_size * p_nrn_size */
/* #define p_syn_R0 0.0 */
/* #define p_syn_S0 0.0 */
/* #define p_syn_GSyn_std 1.0 */
/* #define p_syn_GSyn_mean 0.0 */
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

struct dMdt_params
{
    arpra_uint grp_V;
    arpra_range *pos_1;
    arpra_range *pos_4;
    arpra_range *pos_5;
    arpra_range *pos_25;
    arpra_range *neg_52;
    arpra_range *pos_028;
    arpra_range *pos_032;
    arpra_range *temp1;
    arpra_range *temp2;
    arpra_range *_a;
    arpra_range *_b;
};

void dMdt (arpra_range *y, const void *params,
           const arpra_range *t, const arpra_range **x,
           const arpra_uint x_grp, const arpra_uint x_dim)
{
    const struct dMdt_params *p = (struct dMdt_params *) params;
    const arpra_range *M = &(x[x_grp][x_dim]);
    const arpra_range *V = &(x[p->grp_V][x_dim]);
    const arpra_range *pos_1 = p->pos_1;
    const arpra_range *pos_4 = p->pos_4;
    const arpra_range *pos_5 = p->pos_5;
    const arpra_range *pos_25 = p->pos_25;
    const arpra_range *neg_52 = p->neg_52;
    const arpra_range *pos_028 = p->pos_028;
    const arpra_range *pos_032 = p->pos_032;
    arpra_range *temp1 = p->temp1;
    arpra_range *temp2 = p->temp2;
    arpra_range *_a = p->_a;
    arpra_range *_b = p->_b;

    // Compute M_a
    // M_a = 0.32 * (-52.0 - V) / (exp((-52.0 - V) / 4.0) - 1.0)
    arpra_sub(temp1, neg_52, V);
    arpra_div(_a, temp1, pos_4);
    arpra_exp(_a, _a);
    arpra_sub(_a, _a, pos_1);
    arpra_div(_a, temp1, _a);
    arpra_mul(_a, pos_032, _a);

    // Compute M_b
    // M_b = 0.28 * (V + 25.0) / (exp((V + 25.0) / 5.0) - 1.0)
    arpra_add(temp2, V, pos_25);
    arpra_div(_b, temp2, pos_5);
    arpra_exp(_b, _b);
    arpra_sub(_b, _b, pos_1);
    arpra_div(_b, temp2, _b);
    arpra_mul(_b, pos_028, _b);

    // delta of M
    // dM/dt = (M_a * (1.0 - M) - M_b * M)
    arpra_sub(temp1, pos_1, M);
    arpra_mul(temp1, _a, temp1);
    arpra_mul(temp2, _b, M);
    arpra_sub(y, temp1, temp2);
}

struct dHdt_params
{
    arpra_uint grp_V;
    arpra_range *pos_1;
    arpra_range *pos_4;
    arpra_range *pos_5;
    arpra_range *pos_18;
    arpra_range *neg_25;
    arpra_range *neg_48;
    arpra_range *pos_0128;
    arpra_range *temp1;
    arpra_range *temp2;
    arpra_range *_a;
    arpra_range *_b;
};

void dHdt (arpra_range *y, const void *params,
           const arpra_range *t, const arpra_range **x,
           const arpra_uint x_grp, const arpra_uint x_dim)
{
    const struct dHdt_params *p = (struct dHdt_params *) params;
    const arpra_range *H = &(x[x_grp][x_dim]);
    const arpra_range *V = &(x[p->grp_V][x_dim]);
    const arpra_range *pos_1 = p->pos_1;
    const arpra_range *pos_4 = p->pos_4;
    const arpra_range *pos_5 = p->pos_5;
    const arpra_range *pos_18 = p->pos_18;
    const arpra_range *neg_25 = p->neg_25;
    const arpra_range *neg_48 = p->neg_48;
    const arpra_range *pos_0128 = p->pos_0128;
    arpra_range *temp1 = p->temp1;
    arpra_range *temp2 = p->temp2;
    arpra_range *_a = p->_a;
    arpra_range *_b = p->_b;

    // Compute H_a
    // H_a = 0.128 * exp((-48.0 - V) / 18.0)
    arpra_sub(_a, neg_48, V);
    arpra_div(_a, _a, pos_18);
    arpra_exp(_a, _a);
    arpra_mul(_a, pos_0128, _a);

    // Compute H_b
    // H_b = 4.0 / (exp((-25.0 - V) / 5.0) + 1.0)
    arpra_sub(_b, neg_25, V);
    arpra_div(_b, _b, pos_5);
    arpra_exp(_b, _b);
    arpra_add(_b, _b, pos_1);
    arpra_div(_b, pos_4, _b);

    // delta of H
    // dH/dt = (H_a * (1.0 - H) - H_b * H)
    arpra_sub(temp1, pos_1, H);
    arpra_mul(temp1, _a, temp1);
    arpra_mul(temp2, _b, H);
    arpra_sub(y, temp1, temp2);
}

struct dNdt_params
{
    arpra_uint grp_V;
    arpra_range *pos_1;
    arpra_range *pos_5;
    arpra_range *pos_40;
    arpra_range *neg_50;
    arpra_range *neg_55;
    arpra_range *pos_0032;
    arpra_range *pos_05;
    arpra_range *temp1;
    arpra_range *temp2;
    arpra_range *_a;
    arpra_range *_b;
};

void dNdt (arpra_range *y, const void *params,
           const arpra_range *t, const arpra_range **x,
           const arpra_uint x_grp, const arpra_uint x_dim)
{
    const struct dNdt_params *p = (struct dNdt_params *) params;
    const arpra_range *N = &(x[x_grp][x_dim]);
    const arpra_range *V = &(x[p->grp_V][x_dim]);
    const arpra_range *pos_1 = p->pos_1;
    const arpra_range *pos_5 = p->pos_5;
    const arpra_range *pos_40 = p->pos_40;
    const arpra_range *neg_50 = p->neg_50;
    const arpra_range *neg_55 = p->neg_55;
    const arpra_range *pos_0032 = p->pos_0032;
    const arpra_range *pos_05 = p->pos_05;
    arpra_range *temp1 = p->temp1;
    arpra_range *temp2 = p->temp2;
    arpra_range *_a = p->_a;
    arpra_range *_b = p->_b;

    // Compute N_a
    // N_a = 0.032 * (-50.0 - V) / (exp((-50.0 - V) / 5.0) - 1.0)
    arpra_sub(temp1, neg_50, V);
    arpra_div(_a, temp1, pos_5);
    arpra_exp(_a, _a);
    arpra_sub(_a, _a, pos_1);
    arpra_div(_a, temp1, _a);
    arpra_mul(_a, pos_0032, _a);

    // Compute N_b
    // N_b = 0.5 * exp((-55.0 - V) / 40.0)
    arpra_sub(_b, neg_55, V);
    arpra_div(_b, _b, pos_40);
    arpra_exp(_b, _b);
    arpra_mul(_b, pos_05, _b);

    // delta of N
    // dH/dt = (N_a * (1.0 - N) - N_b * N)
    arpra_sub(temp1, pos_1, N);
    arpra_mul(temp1, _a, temp1);
    arpra_mul(temp2, _b, N);
    arpra_sub(y, temp1, temp2);
}

struct dVdt_params
{
    arpra_uint grp_M;
    arpra_uint grp_H;
    arpra_uint grp_N;
    arpra_uint grp_S;
    arpra_range *GSyn;
    arpra_range *VSyn;
    arpra_range *GL;
    arpra_range *VL;
    arpra_range *GNa;
    arpra_range *VNa;
    arpra_range *GK;
    arpra_range *VK;
    arpra_range *C;
    arpra_uint pre_syn_size;
    arpra_range *I;
    arpra_range *temp1;
};

void dVdt (arpra_range *y, const void *params,
           const arpra_range *t, const arpra_range **x,
           const arpra_uint x_grp, const arpra_uint x_dim)
{
    const struct dVdt_params *p = (struct dVdt_params *) params;
    const arpra_range *V = &(x[x_grp][x_dim]);
    const arpra_range *M = &(x[p->grp_M][x_dim]);
    const arpra_range *H = &(x[p->grp_H][x_dim]);
    const arpra_range *N = &(x[p->grp_N][x_dim]);
    const arpra_range *S = &(x[p->grp_S][x_dim * p->pre_syn_size]);
    const arpra_range *GSyn = &(p->GSyn[x_dim * p->pre_syn_size]);
    const arpra_range *VSyn = p->VSyn;
    const arpra_range *GL = p->GL;
    const arpra_range *VL = p->VL;
    const arpra_range *GNa = p->GNa;
    const arpra_range *VNa = p->VNa;
    const arpra_range *GK = p->GK;
    const arpra_range *VK = p->VK;
    const arpra_range *C = p->C;
    arpra_range *I = p->I;
    arpra_range *temp1 = p->temp1;
    arpra_uint i;

    // Synapse current
    arpra_sub(temp1, VSyn, V);
    for (i = 0; i < p->pre_syn_size; i++) {
        arpra_mul(&(I[i]), temp1, &(GSyn[i]));
        arpra_mul(&(I[i]), &(I[i]), &(S[i]));
    }
    if (p->pre_syn_size > 0) {
        arpra_sum_recursive(y, I, p->pre_syn_size);
        //arpra_sum(y, I, p->pre_syn_size);
    }
    else {
        arpra_set_zero(y);
    }

    // Leak current
    arpra_sub(temp1, V, VL);
    arpra_mul(temp1, temp1, GL);
    arpra_sub(y, y, temp1);

    // Na+ current
    arpra_sub(temp1, V, VNa);
    arpra_mul(temp1, temp1, GNa);
    arpra_mul(temp1, temp1, M);
    arpra_mul(temp1, temp1, M);
    arpra_mul(temp1, temp1, M);
    arpra_mul(temp1, temp1, H);
    arpra_sub(y, y, temp1);

    // K+ current
    arpra_sub(temp1, V, VK);
    arpra_mul(temp1, temp1, GK);
    arpra_mul(temp1, temp1, N);
    arpra_mul(temp1, temp1, N);
    arpra_mul(temp1, temp1, N);
    arpra_mul(temp1, temp1, N);
    arpra_sub(y, y, temp1);

    // delta of membrane potential
    // dV/dt = (I + GL (VL - V) + GNa M^3 H (VNa - V) + GK N^4 (VK - V)) / C
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
    arpra_range *pos_1;
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
    const arpra_range *pos_1 = p->pos_1;
    arpra_range *temp1 = p->temp1;
    arpra_range *temp2 = p->temp2;

    // Sigmoid of threshold difference
    arpra_sub(temp1, VPre, threshold);
    arpra_mul(temp1, temp1, k);
    arpra_exp(temp1, temp1);
    arpra_add(temp1, temp1, pos_1);
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
        grp_nrn_M, grp_nrn_H,
        grp_nrn_N, grp_nrn_V,
        grp_syn_R, grp_syn_S,
    };

    if (argc != 4) {
        printf("Usage: hodgkin_huxley <AA method> <input size> <input frequency>\n");
        exit(0);
    }

    int method = atoi(argv[1]);
    in_size_arg = atoi(argv[2]);
    in_freq_arg = atoi(argv[3]);
    printf("in size: %u\n", in_size_arg);
    printf("in freq: %u\n", in_freq_arg);

    if (method == 0) {
        printf("AA method\n");
        arpra_set_range_method(ARPRA_AA);
    }
    else if (method == 1) {
        printf("mixed IA/AA method\n");
        arpra_set_range_method(ARPRA_MIXED_IAAA);
    }
    else if (method == 2) {
        printf("mixed trimmed IA/AA method\n");
        arpra_set_range_method(ARPRA_MIXED_TRIMMED_IAAA);
    }

    arpra_set_internal_precision(p_prec_internal);

    // Initialise arpra_reduce_small_rel threshold.
    mpfr_t reduce_rel;
    mpfr_init2(reduce_rel, 53);
    mpfr_set_d(reduce_rel, p_reduce_rel, MPFR_RNDN);

    // Allocate system state
    arpra_range *nrn_M = malloc(p_nrn_size * sizeof(arpra_range));
    arpra_range *nrn_H = malloc(p_nrn_size * sizeof(arpra_range));
    arpra_range *nrn_N = malloc(p_nrn_size * sizeof(arpra_range));
    arpra_range *nrn_V = malloc(p_nrn_size * sizeof(arpra_range));
    arpra_range *syn_R = malloc(p_syn_size * sizeof(arpra_range));
    arpra_range *syn_S = malloc(p_syn_size * sizeof(arpra_range));

    // Allocate other arrays
    arpra_range *syn_GSyn = malloc(p_syn_size * sizeof(arpra_range));
    arpra_range *I = malloc(p_in_size * sizeof(arpra_range));
    int *in = malloc(p_in_size * sizeof(int));
    arpra_uint *nrn_M_reduce_epoch = malloc(p_nrn_size * sizeof(arpra_uint));
    arpra_uint *nrn_H_reduce_epoch = malloc(p_nrn_size * sizeof(arpra_uint));
    arpra_uint *nrn_N_reduce_epoch = malloc(p_nrn_size * sizeof(arpra_uint));
    arpra_uint *nrn_V_reduce_epoch = malloc(p_nrn_size * sizeof(arpra_uint));
    arpra_uint *syn_R_reduce_epoch = malloc(p_syn_size * sizeof(arpra_uint));
    arpra_uint *syn_S_reduce_epoch = malloc(p_syn_size * sizeof(arpra_uint));

    mpfr_t in_p0, rand_uf, rand_nf;
    arpra_range nrn_GL, nrn_VL, nrn_GNa, nrn_VNa, nrn_GK, nrn_VK,
        nrn_C, syn_VSyn, syn_thr, syn_a, syn_b, syn_k,
        temp1, temp2, _a, _b, in_V_lo, in_V_hi;
    arpra_range pos_1, pos_4, pos_5, pos_25, neg_52, pos_028, pos_032,
	pos_0128, neg_48, pos_18, neg_25,
	pos_05, pos_0032, neg_50, neg_55, pos_40;

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
        arpra_init2(&(nrn_M[i]), p_prec);
        arpra_init2(&(nrn_H[i]), p_prec);
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
    arpra_init2(&nrn_GNa, p_prec);
    arpra_init2(&nrn_GK, p_prec);
    arpra_init2(&nrn_VL, p_prec);
    arpra_init2(&nrn_VNa, p_prec);
    arpra_init2(&nrn_VK, p_prec);
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
    arpra_init2(&pos_1, p_prec);
    arpra_init2(&pos_4, p_prec);
    arpra_init2(&pos_5, p_prec);
    arpra_init2(&pos_18, p_prec);
    arpra_init2(&pos_25, p_prec);
    arpra_init2(&neg_25, p_prec);
    arpra_init2(&neg_48, p_prec);
    arpra_init2(&pos_40, p_prec);
    arpra_init2(&neg_50, p_prec);
    arpra_init2(&neg_52, p_prec);
    arpra_init2(&neg_55, p_prec);
    arpra_init2(&pos_028, p_prec);
    arpra_init2(&pos_032, p_prec);
    arpra_init2(&pos_0128, p_prec);
    arpra_init2(&pos_0032, p_prec);
    arpra_init2(&pos_05, p_prec);

    // Initialise scratch space
    mpfr_init2(rand_uf, p_rand_prec);
    mpfr_init2(rand_nf, p_rand_prec);
    arpra_init2(&temp1, p_prec);
    arpra_init2(&temp2, p_prec);
    arpra_init2(&_a, p_prec);
    arpra_init2(&_b, p_prec);
    for (i = 0; i < p_in_size; i++) {
        arpra_init2(&(I[i]), p_prec);
    }

    // Set system state
    arpra_set_d(&h, p_h);
    arpra_set_d(&sys_t, p_t0);
    for (i = 0; i < p_nrn_size; i++) {
        arpra_set_d(&(nrn_M[i]), p_nrn_M0);
        arpra_set_d(&(nrn_H[i]), p_nrn_H0);
        arpra_set_d(&(nrn_N[i]), p_nrn_N0);
        arpra_set_d(&(nrn_V[i]), p_nrn_V0);
    }
    for (i = 0; i < p_syn_size; i++) {
        arpra_set_d(&(syn_R[i]), p_syn_R0);
        arpra_set_d(&(syn_S[i]), p_syn_S0);
    }

    // Set Poisson input parameters
    mpfr_set_d(in_p0, -((double) p_in_freq / (double) 1000) * p_h, MPFR_RNDN);
    mpfr_exp(in_p0, in_p0, MPFR_RNDN);
    arpra_set_d(&in_V_lo, p_in_V_lo);
    arpra_set_d(&in_V_hi, p_in_V_hi);

    // Set neuron parameters
    arpra_set_d(&nrn_GL, p_nrn_GL);
    arpra_set_d(&nrn_GNa, p_nrn_GNa);
    arpra_set_d(&nrn_GK, p_nrn_GK);
    arpra_set_d(&nrn_VL, p_nrn_VL);
    arpra_set_d(&nrn_VNa, p_nrn_VNa);
    arpra_set_d(&nrn_VK, p_nrn_VK);
    arpra_set_d(&nrn_C, p_nrn_C);

    // Set synapse parameters
    for (i = 0; i < p_syn_size; i++) {
        mpfr_nrandom(rand_nf, rng_nf, MPFR_RNDN);
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
    arpra_set_d(&pos_1, 1.0);
    arpra_set_d(&pos_4, 4.0);
    arpra_set_d(&pos_5, 5.0);
    arpra_set_d(&pos_18, 18.0);
    arpra_set_d(&pos_25, 25.0);
    arpra_set_d(&neg_25, -25.0);
    arpra_set_d(&pos_40, 40.0);
    arpra_set_d(&neg_48, -48.0);
    arpra_set_d(&neg_50, -50.0);
    arpra_set_d(&neg_52, -52.0);
    arpra_set_d(&neg_55, -55.0);
    arpra_set_d(&pos_028, 0.28);
    arpra_set_d(&pos_032, 0.32);
    arpra_set_d(&pos_0128, 0.128);
    arpra_set_d(&pos_0032, 0.032);
    arpra_set_d(&pos_05, 0.5);

    // Initialise report files
    FILE **f_time_c = malloc(sizeof(FILE *));
    FILE **f_time_r = malloc(sizeof(FILE *));
    FILE **f_time_n = malloc(sizeof(FILE *));
    FILE **f_time_s = malloc(sizeof(FILE *));
    FILE **f_time_d = malloc(sizeof(FILE *));
    file_init("time", 1, f_time_c, f_time_r, f_time_n, f_time_s, f_time_d);

    FILE **f_nrn_M_c = malloc(p_nrn_size * sizeof(FILE *));
    FILE **f_nrn_M_r = malloc(p_nrn_size * sizeof(FILE *));
    FILE **f_nrn_M_n = malloc(p_nrn_size * sizeof(FILE *));
    FILE **f_nrn_M_s = malloc(p_nrn_size * sizeof(FILE *));
    FILE **f_nrn_M_d = malloc(p_nrn_size * sizeof(FILE *));
    file_init("nrn_M", p_nrn_size, f_nrn_M_c, f_nrn_M_r, f_nrn_M_n, f_nrn_M_s, f_nrn_M_d);
    FILE **f_nrn_H_c = malloc(p_nrn_size * sizeof(FILE *));
    FILE **f_nrn_H_r = malloc(p_nrn_size * sizeof(FILE *));
    FILE **f_nrn_H_n = malloc(p_nrn_size * sizeof(FILE *));
    FILE **f_nrn_H_s = malloc(p_nrn_size * sizeof(FILE *));
    FILE **f_nrn_H_d = malloc(p_nrn_size * sizeof(FILE *));
    file_init("nrn_H", p_nrn_size, f_nrn_H_c, f_nrn_H_r, f_nrn_H_n, f_nrn_H_s, f_nrn_H_d);
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
    struct dMdt_params params_nrn_M = {
        .grp_V = grp_nrn_V,
        .pos_1 = &pos_1,
        .pos_4 = &pos_4,
        .pos_5 = &pos_5,
        .pos_25 = &pos_25,
        .neg_52 = &neg_52,
        .pos_028 = &pos_028,
        .pos_032 = &pos_032,
        .temp1 = &temp1,
        .temp2 = &temp2,
        ._a = &_a,
        ._b = &_b,
    };

    struct dHdt_params params_nrn_H = {
        .grp_V = grp_nrn_V,
        .pos_1 = &pos_1,
        .pos_4 = &pos_4,
        .pos_5 = &pos_5,
        .pos_18 = &pos_18,
        .neg_25 = &neg_25,
        .neg_48 = &neg_48,
        .pos_0128 = &pos_0128,
        .temp1 = &temp1,
        .temp2 = &temp2,
        ._a = &_a,
        ._b = &_b,
    };

    struct dNdt_params params_nrn_N = {
        .grp_V = grp_nrn_V,
	.pos_1 = &pos_1,
	.pos_5 = &pos_5,
	.pos_40 = &pos_40,
	.neg_50 = &neg_50,
	.neg_55 = &neg_55,
	.pos_0032 = &pos_0032,
	.pos_05 = &pos_05,
        .temp1 = &temp1,
        .temp2 = &temp2,
        ._a = &_a,
        ._b = &_b,
    };

    struct dVdt_params params_nrn_V = {
        .grp_M = grp_nrn_M,
        .grp_H = grp_nrn_H,
        .grp_N = grp_nrn_N,
        .grp_S = grp_syn_S,
        .GSyn = syn_GSyn,
        .VSyn = &syn_VSyn,
        .GL = &nrn_GL,
        .VL = &nrn_VL,
        .GNa = &nrn_GNa,
        .VNa = &nrn_VNa,
        .GK = &nrn_GK,
        .VK = &nrn_VK,
        .C = &nrn_C,
        .pre_syn_size = p_in_size,
        .I = I,
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
        .pos_1 = &pos_1,
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
    arpra_uint sys_grps = 6;
    arpra_uint sys_dims[6] = {
        p_nrn_size, p_nrn_size,
        p_nrn_size, p_nrn_size,
        p_syn_size, p_syn_size,
    };
    arpra_ode_f sys_f[6] = {
        dMdt, dHdt, dNdt, dVdt, dRdt, dSdt,
    };
    void *sys_params[6] = {
        &params_nrn_M, &params_nrn_H,
        &params_nrn_N, &params_nrn_V,
        &params_syn_R, &params_syn_S,
    };
    arpra_range *sys_x[6] = {
        nrn_M, nrn_H, nrn_N, nrn_V, syn_R, syn_S,
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
            nrn_M_reduce_epoch[j] = ode_system.x[grp_nrn_M][j].nTerms;
            nrn_H_reduce_epoch[j] = ode_system.x[grp_nrn_H][j].nTerms;
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
            reduce_n = ode_system.x[grp_nrn_M][j].nTerms - nrn_M_reduce_epoch[j];
            arpra_reduce_last_n(&(ode_system.x[grp_nrn_M][j]), &(ode_system.x[grp_nrn_M][j]), reduce_n);
            reduce_n = ode_system.x[grp_nrn_H][j].nTerms - nrn_H_reduce_epoch[j];
            arpra_reduce_last_n(&(ode_system.x[grp_nrn_H][j]), &(ode_system.x[grp_nrn_H][j]), reduce_n);
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
                arpra_reduce_small_rel(&(ode_system.x[grp_nrn_M][j]), &(ode_system.x[grp_nrn_M][j]), reduce_rel);
                arpra_reduce_small_rel(&(ode_system.x[grp_nrn_H][j]), &(ode_system.x[grp_nrn_H][j]), reduce_rel);
                arpra_reduce_small_rel(&(ode_system.x[grp_nrn_N][j]), &(ode_system.x[grp_nrn_N][j]), reduce_rel);
                arpra_reduce_small_rel(&(ode_system.x[grp_nrn_V][j]), &(ode_system.x[grp_nrn_V][j]), reduce_rel);
            }
            for (j = 0; j < p_syn_size; j++) {
                arpra_reduce_small_rel(&(ode_system.x[grp_syn_R][j]), &(ode_system.x[grp_syn_R][j]), reduce_rel);
                arpra_reduce_small_rel(&(ode_system.x[grp_syn_S][j]), &(ode_system.x[grp_syn_S][j]), reduce_rel);
            }
        }

        file_write(&sys_t, 1, f_time_c, f_time_r, f_time_n, f_time_s, f_time_d);

        file_write(nrn_M, p_nrn_size, f_nrn_M_c, f_nrn_M_r, f_nrn_M_n, f_nrn_M_s, f_nrn_M_d);
        file_write(nrn_H, p_nrn_size, f_nrn_H_c, f_nrn_H_r, f_nrn_H_n, f_nrn_H_s, f_nrn_H_d);
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
        arpra_clear(&(nrn_M[i]));
        arpra_clear(&(nrn_H[i]));
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
    arpra_clear(&nrn_GNa);
    arpra_clear(&nrn_GK);
    arpra_clear(&nrn_VL);
    arpra_clear(&nrn_VNa);
    arpra_clear(&nrn_VK);
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
    arpra_clear(&pos_1);
    arpra_clear(&pos_4);
    arpra_clear(&pos_5);
    arpra_clear(&pos_18);
    arpra_clear(&pos_25);
    arpra_clear(&neg_25);
    arpra_clear(&pos_40);
    arpra_clear(&neg_48);
    arpra_clear(&neg_50);
    arpra_clear(&neg_52);
    arpra_clear(&neg_55);
    arpra_clear(&pos_028);
    arpra_clear(&pos_032);
    arpra_clear(&pos_0128);
    arpra_clear(&pos_0032);
    arpra_clear(&pos_05);

    // Clear scratch space
    mpfr_clear(rand_uf);
    mpfr_clear(rand_nf);
    arpra_clear(&temp1);
    arpra_clear(&temp2);
    arpra_clear(&_a);
    arpra_clear(&_b);
    for (i = 0; i < p_in_size; i++) {
        arpra_clear(&(I[i]));
    }

    // Free system state
    free(nrn_M);
    free(nrn_H);
    free(nrn_N);
    free(nrn_V);
    free(syn_R);
    free(syn_S);

    // Free other arrays
    free(syn_GSyn);
    free(I);
    free(in);
    free(nrn_M_reduce_epoch);
    free(nrn_H_reduce_epoch);
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

    file_clear(p_nrn_size, f_nrn_M_c, f_nrn_M_r, f_nrn_M_n, f_nrn_M_s, f_nrn_M_d);
    free(f_nrn_M_c);
    free(f_nrn_M_r);
    free(f_nrn_M_n);
    free(f_nrn_M_s);
    free(f_nrn_M_d);
    file_clear(p_nrn_size, f_nrn_H_c, f_nrn_H_r, f_nrn_H_n, f_nrn_H_s, f_nrn_H_d);
    free(f_nrn_H_c);
    free(f_nrn_H_r);
    free(f_nrn_H_n);
    free(f_nrn_H_s);
    free(f_nrn_H_d);
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
