/*
 * morris_lecar_mpfi.c -- Test Morris-Lecar model using MPFI.
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
const double p_h = 0.5;
const double p_t0 = 0.0;
const mp_prec_t p_prec = 53;
const mp_prec_t p_error_prec = 256;
const unsigned long p_sim_steps = 1000;
const unsigned long p_report_step = 20;

// Poisson input parameters (group 1)
const unsigned long p_in1_size = 50;
const double p_in1_freq = 5.0;
const double p_in1_V_lo = -60.0;
const double p_in1_V_hi = 20.0;

// Poisson input parameters (group 2)
const unsigned long p_in2_size = 0;
const double p_in2_freq = 5.0;
const double p_in2_V_lo = -60.0;
const double p_in2_V_hi = 20.0;

// Neuron parameters (group 1)
const unsigned long p_nrn1_size = 1;
const double p_nrn1_N0 = 0.0;
const double p_nrn1_V0 = -60.0;
const int p_nrn1_class = 1;

// Neuron parameters (group 2)
const unsigned long p_nrn2_size = 0;
const double p_nrn2_N0 = 0.0;
const double p_nrn2_V0 = -60.0;
const int p_nrn2_class = 1;

// Neuron parameters (common)
const double p_GL = 2.0;
const double p_GCa = 4.0; // Class 1
//const double p_GCa = 4.4; // Class 2
const double p_GK = 8.0;
const double p_VL = -60.0;
const double p_VCa = 120.0;
const double p_VK = -80.0;
const double p_V1 = -1.2;
const double p_V2 = 18.0;
const double p_V3 = 12.0; // Class 1
//const double p_V3 = 2.0; // Class 2
const double p_V4 = 17.4; // Class 1
//const double p_V4 = 30.0; // Class 2
const double p_phi = 1.0 / 15.0; // Class 1
//const double p_phi = 1.0 / 25.0; // Class 2
const double p_C = 20.0;

// Synapse parameters (excitatory)
const unsigned long p_syn_exc_size = p_in1_size * p_nrn1_size;
const double p_syn_exc_R0 = 0.0;
const double p_syn_exc_S0 = 0.0;
//const double p_syn_exc_GSyn = 40.0;
const double p_syn_exc_GSyn = 3.0;
const double p_syn_exc_VSyn = 0.0;
const double p_syn_exc_thr = -50.0;
const double p_syn_exc_a = 0.25; // in [1/10, 1/2]
const double p_syn_exc_b = 0.15; // in [1/20, 1/4]
const double p_syn_exc_k = 1.0E6;

// Synapse parameters (inhibitory)
const unsigned long p_syn_inh_size = 0;
const double p_syn_inh_R0 = 0.0;
const double p_syn_inh_S0 = 0.0;
const double p_syn_inh_GSyn = 50.0;
const double p_syn_inh_VSyn = -80.0;
const double p_syn_inh_thr = -50.0;
const double p_syn_inh_a = 0.075; // in [1/20, 1/10]
const double p_syn_inh_b = 0.035; // in [1/50, 1/20]
const double p_syn_inh_k = 1.0E6;

// ===================== end of model parameters ======================


int *in1, *in2;
mpfr_t rand_f, temp_error, in1_p0, in2_p0;
mpfr_ptr temp_sum_error1, *temp_sum_error1_ptr, temp_sum_error2, *temp_sum_error2_ptr;
mpfi_t GL, VL, GCa, VCa, GK, VK, V1, V2, V3, V4, phi, C, syn_exc_VSyn, syn_exc_thr,
       syn_exc_a, syn_exc_b, syn_exc_k, syn_inh_VSyn, syn_inh_thr, syn_inh_a,
       syn_inh_b, syn_inh_k, one, two, neg_two, temp_sum, temp1, temp2, M_ss, N_ss,
       in1_V_lo, in1_V_hi, in2_V_lo, in2_V_hi;
mpfi_ptr syn_exc_GSyn, syn_inh_GSyn, I1, I2;
gmp_randstate_t rng_f;
unsigned long rng_f_seed;

// System state variables
mpfi_ptr nrn1_N, nrn2_N, nrn1_V, nrn2_V, syn_exc_R, syn_inh_R, syn_exc_S, syn_inh_S;
mpfi_ptr d_nrn1_N, d_nrn2_N, d_nrn1_V, d_nrn2_V, d_syn_exc_R, d_syn_inh_R, d_syn_exc_S, d_syn_inh_S;

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

void file_write (mpfi_srcptr A, unsigned long grp_size, FILE **f)
{
    unsigned long i, j;

    for (i = 0; i < grp_size; i++) {
        mpfr_out_str(f[i], 10, 40, &(A[i].left), MPFR_RNDN);
        fputs(" ", f[i]);
        mpfr_out_str(f[i], 10, 40, &(A[i].right), MPFR_RNDN);
        fputs("\n", f[i]);
    }
}

void dNdt (const unsigned long idx, int grp)
{
    mpfi_ptr d_N;
    mpfi_srcptr N, V;

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

void dVdt (const unsigned long idx, int grp)
{
    mpfi_ptr d_V;
    mpfi_srcptr N, V;
    mpfi_srcptr S, GSyn, VSyn;
    mpfr_ptr temp_sum_error, *temp_sum_error_ptr;
    mpfi_ptr I;
    unsigned long i, pre_size;

    if (grp == 1) {
        d_V = d_nrn1_V + idx;
        N = nrn1_N + idx;
        V = nrn1_V + idx;
        pre_size = p_in1_size;
        S = syn_exc_S + (idx * pre_size);
        GSyn = syn_exc_GSyn + (idx * pre_size);
        VSyn = syn_exc_VSyn;
        temp_sum_error = temp_sum_error1;
        temp_sum_error_ptr = temp_sum_error1_ptr;
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
    //   temp_sum_error = temp_sum_error2;
    //   temp_sum_error_ptr = temp_sum_error2_ptr;
    //   I = I2;
    //}

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

void dRdt (const unsigned long idx, int grp)
{
    mpfi_ptr d_R;
    mpfi_srcptr R;
    mpfi_srcptr a, b, k;
    mpfi_srcptr VPre, threshold;

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

void dSdt (const unsigned long idx, int grp)
{
    mpfi_ptr d_S;
    mpfi_srcptr R, S;
    mpfi_srcptr a, b;

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
    struct timespec sys_t;
    clock_t run_time;
    unsigned long i, j;

    // Allocate dynamic arrays
    nrn1_N = malloc(p_nrn1_size * sizeof(mpfi_t));
    nrn2_N = malloc(p_nrn2_size * sizeof(mpfi_t));
    nrn1_V = malloc(p_nrn1_size * sizeof(mpfi_t));
    nrn2_V = malloc(p_nrn2_size * sizeof(mpfi_t));
    syn_exc_R = malloc(p_syn_exc_size * sizeof(mpfi_t));
    syn_inh_R = malloc(p_syn_inh_size * sizeof(mpfi_t));
    syn_exc_S = malloc(p_syn_exc_size * sizeof(mpfi_t));
    syn_inh_S = malloc(p_syn_inh_size * sizeof(mpfi_t));
    d_nrn1_N = malloc(p_nrn1_size * sizeof(mpfi_t));
    d_nrn2_N = malloc(p_nrn2_size * sizeof(mpfi_t));
    d_nrn1_V = malloc(p_nrn1_size * sizeof(mpfi_t));
    d_nrn2_V = malloc(p_nrn2_size * sizeof(mpfi_t));
    d_syn_exc_R = malloc(p_syn_exc_size * sizeof(mpfi_t));
    d_syn_inh_R = malloc(p_syn_inh_size * sizeof(mpfi_t));
    d_syn_exc_S = malloc(p_syn_exc_size * sizeof(mpfi_t));
    d_syn_inh_S = malloc(p_syn_inh_size * sizeof(mpfi_t));
    syn_exc_GSyn = malloc(p_syn_exc_size * sizeof(mpfi_t));
    syn_inh_GSyn = malloc(p_syn_inh_size * sizeof(mpfi_t));
    temp_sum_error1 = malloc(p_in1_size * sizeof(mpfr_t));
    temp_sum_error2 = malloc(p_in2_size * sizeof(mpfr_t));
    temp_sum_error1_ptr = malloc(p_in1_size * sizeof(mpfr_ptr));
    temp_sum_error2_ptr = malloc(p_in2_size * sizeof(mpfr_ptr));
    I1 = malloc(p_in1_size * sizeof(mpfi_t));
    I2 = malloc(p_in2_size * sizeof(mpfi_t));
    in1 = malloc(p_in1_size * sizeof(int));
    in2 = malloc(p_in2_size * sizeof(int));

    // Initialise system state
    mpfi_init2(h, p_prec);
    mpfi_init2(t, p_prec);
    for (i = 0; i < p_nrn1_size; i++) {
        mpfi_init2(&(nrn1_N[i]), p_prec);
        mpfi_init2(&(nrn1_V[i]), p_prec);
        mpfi_init2(&(d_nrn1_N[i]), p_prec);
        mpfi_init2(&(d_nrn1_V[i]), p_prec);
    }
    for (i = 0; i < p_nrn2_size; i++) {
        mpfi_init2(&(nrn2_N[i]), p_prec);
        mpfi_init2(&(nrn2_V[i]), p_prec);
        mpfi_init2(&(d_nrn2_N[i]), p_prec);
        mpfi_init2(&(d_nrn2_V[i]), p_prec);
    }
    for (i = 0; i < p_syn_exc_size; i++) {
        mpfi_init2(&(syn_exc_R[i]), p_prec);
        mpfi_init2(&(syn_exc_S[i]), p_prec);
        mpfi_init2(&(d_syn_exc_R[i]), p_prec);
        mpfi_init2(&(d_syn_exc_S[i]), p_prec);
    }
    for (i = 0; i < p_syn_inh_size; i++) {
        mpfi_init2(&(syn_inh_R[i]), p_prec);
        mpfi_init2(&(syn_inh_S[i]), p_prec);
        mpfi_init2(&(d_syn_inh_R[i]), p_prec);
        mpfi_init2(&(d_syn_inh_S[i]), p_prec);
    }

    // Initialise Poisson input parameters (group 1)
    mpfr_init2(in1_p0, p_prec);
    mpfi_init2(in1_V_lo, p_prec);
    mpfi_init2(in1_V_hi, p_prec);

    // Initialise Poisson input parameters (group 2)
    mpfr_init2(in2_p0, p_prec);
    mpfi_init2(in2_V_lo, p_prec);
    mpfi_init2(in2_V_hi, p_prec);

    // Initialise neuron parameters
    mpfi_init2(GL, p_prec);
    mpfi_init2(GCa, p_prec);
    mpfi_init2(GK, p_prec);
    mpfi_init2(VL, p_prec);
    mpfi_init2(VCa, p_prec);
    mpfi_init2(VK, p_prec);
    mpfi_init2(V1, p_prec);
    mpfi_init2(V2, p_prec);
    mpfi_init2(V3, p_prec);
    mpfi_init2(V4, p_prec);
    mpfi_init2(phi, p_prec);
    mpfi_init2(C, p_prec);

    // Initialise excitatory synapse parameters
    for (i = 0; i < p_syn_exc_size; i++) {
        mpfi_init2(&(syn_exc_GSyn[i]), p_prec);
    }
    mpfi_init2(syn_exc_VSyn, p_prec);
    mpfi_init2(syn_exc_thr, p_prec);
    mpfi_init2(syn_exc_a, p_prec);
    mpfi_init2(syn_exc_b, p_prec);
    mpfi_init2(syn_exc_k, p_prec);

    // Initialise inhibitory synapse parameters
    for (i = 0; i < p_syn_inh_size; i++) {
        mpfi_init2(&(syn_inh_GSyn[i]), p_prec);
    }
    mpfi_init2(syn_inh_VSyn, p_prec);
    mpfi_init2(syn_inh_thr, p_prec);
    mpfi_init2(syn_inh_a, p_prec);
    mpfi_init2(syn_inh_b, p_prec);
    mpfi_init2(syn_inh_k, p_prec);

    // Initialise constants
    mpfi_init2(one, p_prec);
    mpfi_init2(two, p_prec);
    mpfi_init2(neg_two, p_prec);

    // Initialise scratch space
    mpfr_init2(rand_f, p_prec);
    mpfr_init2(temp_error, p_error_prec);
    mpfi_init2(temp_sum, p_error_prec);
    mpfi_init2(temp1, p_prec);
    mpfi_init2(temp2, p_prec);
    mpfi_init2(M_ss, p_prec);
    mpfi_init2(N_ss, p_prec);

    // Set system state
    mpfi_set_d(h, p_h);
    mpfi_set_d(t, p_t0);
    for (i = 0; i < p_nrn1_size; i++) {
        mpfi_set_d(&(nrn1_N[i]), p_nrn1_N0);
        mpfi_set_d(&(nrn1_V[i]), p_nrn1_V0);
    }
    for (i = 0; i < p_nrn2_size; i++) {
        mpfi_set_d(&(nrn2_N[i]), p_nrn2_N0);
        mpfi_set_d(&(nrn2_V[i]), p_nrn2_V0);
    }
    for (i = 0; i < p_syn_exc_size; i++) {
        mpfi_set_d(&(syn_exc_R[i]), p_syn_exc_R0);
        mpfi_set_d(&(syn_exc_S[i]), p_syn_exc_S0);
    }
    for (i = 0; i < p_syn_inh_size; i++) {
        mpfi_set_d(&(syn_inh_R[i]), p_syn_inh_R0);
        mpfi_set_d(&(syn_inh_S[i]), p_syn_inh_S0);
    }

    // Set Poisson input parameters (group 1)
    mpfr_set_d(in1_p0, -(p_in1_freq / 1000) * p_h, MPFR_RNDN);
    mpfr_exp(in1_p0, in1_p0, MPFR_RNDN);
    mpfi_set_d(in1_V_lo, p_in1_V_lo);
    mpfi_set_d(in1_V_hi, p_in1_V_hi);
    for (i = 0; i < p_in1_size; i++) {
        mpfr_init2(&(temp_sum_error1[i]), p_error_prec);
        temp_sum_error1_ptr[i] = &(temp_sum_error1[i]);
        mpfi_init2(&(I1[i]), p_prec);
    }

    // Set Poisson input parameters (group 2)
    mpfr_set_d(in2_p0, -(p_in2_freq / 1000) * p_h, MPFR_RNDN);
    mpfr_exp(in2_p0, in2_p0, MPFR_RNDN);
    mpfi_set_d(in2_V_lo, p_in2_V_lo);
    mpfi_set_d(in2_V_hi, p_in2_V_hi);
    for (i = 0; i < p_in2_size; i++) {
        mpfr_init2(&(temp_sum_error2[i]), p_error_prec);
        temp_sum_error2_ptr[i] = &(temp_sum_error2[i]);
        mpfi_init2(&(I2[i]), p_prec);
    }

    // Set neuron parameters
    mpfi_set_d(GL, p_GL);
    mpfi_set_d(GCa, p_GCa);
    mpfi_set_d(GK, p_GK);
    mpfi_set_d(VL, p_VL);
    mpfi_set_d(VCa, p_VCa);
    mpfi_set_d(VK, p_VK);
    mpfi_set_d(V1, p_V1);
    mpfi_set_d(V2, p_V2);
    mpfi_set_d(V3, p_V3);
    mpfi_set_d(V4, p_V4);
    mpfi_set_d(phi, p_phi);
    mpfi_set_d(C, p_C);

    // Set excitatory synapse parameters
    for (i = 0; i < p_syn_exc_size; i++) {
        mpfi_set_d(&(syn_exc_GSyn[i]), p_syn_exc_GSyn);
    }
    mpfi_set_d(syn_exc_VSyn, p_syn_exc_VSyn);
    mpfi_set_d(syn_exc_thr, p_syn_exc_thr);
    mpfi_set_d(syn_exc_a, p_syn_exc_a);
    mpfi_set_d(syn_exc_b, p_syn_exc_b);
    mpfi_set_d(syn_exc_k, p_syn_exc_k);
    mpfi_neg(syn_exc_k, syn_exc_k);

    // Set inhibitory synapse parameters
    for (i = 0; i < p_syn_inh_size; i++) {
        mpfi_set_d(&(syn_inh_GSyn[i]), p_syn_inh_GSyn);
    }
    mpfi_set_d(syn_inh_VSyn, p_syn_inh_VSyn);
    mpfi_set_d(syn_inh_thr, p_syn_inh_thr);
    mpfi_set_d(syn_inh_a, p_syn_inh_a);
    mpfi_set_d(syn_inh_b, p_syn_inh_b);
    mpfi_set_d(syn_inh_k, p_syn_inh_k);
    mpfi_neg(syn_inh_k, syn_inh_k);

    // Set constants
    mpfi_set_d(one, 1.0);
    mpfi_set_d(two, 2.0);
    mpfi_set_d(neg_two, -2.0);

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

    // Initialise RNG
    gmp_randinit_default(rng_f);
    clock_gettime(CLOCK_REALTIME, &sys_t);
    //rng_f_seed = 707135875931353ul;
    rng_f_seed = sys_t.tv_sec + sys_t.tv_nsec;
    gmp_randseed_ui(rng_f, rng_f_seed);
    printf("GMP rand seed: %lu\n", rng_f_seed);


    // Begin simulation loop
    // =====================

    run_time = clock();

    for (i = 0; i < p_sim_steps; i++) {
        if (i % p_report_step == 0) printf("%lu\n", i);

        // Event(s) occur if urandom >= e^-rate
        for (j = 0; j < p_in1_size; j++) {
            mpfr_urandom(rand_f, rng_f, MPFR_RNDN);
            in1[j] = mpfr_greaterequal_p(rand_f, in1_p0);
            fprintf(stderr, "%s", (in1[j] ? "\x1B[31m\xE2\x96\xA3\x1B[0m" : "\xE2\x96\xA3"));
        }
        fprintf(stderr, "  ");
        for (j = 0; j < p_in2_size; j++) {
            mpfr_urandom(rand_f, rng_f, MPFR_RNDN);
            in2[j] = mpfr_greaterequal_p(rand_f, in2_p0);
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
        mpfi_add(t, t, h);
        for (j = 0; j < p_nrn1_size; j++) {
            mpfi_mul(&(d_nrn1_N[j]), &(d_nrn1_N[j]), h);
            mpfi_add(&(nrn1_N[j]), &(nrn1_N[j]), &(d_nrn1_N[j]));
            mpfi_mul(&(d_nrn1_V[j]), &(d_nrn1_V[j]), h);
            mpfi_add(&(nrn1_V[j]), &(nrn1_V[j]), &(d_nrn1_V[j]));
        }
        for (j = 0; j < p_nrn2_size; j++) {
            mpfi_mul(&(d_nrn2_N[j]), &(d_nrn2_N[j]), h);
            mpfi_add(&(nrn2_N[j]), &(nrn2_N[j]), &(d_nrn2_N[j]));
            mpfi_mul(&(d_nrn2_V[j]), &(d_nrn2_V[j]), h);
            mpfi_add(&(nrn2_V[j]), &(nrn2_V[j]), &(d_nrn2_V[j]));
        }
        for (j = 0; j < p_syn_exc_size; j++) {
            mpfi_mul(&(d_syn_exc_R[j]), &(d_syn_exc_R[j]), h);
            mpfi_add(&(syn_exc_R[j]), &(syn_exc_R[j]), &(d_syn_exc_R[j]));
            mpfi_mul(&(d_syn_exc_S[j]), &(d_syn_exc_S[j]), h);
            mpfi_add(&(syn_exc_S[j]), &(syn_exc_S[j]), &(d_syn_exc_S[j]));
        }
        for (j = 0; j < p_syn_inh_size; j++) {
            mpfi_mul(&(d_syn_inh_R[j]), &(d_syn_inh_R[j]), h);
            mpfi_add(&(syn_inh_R[j]), &(syn_inh_R[j]), &(d_syn_inh_R[j]));
            mpfi_mul(&(d_syn_inh_S[j]), &(d_syn_inh_S[j]), h);
            mpfi_add(&(syn_inh_S[j]), &(syn_inh_S[j]), &(d_syn_inh_S[j]));
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
    mpfi_clear(h);
    mpfi_clear(t);
    for (i = 0; i < p_nrn1_size; i++) {
        mpfi_clear(&(nrn1_N[i]));
        mpfi_clear(&(nrn1_V[i]));
        mpfi_clear(&(d_nrn1_N[i]));
        mpfi_clear(&(d_nrn1_V[i]));
    }
    for (i = 0; i < p_nrn2_size; i++) {
        mpfi_clear(&(nrn2_N[i]));
        mpfi_clear(&(nrn2_V[i]));
        mpfi_clear(&(d_nrn2_N[i]));
        mpfi_clear(&(d_nrn2_V[i]));
    }
    for (i = 0; i < p_syn_exc_size; i++) {
        mpfi_clear(&(syn_exc_R[i]));
        mpfi_clear(&(syn_exc_S[i]));
        mpfi_clear(&(d_syn_exc_R[i]));
        mpfi_clear(&(d_syn_exc_S[i]));
    }
    for (i = 0; i < p_syn_inh_size; i++) {
        mpfi_clear(&(syn_inh_R[i]));
        mpfi_clear(&(syn_inh_S[i]));
        mpfi_clear(&(d_syn_inh_R[i]));
        mpfi_clear(&(d_syn_inh_S[i]));
    }

    // Clear Poisson input parameters (group 1)
    mpfr_clear(in1_p0);
    mpfi_clear(in1_V_lo);
    mpfi_clear(in1_V_hi);
    for (i = 0; i < p_in1_size; i++) {
        mpfr_clear(&(temp_sum_error1[i]));
        mpfi_clear(&(I1[i]));
    }

    // Clear Poisson input parameters (group 2)
    mpfr_clear(in2_p0);
    mpfi_clear(in2_V_lo);
    mpfi_clear(in2_V_hi);
    for (i = 0; i < p_in2_size; i++) {
        mpfr_clear(&(temp_sum_error2[i]));
        mpfi_clear(&(I2[i]));
    }

    // Clear neuron parameters
    mpfi_clear(GL);
    mpfi_clear(GCa);
    mpfi_clear(GK);
    mpfi_clear(VL);
    mpfi_clear(VCa);
    mpfi_clear(VK);
    mpfi_clear(V1);
    mpfi_clear(V2);
    mpfi_clear(V3);
    mpfi_clear(V4);
    mpfi_clear(phi);
    mpfi_clear(C);

    // Clear excitatory synapse parameters
    for (i = 0; i < p_syn_exc_size; i++) {
        mpfi_clear(&(syn_exc_GSyn[i]));
    }
    mpfi_clear(syn_exc_VSyn);
    mpfi_clear(syn_exc_thr);
    mpfi_clear(syn_exc_a);
    mpfi_clear(syn_exc_b);
    mpfi_clear(syn_exc_k);

    // Clear inhibitory synapse parameters
    for (i = 0; i < p_syn_inh_size; i++) {
        mpfi_clear(&(syn_inh_GSyn[i]));
    }
    mpfi_clear(syn_inh_VSyn);
    mpfi_clear(syn_inh_thr);
    mpfi_clear(syn_inh_a);
    mpfi_clear(syn_inh_b);
    mpfi_clear(syn_inh_k);

    // Clear constants
    mpfi_clear(one);
    mpfi_clear(two);
    mpfi_clear(neg_two);

    // Clear scratch space
    mpfr_clear(rand_f);
    mpfr_clear(temp_error);
    mpfi_clear(temp_sum);
    mpfi_clear(temp1);
    mpfi_clear(temp2);
    mpfi_clear(M_ss);
    mpfi_clear(N_ss);

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
    free(temp_sum_error1);
    free(temp_sum_error2);
    free(temp_sum_error1_ptr);
    free(temp_sum_error2_ptr);
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

    gmp_randclear(rng_f);
    mpfr_free_cache();

    return 0;
}
