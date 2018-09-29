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
 * Synapse parameters
 * ------------------
 * GSyn   Maximum synapse conductance (mS/cm^2)
 * VSyn   Synapse conductance equilibrium potential (mV)
 * thr    Neuronal spike threshold (mV)
 * a      Transmitter release/bind rise factor
 * b      Transmitter release/bind decay factor
 * k      Negated steepness of activation function
 */




// =========== DOCUMENT MODE VARIABLES ALSO

// N: Fraction of open K+ channels
// V: Membrane potential
// R: Transmitter release
// S: Transmitter binding





// General parameters
const double p_h0 = 0.5;
const double p_t0 = 0.0;
const double p_reduce_ratio = 0.3;
const arpra_precision p_prec = 53;
const arpra_uint p_sim_steps = 10;
const arpra_uint p_report_step = 20;
const arpra_uint p_reduce_step = 50;

// Poisson input parameters
const arpra_uint p_in1_size = 5;
const arpra_uint p_in2_size = 0;
const double p_in_V_lo = -60.0;
const double p_in_V_hi = 20.0;

// Neuron parameters (group 1)
const arpra_uint p_nrn1_size = 5;
const double p_nrn1_N0 = 0.0;
const double p_nrn1_V0 = -60.0;
const int p_nrn1_class = 1;

// Neuron parameters (group 2)
const arpra_uint p_nrn2_size = 0;
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
const arpra_uint p_syn_exc_size = p_nrn1_size * p_nrn1_size;
const double p_syn_exc_R0 = 0.0;
const double p_syn_exc_S0 = 0.0;
const double p_syn_exc_GSyn = 0.1;
const double p_syn_exc_VSyn = 0.0;
const double p_syn_exc_thr = -50.0;
const double p_syn_exc_a = 0.25; // in [1/10, 1/2]
const double p_syn_exc_b = 0.15; // in [1/20, 1/4]
const double p_syn_exc_k = -1.0E6;

// Synapse parameters (inhibitory)
const arpra_uint p_syn_inh_size = 0;
const double p_syn_inh_R0 = 0.0;
const double p_syn_inh_S0 = 0.0;
const double p_syn_inh_GSyn = 0.1;
const double p_syn_inh_VSyn = -80.0;
const double p_syn_inh_thr = -50.0;
const double p_syn_inh_a = 0.075; // in [1/20, 1/10]
const double p_syn_inh_b = 0.035; // in [1/50, 1/20]
const double p_syn_inh_k = -1.0E6;

// ===================== end of model parameters ======================


int *in1, *in2;
arpra_range in_V_lo, in_V_hi, GL, VL, GCa, VCa, GK, VK, V1, V2, V3, V4,
    phi, C, *syn_exc_GSyn, syn_exc_VSyn, syn_exc_thr, syn_exc_a,
    syn_exc_b, syn_exc_k, *syn_inh_GSyn, syn_inh_VSyn, syn_inh_thr,
    syn_inh_a, syn_inh_b, syn_inh_k, one, two, neg_two, temp1, temp2,
    M_ss, N_ss, *I1, *I2;

// State memory offsets
const arpra_uint nrn1_N_offset = 0;
const arpra_uint nrn2_N_offset = nrn1_N_offset + p_nrn1_size;
const arpra_uint nrn1_V_offset = nrn2_N_offset + p_nrn2_size;
const arpra_uint nrn2_V_offset = nrn1_V_offset + p_nrn1_size;
const arpra_uint syn_exc_R_offset = nrn2_V_offset + p_nrn2_size;
const arpra_uint syn_inh_R_offset = syn_exc_R_offset + p_syn_exc_size;
const arpra_uint syn_exc_S_offset = syn_inh_R_offset + p_syn_inh_size;
const arpra_uint syn_inh_S_offset = syn_exc_S_offset + p_syn_exc_size;
const arpra_uint dimensions = syn_inh_S_offset + p_syn_inh_size;

// State indexing macros
#define nrn1_N (x + nrn1_N_offset)
#define nrn2_N (x + nrn2_N_offset)
#define nrn1_V (x + nrn1_V_offset)
#define nrn2_V (x + nrn2_V_offset)
#define syn_exc_R (x + syn_exc_R_offset)
#define syn_inh_R (x + syn_inh_R_offset)
#define syn_exc_S (x + syn_exc_S_offset)
#define syn_inh_S (x + syn_inh_S_offset)

// DEBUG: print MPFR numbers to stderr
void debug (const arpra_mpfr x) {
    mpfr_out_str(stderr, 10, 80, &x, MPFR_RNDN);
    fputs("\n", stderr);
}


void file_init (char *grp, char *var, arpra_uint grp_size,
                FILE **c, FILE **r, FILE **n, FILE **s, FILE **d)
{
    char fname[20];
    arpra_uint i;

    for (i = 0; i < grp_size; i++) {
        sprintf(fname, "%s_%s_%03u_c.dat", grp, var, (unsigned) i);
        c[i] = fopen(fname, "w");
        sprintf(fname, "%s_%s_%03u_r.dat", grp, var, (unsigned) i);
        r[i] = fopen(fname, "w");
        sprintf(fname, "%s_%s_%03u_n.dat", grp, var, (unsigned) i);
        n[i] = fopen(fname, "w");
        sprintf(fname, "%s_%s_%03u_s.dat", grp, var, (unsigned) i);
        s[i] = fopen(fname, "w");
        sprintf(fname, "%s_%s_%03u_d.dat", grp, var, (unsigned) i);
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
        fprintf(n[i], "%u\n", (unsigned) A[i].nTerms);
        for (j = 0; j < A[i].nTerms; j++) {
            fprintf(s[i], "%u ", (unsigned) A[i].symbols[j]);
            mpfr_out_str(d[i], 10, 80, &(A[i].deviations[j]), MPFR_RNDN);
            fputc(' ', d[i]);
        }
        fputc('\n', s[i]);
        fputc('\n', d[i]);
    }
}

void dxdt (arpra_range *out,
           const arpra_range *t, const arpra_range *x,
           const arpra_uint x_idx, const void *params)
{
    arpra_uint idx;

    // dN/dt
    if (x_idx < nrn1_V_offset) {
        const arpra_range *N, *V;

        if (x_idx < nrn2_N_offset) {
            idx = x_idx;
            N = nrn1_N + idx;
            V = nrn1_V + idx;
        }
        //else {
        //    idx = x_idx % nrn2_N_offset;
        //    N = nrn2_N + idx;
        //    V = nrn2_V + idx;
        //}

        // K+ channel activation steady-state
        // N_ss = 1 / (1 + exp(-2 (V - V3) / V4))
        arpra_sub(&temp1, V, &V3);
        arpra_mul(&N_ss, &neg_two, &temp1);
        arpra_div(&N_ss, &N_ss, &V4);
        arpra_exp(&N_ss, &N_ss);
        arpra_add(&N_ss, &one, &N_ss);
        arpra_inv(&N_ss, &N_ss);

        // tau of K+ channel activation
        // tau = 1 / (phi ((p + q) / 2))
        // p = exp(-(V - V3) / (2 V4))
        // q = exp( (V - V3) / (2 V4))
        arpra_mul(&temp2, &two, &V4);
        arpra_div(&temp2, &temp1, &temp2);
        arpra_neg(&temp1, &temp2);
        arpra_exp(&temp1, &temp1);
        arpra_exp(&temp2, &temp2);
        arpra_add(&temp1, &temp1, &temp2);
        arpra_div(&temp1, &temp1, &two);
        arpra_mul(&temp1, &phi, &temp1);
        arpra_inv(&temp1, &temp1);

        // delta of K+ channel activation
        // dN/dt = (N_ss - N) / tau
        arpra_sub(out, &N_ss, N);
        arpra_div(out, out, &temp1);
    }

    // dV/dt
    else if (x_idx < syn_exc_R_offset) {
        const arpra_range *N, *V;
        const arpra_range *S, *GSyn, *VSyn;
        arpra_range *I;
        arpra_uint i, pre_size;

        if (x_idx < nrn2_V_offset) {
            idx = x_idx % nrn1_V_offset;
            N = nrn1_N + idx;
            V = nrn1_V + idx;
            pre_size = p_in1_size;
            S = syn_exc_S + (idx * pre_size);
            GSyn = syn_exc_GSyn + (idx * pre_size);
            VSyn = &syn_exc_VSyn;
            I = I1;
        }
        //else {
        //    idx = x_idx % nrn2_V_offset;
        //    N = nrn2_N + idx;
        //    V = nrn2_V + idx;
        //    pre_size = p_in2_size;
        //    S = syn_exc_S + (idx * pre_size);
        //    GSyn = syn_exc_GSyn + (idx * pre_size);
        //    VSyn = &syn_exc_VSyn;
        //    I = I2;
        //}

        // Ca++ channel activation steady-state
        // M_ss = 1 / (1 + exp(-2 (V - V1) / V2))
        arpra_sub(&M_ss, V, &V1);
        arpra_mul(&M_ss, &neg_two, &M_ss);
        arpra_div(&M_ss, &M_ss, &V2);
        arpra_exp(&M_ss, &M_ss);
        arpra_add(&M_ss, &one, &M_ss);
        arpra_inv(&M_ss, &M_ss);

        // Synapse current
        arpra_sub(&temp1, VSyn, V);
        for (i = 0; i < pre_size; i++) {
            arpra_mul(&(I[i]), &temp1, &(GSyn[i]));
            arpra_mul(&(I[i]), &(I[i]), &(S[i]));
        }
        arpra_sum(out, I, pre_size);




        // ======== TEMP DEBUG ==========
        // NOTE: bifurcation at sum(I) = 80.0
        //arpra_set_d(out, 80.0);
        if (idx == 0) {            
            //fprintf(stderr, "VSyn: "); debug(VSyn->centre);
            //fprintf(stderr, "V - VSyn: "); debug(temp1.centre);
            //fprintf(stderr, "I[0]: "); debug(I[0].centre);
            //fprintf(stderr, "sum(I): "); debug(out->centre);
        }



        // Leak current
        arpra_sub(&temp1, &VL, V);
        arpra_mul(&temp1, &temp1, &GL);
        arpra_add(out, out, &temp1);

        // Ca++ current
        arpra_sub(&temp1, &VCa, V);
        arpra_mul(&temp1, &temp1, &GCa);
        arpra_mul(&temp1, &temp1, &M_ss);
        arpra_add(out, out, &temp1);

        // K+ current
        arpra_sub(&temp1, &VK, V);
        arpra_mul(&temp1, &temp1, &GK);
        arpra_mul(&temp1, &temp1, N);
        arpra_add(out, out, &temp1);

        // delta of membrane potential
        // dV/dt = (I + GL (VL - V) + GCa M (VCa - V) + GK N (VK - V)) / C
        arpra_div(out, out, &C);
    }

    // dR/dt
    else if (x_idx < syn_exc_S_offset) {
        const arpra_range *R;
        const arpra_range *a, *b, *k;
        const arpra_range *VPre, *threshold;

        if (x_idx < syn_inh_R_offset) {
            idx = x_idx % syn_exc_R_offset;
            R = syn_exc_R + idx;
            a = &syn_exc_a;
            b = &syn_exc_b;
            k = &syn_exc_k;
            VPre = in1[idx % p_in1_size] ? &in_V_hi : &in_V_lo;
            threshold = &syn_exc_thr;
        }
        //else {
        //    idx = x_idx % syn_inh_R_offset;
        //    R = syn_inh_R + idx;
        //    a = &syn_inh_a;
        //    b = &syn_inh_b;
        //    k = &syn_inh_k;
        //    VPre = in2[idx % p_in2_size] ? &in_V_hi : &in_V_lo;
        //    threshold = &syn_inh_thr;
        //}

        // Sigmoid of threshold difference
        arpra_sub(&temp1, VPre, threshold);
        arpra_mul(&temp1, &temp1, k);
        arpra_exp(&temp1, &temp1);
        arpra_add(&temp1, &temp1, &one);
        arpra_inv(&temp1, &temp1);

        // Presynaptic transmitter release rise
        arpra_mul(&temp1, a, &temp1);

        // Presynaptic transmitter release decay
        arpra_mul(&temp2, b, R);

        // delta of presynaptic transmitter release
        // dR/dt = a t - b R
        // t = 1 / (1 + e^(k(V - threshold)))
        arpra_sub(out, &temp1, &temp2);
    }

    // dS/dt
    else {
        const arpra_range *R, *S;
        const arpra_range *a, *b;

        if (x_idx < syn_inh_S_offset) {
            idx = x_idx % syn_exc_S_offset;
            R = syn_exc_R + idx;
            S = syn_exc_S + idx;
            a = &syn_exc_a;
            b = &syn_exc_b;
        }
        //else {
        //    idx = x_idx % syn_inh_S_offset;
        //    R = syn_inh_R + idx;
        //    S = syn_inh_S + idx;
        //    a = &syn_inh_a;
        //    b = &syn_inh_b;
        //}

        // Postsynaptic transmitter binding rise
        arpra_mul(&temp1, a, R);

        // Postsynaptic transmitter binding decay
        arpra_mul(&temp2, b, S);

        // delta of postsynaptic transmitter binding
        // dS/dt = a R - b S
        arpra_sub(out, &temp1, &temp2);
    }
}

int main (int argc, char *argv[])
{
    arpra_range h, t, *x;
    arpra_uint *reduce_epoch;
    clock_t run_time;
    arpra_uint i, j;

    // Allocate dynamic arrays
    x = malloc(dimensions * sizeof(arpra_range));
    syn_exc_GSyn = malloc(p_syn_exc_size * sizeof(arpra_range));
    syn_inh_GSyn = malloc(p_syn_inh_size * sizeof(arpra_range));
    I1 = malloc(p_in1_size * sizeof(arpra_range));
    I2 = malloc(p_in2_size * sizeof(arpra_range));
    in1 = malloc(p_in1_size * sizeof(int));
    in2 = malloc(p_in2_size * sizeof(int));
    reduce_epoch = malloc(dimensions * sizeof(arpra_uint));

    // Initialise system state
    arpra_init2(&h, p_prec);
    arpra_init2(&t, p_prec);
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

    // Initialise Poisson input parameters
    arpra_init2(&in_V_lo, p_prec);
    arpra_init2(&in_V_hi, p_prec);

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
    arpra_set_d(&h, p_h0);
    arpra_set_d(&t, p_t0);
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

    // Set Poisson input parameters
    arpra_set_d(&in_V_lo, p_in_V_lo);
    arpra_set_d(&in_V_hi, p_in_V_hi);

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
        arpra_set_d(&(syn_exc_GSyn[i]), p_syn_exc_GSyn);
    }
    arpra_set_d(&syn_exc_VSyn, p_syn_exc_VSyn);
    arpra_set_d(&syn_exc_thr, p_syn_exc_thr);
    arpra_set_d(&syn_exc_a, p_syn_exc_a);
    arpra_set_d(&syn_exc_b, p_syn_exc_b);
    arpra_set_d(&syn_exc_k, p_syn_exc_k);

    // Set inhibitory synapse parameters
    for (i = 0; i < p_syn_inh_size; i++) {
        arpra_set_d(&(syn_inh_GSyn[i]), p_syn_inh_GSyn);
    }
    arpra_set_d(&syn_inh_VSyn, p_syn_inh_VSyn);
    arpra_set_d(&syn_inh_thr, p_syn_inh_thr);
    arpra_set_d(&syn_inh_a, p_syn_inh_a);
    arpra_set_d(&syn_inh_b, p_syn_inh_b);
    arpra_set_d(&syn_inh_k, p_syn_inh_k);

    // Set constants
    arpra_set_d(&one, 1.0);
    arpra_set_d(&two, 2.0);
    arpra_set_d(&neg_two, -2.0);

    // Initialise report files
    FILE **f_nrn1_N_c = malloc(p_nrn1_size * sizeof(FILE *));
    FILE **f_nrn1_N_r = malloc(p_nrn1_size * sizeof(FILE *));
    FILE **f_nrn1_N_n = malloc(p_nrn1_size * sizeof(FILE *));
    FILE **f_nrn1_N_s = malloc(p_nrn1_size * sizeof(FILE *));
    FILE **f_nrn1_N_d = malloc(p_nrn1_size * sizeof(FILE *));
    file_init("nrn1", "N", p_nrn1_size, f_nrn1_N_c, f_nrn1_N_r, f_nrn1_N_n, f_nrn1_N_s, f_nrn1_N_d);
    FILE **f_nrn1_V_c = malloc(p_nrn1_size * sizeof(FILE *));
    FILE **f_nrn1_V_r = malloc(p_nrn1_size * sizeof(FILE *));
    FILE **f_nrn1_V_n = malloc(p_nrn1_size * sizeof(FILE *));
    FILE **f_nrn1_V_s = malloc(p_nrn1_size * sizeof(FILE *));
    FILE **f_nrn1_V_d = malloc(p_nrn1_size * sizeof(FILE *));
    file_init("nrn1", "V", p_nrn1_size, f_nrn1_V_c, f_nrn1_V_r, f_nrn1_V_n, f_nrn1_V_s, f_nrn1_V_d);

    //FILE **f_nrn2_N_c = malloc(p_nrn2_size * sizeof(FILE *));
    //FILE **f_nrn2_N_r = malloc(p_nrn2_size * sizeof(FILE *));
    //FILE **f_nrn2_N_n = malloc(p_nrn2_size * sizeof(FILE *));
    //FILE **f_nrn2_N_s = malloc(p_nrn2_size * sizeof(FILE *));
    //FILE **f_nrn2_N_d = malloc(p_nrn2_size * sizeof(FILE *));
    //file_init("nrn2", "N", p_nrn2_size, f_nrn2_N_c, f_nrn2_N_r, f_nrn2_N_n, f_nrn2_N_s, f_nrn2_N_d);
    //FILE **f_nrn2_V_c = malloc(p_nrn2_size * sizeof(FILE *));
    //FILE **f_nrn2_V_r = malloc(p_nrn2_size * sizeof(FILE *));
    //FILE **f_nrn2_V_n = malloc(p_nrn2_size * sizeof(FILE *));
    //FILE **f_nrn2_V_s = malloc(p_nrn2_size * sizeof(FILE *));
    //FILE **f_nrn2_V_d = malloc(p_nrn2_size * sizeof(FILE *));
    //file_init("nrn2", "V", p_nrn2_size, f_nrn2_V_c, f_nrn2_V_r, f_nrn2_V_n, f_nrn2_V_s, f_nrn2_V_d);

    FILE **f_syn_exc_R_c = malloc(p_syn_exc_size * sizeof(FILE *));
    FILE **f_syn_exc_R_r = malloc(p_syn_exc_size * sizeof(FILE *));
    FILE **f_syn_exc_R_n = malloc(p_syn_exc_size * sizeof(FILE *));
    FILE **f_syn_exc_R_s = malloc(p_syn_exc_size * sizeof(FILE *));
    FILE **f_syn_exc_R_d = malloc(p_syn_exc_size * sizeof(FILE *));
    file_init("syn_exc", "R", p_syn_exc_size, f_syn_exc_R_c, f_syn_exc_R_r, f_syn_exc_R_n, f_syn_exc_R_s, f_syn_exc_R_d);
    FILE **f_syn_exc_S_c = malloc(p_syn_exc_size * sizeof(FILE *));
    FILE **f_syn_exc_S_r = malloc(p_syn_exc_size * sizeof(FILE *));
    FILE **f_syn_exc_S_n = malloc(p_syn_exc_size * sizeof(FILE *));
    FILE **f_syn_exc_S_s = malloc(p_syn_exc_size * sizeof(FILE *));
    FILE **f_syn_exc_S_d = malloc(p_syn_exc_size * sizeof(FILE *));
    file_init("syn_exc", "S", p_syn_exc_size, f_syn_exc_S_c, f_syn_exc_S_r, f_syn_exc_S_n, f_syn_exc_S_s, f_syn_exc_S_d);

    //FILE **f_syn_inh_R_c = malloc(p_syn_inh_size * sizeof(FILE *));
    //FILE **f_syn_inh_R_r = malloc(p_syn_inh_size * sizeof(FILE *));
    //FILE **f_syn_inh_R_n = malloc(p_syn_inh_size * sizeof(FILE *));
    //FILE **f_syn_inh_R_s = malloc(p_syn_inh_size * sizeof(FILE *));
    //FILE **f_syn_inh_R_d = malloc(p_syn_inh_size * sizeof(FILE *));
    //file_init("syn_inh", "R", p_syn_inh_size, f_syn_inh_R_c, f_syn_inh_R_r, f_syn_inh_R_n, f_syn_inh_R_s, f_syn_inh_R_d);
    //FILE **f_syn_inh_S_c = malloc(p_syn_inh_size * sizeof(FILE *));
    //FILE **f_syn_inh_S_r = malloc(p_syn_inh_size * sizeof(FILE *));
    //FILE **f_syn_inh_S_n = malloc(p_syn_inh_size * sizeof(FILE *));
    //FILE **f_syn_inh_S_s = malloc(p_syn_inh_size * sizeof(FILE *));
    //FILE **f_syn_inh_S_d = malloc(p_syn_inh_size * sizeof(FILE *));
    //file_init("syn_inh", "S", p_syn_inh_size, f_syn_inh_S_c, f_syn_inh_S_r, f_syn_inh_S_n, f_syn_inh_S_s, f_syn_inh_S_d);

    // ODE system
    arpra_ode_system ode_system;
    ode_system.f = dxdt;
    ode_system.t = &t;
    ode_system.x = x;
    ode_system.dims = dimensions;
    ode_system.params = NULL;

    // ODE stepper
    arpra_ode_stepper ode_stepper;
    //arpra_ode_stepper_init(&ode_stepper, &ode_system, arpra_ode_euler);
    //arpra_ode_stepper_init(&ode_stepper, &ode_system, arpra_ode_trapezoidal);
    arpra_ode_stepper_init(&ode_stepper, &ode_system, arpra_ode_bogsham32);
    //arpra_ode_stepper_init(&ode_stepper, &ode_system, arpra_ode_dopri54);
    //arpra_ode_stepper_init(&ode_stepper, &ode_system, arpra_ode_dopri87);


    // Begin simulation loop
    // =====================

    run_time = clock();

    for (i = 0; i < p_sim_steps; i++) {
        if (i % p_report_step == 0) printf("%u\n", i);

        for (j = 0; j < dimensions; j++) {
            reduce_epoch[j] = ode_system.x[j].nTerms;
        }




        // ===== COMPUTE POISSON IN1 AND IN2
        for (j = 0; j < p_in1_size; j++) {
            in1[j] = 1;
        }
        for (j = 0; j < p_in2_size; j++) {
            in2[j] = 0;
        }





        arpra_ode_stepper_step(&ode_stepper, &h);

        for (j = 0; j < dimensions; j++) {
            arpra_reduce_last_n(&(ode_system.x[j]), (ode_system.x[j].nTerms - reduce_epoch[j]));
            if (i % p_reduce_step == 0) {
                arpra_reduce_small(&(ode_system.x[j]), p_reduce_ratio);
            }
        }

        file_write(nrn1_N, p_nrn1_size, f_nrn1_N_c, f_nrn1_N_r, f_nrn1_N_n, f_nrn1_N_s, f_nrn1_N_d);
        file_write(nrn1_V, p_nrn1_size, f_nrn1_V_c, f_nrn1_V_r, f_nrn1_V_n, f_nrn1_V_s, f_nrn1_V_d);

        //file_write(nrn2_N, p_nrn2_size, f_nrn2_N_c, f_nrn2_N_r, f_nrn2_N_n, f_nrn2_N_s, f_nrn2_N_d);
        //file_write(nrn2_V, p_nrn2_size, f_nrn2_V_c, f_nrn2_V_r, f_nrn2_V_n, f_nrn2_V_s, f_nrn2_V_d);

        file_write(syn_exc_R, p_syn_exc_size, f_syn_exc_R_c, f_syn_exc_R_r, f_syn_exc_R_n, f_syn_exc_R_s, f_syn_exc_R_d);
        file_write(syn_exc_S, p_syn_exc_size, f_syn_exc_S_c, f_syn_exc_S_r, f_syn_exc_S_n, f_syn_exc_S_s, f_syn_exc_S_d);

        //file_write(syn_inh_R, p_syn_inh_size, f_syn_inh_R_c, f_syn_inh_R_r, f_syn_inh_R_n, f_syn_inh_R_s, f_syn_inh_R_d);
        //file_write(syn_inh_S, p_syn_inh_size, f_syn_inh_S_c, f_syn_inh_S_r, f_syn_inh_S_n, f_syn_inh_S_s, f_syn_inh_S_d);
    }

    run_time = clock() - run_time;
    printf("Finished in %f seconds.\n", ((float) run_time) / CLOCKS_PER_SEC);

    // End simulation loop
    // ===================


    // Clear system state
    arpra_clear(&h);
    arpra_clear(&t);
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

    // Clear Poisson input parameters
    arpra_clear(&in_V_lo);
    arpra_clear(&in_V_hi);

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

    // Free dynamic arrays
    free(x);
    free(syn_exc_GSyn);
    free(syn_inh_GSyn);
    free(I1);
    free(I2);
    free(in1);
    free(in2);
    free(reduce_epoch);

    // Clear report files
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

    //file_clear(p_nrn2_size, f_nrn2_N_c, f_nrn2_N_r, f_nrn2_N_n, f_nrn2_N_s, f_nrn2_N_d);
    //free(f_nrn2_N_c);
    //free(f_nrn2_N_r);
    //free(f_nrn2_N_n);
    //free(f_nrn2_N_s);
    //free(f_nrn2_N_d);
    //file_clear(p_nrn2_size, f_nrn2_V_c, f_nrn2_V_r, f_nrn2_V_n, f_nrn2_V_s, f_nrn2_V_d);
    //free(f_nrn2_V_c);
    //free(f_nrn2_V_r);
    //free(f_nrn2_V_n);
    //free(f_nrn2_V_s);
    //free(f_nrn2_V_d);

    file_clear(p_syn_exc_size, f_syn_exc_R_c, f_syn_exc_R_r, f_syn_exc_R_n, f_syn_exc_R_s, f_syn_exc_R_d);
    free(f_syn_exc_R_c);
    free(f_syn_exc_R_r);
    free(f_syn_exc_R_n);
    free(f_syn_exc_R_s);
    free(f_syn_exc_R_d);
    file_clear(p_syn_exc_size, f_syn_exc_S_c, f_syn_exc_S_r, f_syn_exc_S_n, f_syn_exc_S_s, f_syn_exc_S_d);
    free(f_syn_exc_S_c);
    free(f_syn_exc_S_r);
    free(f_syn_exc_S_n);
    free(f_syn_exc_S_s);
    free(f_syn_exc_S_d);

    //file_clear(p_syn_inh_size, f_syn_inh_R_c, f_syn_inh_R_r, f_syn_inh_R_n, f_syn_inh_R_s, f_syn_inh_R_d);
    //free(f_syn_inh_R_c);
    //free(f_syn_inh_R_r);
    //free(f_syn_inh_R_n);
    //free(f_syn_inh_R_s);
    //free(f_syn_inh_R_d);
    //file_clear(p_syn_inh_size, f_syn_inh_S_c, f_syn_inh_S_r, f_syn_inh_S_n, f_syn_inh_S_s, f_syn_inh_S_d);
    //free(f_syn_inh_S_c);
    //free(f_syn_inh_S_r);
    //free(f_syn_inh_S_n);
    //free(f_syn_inh_S_s);
    //free(f_syn_inh_S_d);

    arpra_ode_stepper_clear(&ode_stepper);

    mpfr_free_cache();
    return 0;
}
