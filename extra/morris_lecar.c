/*
 * morris_lecar.c -- Test Morris-Lecar model.
 *
 * Copyright 2016-2017 James Paul Turner.
 *
 * This file is part of the MPFA library.
 *
 * The MPFA library is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * The MPFA library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with the MPFA library. If not, see <http://www.gnu.org/licenses/>.
 */

#include <stdlib.h>
#include <stdio.h>
#include <mpfa.h>

// General parameters
const static mpfa_uint_t sim_steps = 500;
const static double step_size = 1.0;

const static mpfa_uint_t n_grp1 = 10;
const static mpfa_uint_t n_grp2 = 10;

const static mpfa_uint_t reduce_step = 50;
const static double reduce_ratio = 0.3;

const static mpfa_prec_t prec = 53;

// Neuron parameters
static mpfa_t gL;    // Maximum leak conductance (mmho/cm^2)
static mpfa_t VL;    // Equilibrium potential of leak conductance (mV)
static mpfa_t gCa;   // Maximum Ca++ conductance (mmho/cm^2)
static mpfa_t VCa;   // Equilibrium potential of Ca++ conductance (mV)
static mpfa_t gK;    // Maximum K+ conductance (mmho/cm^2)
static mpfa_t VK;    // Equilibrium potential of K+ conductance (mV)
static mpfa_t V1;    // Potential at which M_ss(V) = 0.5 (mV)
static mpfa_t V2;    // Reciprocal of voltage dependence slope of M_ss(V) (mV)
static mpfa_t V3;    // Potential at which N_ss(V) = 0.5 (mV)
static mpfa_t V4;    // Reciprocal of voltage dependence slope of N_ss(V) (mV)
static mpfa_t phi;   // (s^-1)
static mpfa_t C;     // Membrane capacitance (uF/cm^2)

// Constants
static mpfa_t one, two, neg_two;


void M_ss (mpfa_ptr out, mpfa_srcptr V)
{
    // Compute Ca++ channel activation steady-state
    // M_ss = 1 / (1 + exp(-2 (V - V1) / V2))
    mpfa_sub(out, V, V1);
    mpfa_mul(out, neg_two, out);
    mpfa_div(out, out, V2);
    mpfa_exp(out, out);
    mpfa_add(out, one, out);
    mpfa_div(out, one, out);
}


void d_V (mpfa_ptr out, mpfa_srcptr V, mpfa_srcptr M, mpfa_srcptr N, mpfa_srcptr I)
{
    mpfa_t temp;
    mpfa_init2(temp, prec);

    // Compute leak current
    mpfa_sub(temp, V, VL);
    mpfa_mul(temp, temp, gL);
    mpfa_sub(out, I, temp);

    // Compute Ca++ current
    mpfa_sub(temp, V, VCa);
    mpfa_mul(temp, temp, gCa);
    mpfa_mul(temp, temp, M);
    mpfa_sub(out, out, temp);

    // Compute K+ current
    mpfa_sub(temp, V, VK);
    mpfa_mul(temp, temp, gK);
    mpfa_mul(temp, temp, N);
    mpfa_sub(out, out, temp);

    // Compute delta of membrane potential
    // dV / dt = (I - gL (V - VL) - gCa M (V - VCa) - gK N (V - VK)) / C
    mpfa_div(out, out, C);

    mpfa_clear(temp);
}


void d_N (mpfa_ptr out, mpfa_srcptr N, mpfa_srcptr V)
{
    mpfa_t temp1, temp2, temp3;
    mpfa_inits2(prec, temp1, temp2, temp3, NULL);

    // Compute K+ channel activation steady-state
    // N_ss = 1 / (1 + exp(-2 (V - V3) / V4))
    mpfa_sub(temp2, V, V3);
    mpfa_mul(temp1, neg_two, temp2);
    mpfa_div(temp1, temp1, V4);
    mpfa_exp(temp1, temp1);
    mpfa_add(temp1, one, temp1);
    mpfa_div(temp1, one, temp1);

    // Compute tau of K+ channel activation
    // tau = 1 / (phi ((p + q) / 2))
    // p = exp(-(V - V3) / (2 V4))
    // q = exp( (V - V3) / (2 V4))
    mpfa_mul(temp3, two, V4);
    mpfa_div(temp3, temp2, temp3);
    mpfa_neg(temp2, temp3);
    mpfa_exp(temp2, temp2);
    mpfa_exp(temp3, temp3);
    mpfa_add(temp2, temp2, temp3);
    mpfa_div(temp2, temp2, two);
    mpfa_mul(temp2, phi, temp2);
    mpfa_div(temp2, one, temp2);

    // Compute delta of K+ channel activation
    // dN / dt = (N_ss - N) / tau
    mpfa_sub(out, temp1, N);
    mpfa_div(out, out, temp2);

    mpfa_clears(temp1, temp2, temp3, NULL);
}


void file_init (char *grp, char *var, mpfa_uint_t num,
                FILE **c, FILE **r, FILE **n, FILE **s, FILE **d)
{
    mpfa_uint_t j;
    char fname[20];

    for (j = 0; j < num; j++) {
        sprintf(fname, "%s_%03u_%s_c.dat", grp, (unsigned) j, var);
        c[j] = fopen(fname, "w");
        sprintf(fname, "%s_%03u_%s_r.dat", grp, (unsigned) j, var);
        r[j] = fopen(fname, "w");
        sprintf(fname, "%s_%03u_%s_n.dat", grp, (unsigned) j, var);
        n[j] = fopen(fname, "w");
        sprintf(fname, "%s_%03u_%s_s.dat", grp, (unsigned) j, var);
        s[j] = fopen(fname, "w");
        sprintf(fname, "%s_%03u_%s_d.dat", grp, (unsigned) j, var);
        d[j] = fopen(fname, "w");
    }
}


void file_clear (mpfa_uint_t num, FILE **c, FILE **r, FILE **n, FILE **s, FILE **d)
{
    mpfa_uint_t j;

    for (j = 0; j < num; j++) {
        fclose(c[j]);
        fclose(r[j]);
        fclose(n[j]);
        fclose(s[j]);
        fclose(d[j]);
    }
}


void file_write (mpfa_srcptr A, mpfa_uint_t i,
                 FILE **c, FILE **r, FILE **n, FILE **s, FILE **d)
{
    mpfa_uint_t j;

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


int main (int argc, char *argv[])
{
    unsigned int i, j, N_mark, V_mark;
    mpfa_t dN, dV, dt, M;

    // Init parameters
    mpfa_inits2(prec,
                dN, dV, dt, M,
                gL, gCa, gK,
                VL, VCa, VK,
                V1, V2, V3, V4,
                phi, C,
                one, two, neg_two,
                NULL);

    mpfa_set_d(dt, step_size);

    mpfa_set_d(gL, 2.0);
    mpfa_set_d(gCa, 4.0); // Class 1 excitability
    //mpfa_set_d(gCa, 4.4); // Class 2 excitability
    mpfa_set_d(gK, 8.0);
    mpfa_set_d(VL, -60.0);
    mpfa_set_d(VCa, 120.0);
    mpfa_set_d(VK, -80.0);
    mpfa_set_d(V1, -1.2);
    mpfa_set_d(V2, 18.0);
    mpfa_set_d(V3, 12.0); // Class 1 excitability
    //mpfa_set_d(V3, 2.0); // Class 2 excitability
    mpfa_set_d(V4, 17.4); // Class 1 excitability
    //mpfa_set_d(V4, 30.0); // Class 2 excitability
    mpfa_set_d(phi, 1.0 / 15.0); // Class 1 excitability
    //mpfa_set_d(phi, 1.0 / 25.0); // Class 2 excitability
    mpfa_set_d(C, 20.0);

    mpfa_set_d(one, 1.0);
    mpfa_set_d(two, 2.0);
    mpfa_set_d(neg_two, -2.0);

    // Init neuron group 1 state variables
    mpfa_ptr nrn1_N = malloc(n_grp1 * sizeof(mpfa_t));  // Fraction of open K+ channels
    mpfa_ptr nrn1_V = malloc(n_grp1 * sizeof(mpfa_t));  // Membrane potential (mV)
    mpfa_ptr nrn1_I = malloc(n_grp1 * sizeof(mpfa_t));  // Applied current (uA/cm^2)

    for (j = 0; j < n_grp1; j++) {
        mpfa_init2(&(nrn1_N[j]), prec);
        mpfa_init2(&(nrn1_V[j]), prec);
        mpfa_init2(&(nrn1_I[j]), prec);
        mpfa_set_d(&(nrn1_N[j]), 0.0);
        mpfa_set_d(&(nrn1_V[j]), -60.0);
        mpfa_set_d(&(nrn1_I[j]), 80.0);
    }

    // Init neuron group 1 files
    FILE **f_nrn1_N_c = malloc(n_grp1 * sizeof(FILE *));
    FILE **f_nrn1_N_r = malloc(n_grp1 * sizeof(FILE *));
    FILE **f_nrn1_N_n = malloc(n_grp1 * sizeof(FILE *));
    FILE **f_nrn1_N_s = malloc(n_grp1 * sizeof(FILE *));
    FILE **f_nrn1_N_d = malloc(n_grp1 * sizeof(FILE *));
    file_init("nrn1", "N", n_grp1, f_nrn1_N_c, f_nrn1_N_r, f_nrn1_N_n, f_nrn1_N_s, f_nrn1_N_d);

    FILE **f_nrn1_V_c = malloc(n_grp1 * sizeof(FILE *));
    FILE **f_nrn1_V_r = malloc(n_grp1 * sizeof(FILE *));
    FILE **f_nrn1_V_n = malloc(n_grp1 * sizeof(FILE *));
    FILE **f_nrn1_V_s = malloc(n_grp1 * sizeof(FILE *));
    FILE **f_nrn1_V_d = malloc(n_grp1 * sizeof(FILE *));
    file_init("nrn1", "V", n_grp1, f_nrn1_V_c, f_nrn1_V_r, f_nrn1_V_n, f_nrn1_V_s, f_nrn1_V_d);

    // Init neuron group 2 state variables
    mpfa_ptr nrn2_N = malloc(n_grp2 * sizeof(mpfa_t));  // Fraction of open K+ channels
    mpfa_ptr nrn2_V = malloc(n_grp2 * sizeof(mpfa_t));  // Membrane potential (mV)
    mpfa_ptr nrn2_I = malloc(n_grp2 * sizeof(mpfa_t));  // Applied current (uA/cm^2)

    for (j = 0; j < n_grp2; j++) {
        mpfa_init2(&(nrn2_N[j]), prec);
        mpfa_init2(&(nrn2_V[j]), prec);
        mpfa_init2(&(nrn2_I[j]), prec);
        mpfa_set_d(&(nrn2_N[j]), 0.0);
        mpfa_set_d(&(nrn2_V[j]), -60.0);
        mpfa_set_d(&(nrn2_I[j]), 80.0);
    }

    // Init neuron group 2 files
    FILE **f_nrn2_N_c = malloc(n_grp2 * sizeof(FILE *));
    FILE **f_nrn2_N_r = malloc(n_grp2 * sizeof(FILE *));
    FILE **f_nrn2_N_n = malloc(n_grp2 * sizeof(FILE *));
    FILE **f_nrn2_N_s = malloc(n_grp2 * sizeof(FILE *));
    FILE **f_nrn2_N_d = malloc(n_grp2 * sizeof(FILE *));
    file_init("nrn2", "N", n_grp2, f_nrn2_N_c, f_nrn2_N_r, f_nrn2_N_n, f_nrn2_N_s, f_nrn2_N_d);

    FILE **f_nrn2_V_c = malloc(n_grp2 * sizeof(FILE *));
    FILE **f_nrn2_V_r = malloc(n_grp2 * sizeof(FILE *));
    FILE **f_nrn2_V_n = malloc(n_grp2 * sizeof(FILE *));
    FILE **f_nrn2_V_s = malloc(n_grp2 * sizeof(FILE *));
    FILE **f_nrn2_V_d = malloc(n_grp2 * sizeof(FILE *));
    file_init("nrn2", "V", n_grp2, f_nrn2_V_c, f_nrn2_V_r, f_nrn2_V_n, f_nrn2_V_s, f_nrn2_V_d);


    // Begin simulation loop
    // =====================

    for (i = 0; i < sim_steps; i++) {
        printf("%u\n", i);


        // Neuron group 1
        // ==============

        for (j = 0; j < n_grp1; j++) {

            /* // (nu) / (1 - nu)
               mpfr_mul_si(temp, u, (n - 1), MPFR_RNDU);
               mpfr_si_sub(error, 1, temp, MPFR_RNDD);
               mpfr_div(error, temp, error, MPFR_RNDU);
            */

            M_ss(M, &(nrn1_V[j]));

            N_mark = nrn1_N[j].nTerms;
            V_mark = nrn1_V[j].nTerms;

            d_N(dN, &(nrn1_N[j]), &(nrn1_V[j]));
            d_V(dV, &(nrn1_V[j]), M, &(nrn1_N[j]), &(nrn1_I[j]));

            mpfa_mul(dN, dN, dt);
            mpfa_add(&(nrn1_N[j]), &(nrn1_N[j]), dN);
            mpfa_mul(dV, dV, dt);
            mpfa_add(&(nrn1_V[j]), &(nrn1_V[j]), dV);

            mpfa_reduce_last_n(&(nrn1_N[j]), (nrn1_N[j].nTerms - N_mark));
            mpfa_reduce_last_n(&(nrn1_V[j]), (nrn1_V[j].nTerms - V_mark));

            if (i % reduce_step == 0) {
                mpfa_reduce_small(&(nrn1_N[j]), reduce_ratio);
                mpfa_reduce_small(&(nrn1_V[j]), reduce_ratio);
            }

            file_write(nrn1_N, j, f_nrn1_N_c, f_nrn1_N_r, f_nrn1_N_n, f_nrn1_N_s, f_nrn1_N_d);
            file_write(nrn1_V, j, f_nrn1_V_c, f_nrn1_V_r, f_nrn1_V_n, f_nrn1_V_s, f_nrn1_V_d);
        }


        // Neuron group 2
        // ==============

        for (j = 0; j < n_grp2; j++) {

            /* // (nu) / (1 - nu)
               mpfr_mul_si(temp, u, (n - 1), MPFR_RNDU);
               mpfr_si_sub(error, 1, temp, MPFR_RNDD);
               mpfr_div(error, temp, error, MPFR_RNDU);
            */

            M_ss(M, &(nrn2_V[j]));

            N_mark = nrn2_N[j].nTerms;
            V_mark = nrn2_V[j].nTerms;

            d_N(dN, &(nrn2_N[j]), &(nrn2_V[j]));
            d_V(dV, &(nrn2_V[j]), M, &(nrn2_N[j]), &(nrn2_I[j]));

            mpfa_mul(dN, dN, dt);
            mpfa_add(&(nrn2_N[j]), &(nrn2_N[j]), dN);
            mpfa_mul(dV, dV, dt);
            mpfa_add(&(nrn2_V[j]), &(nrn2_V[j]), dV);

            mpfa_reduce_last_n(&(nrn2_N[j]), (nrn2_N[j].nTerms - N_mark));
            mpfa_reduce_last_n(&(nrn2_V[j]), (nrn2_V[j].nTerms - V_mark));

            if (i % reduce_step == 0) {
                mpfa_reduce_small(&(nrn2_N[j]), reduce_ratio);
                mpfa_reduce_small(&(nrn2_V[j]), reduce_ratio);
            }

            file_write(nrn2_N, j, f_nrn2_N_c, f_nrn2_N_r, f_nrn2_N_n, f_nrn2_N_s, f_nrn2_N_d);
            file_write(nrn2_V, j, f_nrn2_V_c, f_nrn2_V_r, f_nrn2_V_n, f_nrn2_V_s, f_nrn2_V_d);
        }
    }

    // End simulation loop
    // ===================


    // Clear neuron group 1 state variables
    for (j = 0; j < n_grp1; j++) {
        mpfa_clear(&(nrn1_N[j]));
        mpfa_clear(&(nrn1_V[j]));
        mpfa_clear(&(nrn1_I[j]));
    }
    free(nrn1_N);
    free(nrn1_V);
    free(nrn1_I);

    // Clear neuron group 1 files
    file_clear(n_grp1, f_nrn1_N_c, f_nrn1_N_r, f_nrn1_N_n, f_nrn1_N_s, f_nrn1_N_d);
    free(f_nrn1_N_c);
    free(f_nrn1_N_r);
    free(f_nrn1_N_n);
    free(f_nrn1_N_s);
    free(f_nrn1_N_d);

    file_clear(n_grp1, f_nrn1_V_c, f_nrn1_V_r, f_nrn1_V_n, f_nrn1_V_s, f_nrn1_V_d);
    free(f_nrn1_V_c);
    free(f_nrn1_V_r);
    free(f_nrn1_V_n);
    free(f_nrn1_V_s);
    free(f_nrn1_V_d);

    // Clear neuron group 2 state variables
    for (j = 0; j < n_grp2; j++) {
        mpfa_clear(&(nrn2_N[j]));
        mpfa_clear(&(nrn2_V[j]));
        mpfa_clear(&(nrn2_I[j]));
    }
    free(nrn2_N);
    free(nrn2_V);
    free(nrn2_I);

    // Clear neuron group 2 files
    file_clear(n_grp2, f_nrn2_N_c, f_nrn2_N_r, f_nrn2_N_n, f_nrn2_N_s, f_nrn2_N_d);
    free(f_nrn2_N_c);
    free(f_nrn2_N_r);
    free(f_nrn2_N_n);
    free(f_nrn2_N_s);
    free(f_nrn2_N_d);

    file_clear(n_grp2, f_nrn2_V_c, f_nrn2_V_r, f_nrn2_V_n, f_nrn2_V_s, f_nrn2_V_d);
    free(f_nrn2_V_c);
    free(f_nrn2_V_r);
    free(f_nrn2_V_n);
    free(f_nrn2_V_s);
    free(f_nrn2_V_d);

    // Clear parameters
    mpfa_clears(dN, dV, dt, M,
                gL, gCa, gK,
                VL, VCa, VK,
                V1, V2, V3, V4,
                phi, C,
                one, two, neg_two,
                NULL);

    mpfr_free_cache();
    return 0;
}
