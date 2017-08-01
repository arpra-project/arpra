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

#include "mpfa.h"
#include <stdlib.h>

//#define M_DYNAMICS // Define if M steady-state is instantaneous

// Intermediate constants
static mpfa_t one, two, neg_two;


void f_V (mpfa_ptr out, mpfa_srcptr V, mpfa_srcptr M, mpfa_srcptr N,
          mpfa_srcptr gL, mpfa_srcptr gCa, mpfa_srcptr gK,
          mpfa_srcptr VL, mpfa_srcptr VCa, mpfa_srcptr VK,
          mpfa_srcptr I, mpfa_srcptr C)
{
    mpfa_t temp;
    mpfa_init(temp);

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


void f_channel (mpfa_ptr out, mpfa_srcptr A, mpfa_srcptr V,
                mpfa_srcptr Va, mpfa_srcptr Vb, mpfa_srcptr phi)
{
    mpfa_t temp1, temp2, temp3;
    mpfa_inits(temp1, temp2, temp3, NULL);

    // Compute channel activation steady-state
    // A_ss = 1 / (1 + exp(-2 (V - Va) / Vb))
    mpfa_sub(temp2, V, Va);
    mpfa_mul(temp1, neg_two, temp2);
    mpfa_div(temp1, temp1, Vb);
    mpfa_exp(temp1, temp1);
    mpfa_add(temp1, one, temp1);
    mpfa_div(temp1, one, temp1);

    // Compute tau of channel activation
    // A_tau = 1 / (phi ((p + q) / 2))
    // p = exp(-(V - Va) / (2 Vb))
    // q = exp( (V - Va) / (2 Vb))
    mpfa_mul(temp3, two, Vb);
    mpfa_div(temp3, temp2, temp3);
    mpfa_neg(temp2, temp3);
    mpfa_exp(temp2, temp2);
    mpfa_exp(temp3, temp3);
    mpfa_add(temp2, temp2, temp3);
    mpfa_div(temp2, temp2, two);
    mpfa_mul(temp2, phi, temp2);
    mpfa_div(temp2, one, temp2);

    // Compute delta of channel activation
    // dA / dt = (A_ss - A) / A_tau
    mpfa_sub(out, temp1, A);
    mpfa_div(out, out, temp2);

    mpfa_clears(temp1, temp2, temp3, NULL);
}


void M_ss (mpfa_ptr out, mpfa_srcptr V, mpfa_srcptr V1, mpfa_srcptr V2)
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


void file_init (char *var, mpfa_uint_t num, FILE **c, FILE **r, FILE **n, FILE **s, FILE **d) {
    mpfa_uint_t j;
    char fname[20];

    for (j = 0; j < num; j++) {
        sprintf(fname, "%03u_%s_c.dat", j, var);
        c[j] = fopen(fname, "w");
        sprintf(fname, "%03u_%s_r.dat", j, var);
        r[j] = fopen(fname, "w");
        sprintf(fname, "%03u_%s_n.dat", j, var);
        n[j] = fopen(fname, "w");
        sprintf(fname, "%03u_%s_s.dat", j, var);
        s[j] = fopen(fname, "w");
        sprintf(fname, "%03u_%s_d.dat", j, var);
        d[j] = fopen(fname, "w");
    }
}


void file_clear (mpfa_uint_t num, FILE **c, FILE **r, FILE **n, FILE **s, FILE **d) {
    mpfa_uint_t j;

    for (j = 0; j < num; j++) {
        fclose(c[j]);
        fclose(r[j]);
        fclose(n[j]);
        fclose(s[j]);
        fclose(d[j]);
    }
}


void file_write (mpfa_srcptr a, mpfa_uint_t num, FILE **c, FILE **r, FILE **n, FILE **s, FILE **d) {
    mpfa_uint_t j, k;

    for (j = 0; j < num; j++) {
        mpfr_out_str(c[j], 10, 80, &(a[j].centre), MPFR_RNDN);
        fputc('\n', c[j]);
        mpfr_out_str(r[j], 10, 80, &(a[j].radius), MPFR_RNDN);
        fputc('\n', r[j]);
        fprintf(n[j], "%u\n", a[j].nTerms);
        for (k = 0; k < a[j].nTerms; k++) {
            fprintf(s[j], "%u ", a[j].symbols[k]);
            mpfr_out_str(d[j], 10, 80, &(a[j].deviations[k]), MPFR_RNDN);
            fputc(' ', d[j]);
        }
        fputc('\n', s[j]);
        fputc('\n', d[j]);
    }
}


int main (int argc, char *argv[])
{
    const mpfa_uint_t sim_steps = 500;
    const mpfa_uint_t n_grp1 = 10;
    const mpfa_uint_t n_grp2 = 10;
    const mpfa_uint_t n_grp1_grp2 = n_grp1 * n_grp2;
    const mpfa_uint_t condense_step = 50;
    const double condense_ratio = 0.3;

    size_t size;
    mpfa_uint_t i, j;
    mpfa_uint_t M_term, N_term, V_term;

    // Init state variables
    size = n_grp1 * sizeof(mpfa_t);
    mpfa_ptr M = malloc(size);  // Fraction of open Ca++ channels
    mpfa_ptr N = malloc(size);  // Fraction of open K+ channels
    mpfa_ptr V = malloc(size);  // Membrane potential (mV)
    mpfa_ptr I = malloc(size);  // Applied current (uA/cm^2)

    for (j = 0; j < n_grp1; j++) {
        mpfa_init(&(M[j]));
        mpfa_init(&(N[j]));
        mpfa_init(&(V[j]));
        mpfa_init(&(I[j]));

        mpfa_set_d(&(M[j]), 0.0);
        mpfa_set_d(&(N[j]), 0.0);
        mpfa_set_d(&(V[j]), -60.0);
        mpfa_set_d(&(I[j]), 80.0);
    }

    // Init file handles
    FILE **f_M_c = malloc(n_grp1 * sizeof(FILE *)); // centre
    FILE **f_M_r = malloc(n_grp1 * sizeof(FILE *)); // radius
    FILE **f_M_n = malloc(n_grp1 * sizeof(FILE *)); // nTerms
    FILE **f_M_s = malloc(n_grp1 * sizeof(FILE *)); // symbols
    FILE **f_M_d = malloc(n_grp1 * sizeof(FILE *)); // deviations
    file_init("m", n_grp1, f_M_c, f_M_r, f_M_n, f_M_s, f_M_d);

    FILE **f_N_c = malloc(n_grp1 * sizeof(FILE *)); // centre
    FILE **f_N_r = malloc(n_grp1 * sizeof(FILE *)); // radius
    FILE **f_N_n = malloc(n_grp1 * sizeof(FILE *)); // nterms
    FILE **f_N_s = malloc(n_grp1 * sizeof(FILE *)); // symbols
    FILE **f_N_d = malloc(n_grp1 * sizeof(FILE *)); // deviations
    file_init("n", n_grp1, f_N_c, f_N_r, f_N_n, f_N_s, f_N_d);

    FILE **f_V_c = malloc(n_grp1 * sizeof(FILE *)); // centre
    FILE **f_V_r = malloc(n_grp1 * sizeof(FILE *)); // radius
    FILE **f_V_n = malloc(n_grp1 * sizeof(FILE *)); // nterms
    FILE **f_V_s = malloc(n_grp1 * sizeof(FILE *)); // symbols
    FILE **f_V_d = malloc(n_grp1 * sizeof(FILE *)); // deviations
    file_init("v", n_grp1, f_V_c, f_V_r, f_V_n, f_V_s, f_V_d);

    // Init parameters
    mpfa_t dM;    // Delta M
    mpfa_t dN;    // Delta N
    mpfa_t dV;    // Delta V
    mpfa_t dt;    // Delta time
    mpfa_t gL;    // Maximum leak conductance (mmho/cm^2)
    mpfa_t VL;    // Equilibrium potential of leak conductance (mV)
    mpfa_t gCa;   // Maximum Ca++ conductance (mmho/cm^2)
    mpfa_t VCa;   // Equilibrium potential of Ca++ conductance (mV)
    mpfa_t gK;    // Maximum K+ conductance (mmho/cm^2)
    mpfa_t VK;    // Equilibrium potential of K+ conductance (mV)
    mpfa_t V1;    // Potential at which Mss(V) = 0.5 (mV)
    mpfa_t V2;    // Reciprocal of voltage dependence slope of Mss(V) (mV)
    mpfa_t V3;    // Potential at which Nss(V) = 0.5 (mV)
    mpfa_t V4;    // Reciprocal of voltage dependence slope of Nss(V) (mV)
    mpfa_t M_phi; // (s^-1)
    mpfa_t N_phi; // (s^-1)
    mpfa_t C;     // Membrane capacitance (uF/cm^2)

    mpfa_inits(dM, dN, dV, dt,
               gL, gCa, gK,
               VL, VCa, VK,
               V1, V2, V3, V4,
               M_phi, N_phi, C,
               one, two, neg_two,
               NULL);

    mpfa_set_d(dt, 1.0);
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
    mpfa_set_d(M_phi, 1.0);
    mpfa_set_d(N_phi, 1.0 / 15.0); // Class 1 excitability
    //mpfa_set_d(N_phi, 1.0 / 25.0); // Class 2 excitability
    mpfa_set_d(C, 20.0);
    mpfa_set_d(one, 1.0);
    mpfa_set_d(two, 2.0);
    mpfa_set_d(neg_two, -2.0);


    // Simulation loop
    // ===============

    for (i = 0; i < sim_steps; i++) {
        printf("%u\n", i);

        for (j = 0; j < n_grp1; j++) {

            /* // (nu) / (1 - nu)
               mpfr_mul_si(temp, u, (n - 1), MPFR_RNDU);
               mpfr_si_sub(error, 1, temp, MPFR_RNDD);
               mpfr_div(error, temp, error, MPFR_RNDU);
            */

            M_term = M[j].nTerms;
            N_term = N[j].nTerms;
            V_term = V[j].nTerms;

#ifdef M_DYNAMICS // If we need M dynamics
            f_channel(dM, &(M[j]), &(V[j]), V1, V2, M_phi);
#else // Else M steady-state is instantaneous
            M_ss(&(M[j]), &(V[j]), V1, V2);
#endif // M_DYNAMICS
            f_channel(dN, &(N[j]), &(V[j]), V3, V4, N_phi);
            f_V(dV, &(V[j]), &(M[j]), &(N[j]), gL, gCa, gK, VL, VCa, VK, I, C);

#ifdef M_DYNAMICS // If we need M dynamics
            mpfa_mul(dM, dM, dt);
            mpfa_add(&(M[j]), &(M[j]), dM);
#endif // M_DYNAMICS
            mpfa_mul(dN, dN, dt);
            mpfa_add(&(N[j]), &(N[j]), dN);
            mpfa_mul(dV, dV, dt);
            mpfa_add(&(V[j]), &(V[j]), dV);

            mpfa_condense_last_n(&(M[j]), (M[j].nTerms - M_term));
            mpfa_condense_last_n(&(N[j]), (N[j].nTerms - N_term));
            mpfa_condense_last_n(&(V[j]), (V[j].nTerms - V_term));

            if (i % condense_step == 0) {
                mpfa_condense_small(&(M[j]), condense_ratio);
                mpfa_condense_small(&(N[j]), condense_ratio);
                mpfa_condense_small(&(V[j]), condense_ratio);
            }
        }

        file_write(M, n_grp1, f_M_c, f_M_r, f_M_n, f_M_s, f_M_d);
        file_write(N, n_grp1, f_N_c, f_N_r, f_N_n, f_N_s, f_N_d);
        file_write(V, n_grp1, f_V_c, f_V_r, f_V_n, f_V_s, f_V_d);
    }

    // End simulation loop
    // ===================


    // Clear state variables
    for (j = 0; j < n_grp1; j++) {
        mpfa_clear(&(M[j]));
        mpfa_clear(&(N[j]));
        mpfa_clear(&(V[j]));
        mpfa_clear(&(I[j]));
    }

    free(V);
    free(M);
    free(N);
    free(I);

    // Clear parameters
    mpfa_clears(dM, dN, dV, dt,
                gL, gCa, gK,
                VL, VCa, VK,
                V1, V2, V3, V4,
                M_phi, N_phi, C,
                one, two, neg_two,
                NULL);

    // Clear file handles
    file_clear(n_grp1, f_M_c, f_M_r, f_M_n, f_M_s, f_M_d);
    free(f_M_c);
    free(f_M_r);
    free(f_M_n);
    free(f_M_s);
    free(f_M_d);

    file_clear(n_grp1, f_N_c, f_N_r, f_N_n, f_N_s, f_N_d);
    free(f_N_c);
    free(f_N_r);
    free(f_N_n);
    free(f_N_s);
    free(f_N_d);

    file_clear(n_grp1, f_V_c, f_V_r, f_V_n, f_V_s, f_V_d);
    free(f_V_c);
    free(f_V_r);
    free(f_V_n);
    free(f_V_s);
    free(f_V_d);

    mpfr_free_cache();
    return 0;
}
