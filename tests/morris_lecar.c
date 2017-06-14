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


void f_A (mpfa_ptr out, mpfa_srcptr A, mpfa_srcptr V,
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


int main (int argc, char *argv[])
{
    unsigned i;
    const unsigned sim_time = 100;

    mpfa_t V;    // Membrane potential (mV)
    mpfa_t M;    // Fraction of open Ca++ channels
    mpfa_t N;    // Fraction of open K+ channels
    mpfa_t I;    // Applied current (uA/cm^2)
    mpfa_t C;    // Membrane capacitance (uF/cm^2)

    mpfa_t dt;   // Delta time
    mpfa_t dV;   // Delta V
    mpfa_t dM;   // Delta M
    mpfa_t dN;   // Delta N

    mpfa_t gL;   // Maximum leak conductance (mmho/cm^2)
    mpfa_t gCa;  // Maximum Ca++ conductance (mmho/cm^2)
    mpfa_t gK;   // Maximum K+ conductance (mmho/cm^2)

    mpfa_t VL;   // Equilibrium potential of leak conductance (mV)
    mpfa_t VCa;  // Equilibrium potential of Ca++ conductance (mV)
    mpfa_t VK;   // Equilibrium potential of K+ conductance (mV)

    mpfa_t V1;   // Potential at which Mss(V) = 0.5 (mV)
    mpfa_t V2;   // Reciprocal of voltage dependence slope of Mss(V) (mV)
    mpfa_t V3;   // Potential at which Nss(V) = 0.5 (mV)
    mpfa_t V4;   // Reciprocal of voltage dependence slope of Nss(V) (mV)

    mpfa_t M_phi; // (s^-1)
    mpfa_t N_phi; // (s^-1)

    mpfa_inits(V, M, N, I, C,
               dt, dV, dM, dN,
               gL, gCa, gK,
               VL, VCa, VK,
               V1, V2, V3, V4,
               M_phi, N_phi,
               one, two, neg_two,
               NULL);

    // Initialise Variables and parameters
    mpfa_set_d(V, -60.0);
    mpfa_set_d(M, 0.0);
    mpfa_set_d(N, 0.0);
    mpfa_set_d(I, 80.0);
    mpfa_set_d(C, 20.0);

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

    mpfa_set_d(dt, 1.0);

    // Initialise constants
    mpfa_set_d(one, 1.0);
    mpfa_set_d(two, 2.0);
    mpfa_set_d(neg_two, -2.0);

    for (i = 0; i < sim_time; i++) {
        /* // (nu) / (1 - nu)
           mpfr_mul_si(temp, u, (n - 1), MPFR_RNDU);
           mpfr_si_sub(error, 1, temp, MPFR_RNDD);
           mpfr_div(error, temp, error, MPFR_RNDU);
        */

#ifdef M_DYNAMICS // If we need M dynamics
        f_A(dM, M, V, V1, V2, M_phi);
#else // Else M steady-state is instantaneous
        M_ss(M, V, V1, V2);
#endif
        f_A(dN, N, V, V3, V4, N_phi);
        f_V(dV, V, M, N, gL, gCa, gK, VL, VCa, VK, I, C);

#ifdef M_DYNAMICS // If we need M dynamics
        mpfa_mul(dM, dM, dt);
        mpfa_add(M, M, dM);
        //printf("Mc: "); mpfr_out_str (stdout, 10, 100, &(M->centre), MPFR_RNDN); putchar('\n');
        //printf("Mr: "); mpfr_out_str (stdout, 10, 100, &(M->radius), MPFR_RNDN); putchar('\n');
#endif

        mpfa_mul(dN, dN, dt);
        mpfa_add(N, N, dN);
        //printf("Nc: "); mpfr_out_str (stdout, 10, 100, &(N->centre), MPFR_RNDN); putchar('\n');
        //printf("Nr: "); mpfr_out_str (stdout, 10, 100, &(N->radius), MPFR_RNDN); putchar('\n');

        mpfa_mul(dV, dV, dt);
        mpfa_add(V, V, dV);
        printf("Vc: "); mpfr_out_str (stdout, 10, 100, &(V->centre), MPFR_RNDN); putchar('\n');
        //printf("Vr: "); mpfr_out_str (stdout, 10, 100, &(V->radius), MPFR_RNDN); putchar('\n');
    }

    mpfa_clears(V, M, N, I, C,
                dt, dV, dM, dN,
                gL, gCa, gK,
                VL, VCa, VK,
                V1, V2, V3, V4,
                M_phi, N_phi,
                one, two, neg_two,
                NULL);

    return 0;
}
