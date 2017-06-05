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

static mpfa_t V;    // Membrane potential (mV)
static mpfa_t M;    // Fraction of open Ca++ channels
static mpfa_t N;    // Fraction of open K+ channels
static mpfa_t I;    // Applied current (uA/cm^2)
static mpfa_t C;    // Membrane capacitance (uF/cm^2)
static mpfa_t gL;   // Maximum leak conductance (mmho/cm^2)
static mpfa_t gCa;  // Maximum Ca++ conductance (mmho/cm^2)
static mpfa_t gK;   // Maximum K+ conductance (mmho/cm^2)
static mpfa_t VL;   // Equilibrium potential of leak conductance (mV)
static mpfa_t VCa;  // Equilibrium potential of Ca++ conductance (mV)
static mpfa_t VK;   // Equilibrium potential of K+ conductance (mV)
static mpfa_t V1;   // Potential at which Mss(V) = 0.5 (mV)
static mpfa_t V2;   // Reciprocal of voltage dependence slope of Mss(V) (mV)
static mpfa_t V3;   // Potential at which Nss(V) = 0.5 (mV)
static mpfa_t V4;   // Reciprocal of voltage dependence slope of Nss(V) (mV)
static mpfa_t phi;  // Reference frequency
static mpfa_t t1, t2, t3; // Intermediate variables

void f_V (mpfa_ptr out, mpfa_srcptr V, mpfa_srcptr M, mpfa_srcptr N) {
    mpfa_sub(t1, V, VL);
    mpfa_mul(t1, t1, gL);

    mpfa_sub(t2, V, VCa);
    mpfa_mul(t2, t2, gCa);

    mpfa_sub(t3, V, VK);
    mpfa_mul(t3, t3, gK);

    mpfa_sub(out, I, t1);
    mpfa_sub(out, out, t2);
    mpfa_sub(out, out, t3);
    mpfa_div(out, out, C);
}

void f_M (mpfa_ptr out, mpfa_srcptr M, mpfa_srcptr V) {
    // compute steady-state M
    mpfa_sub(t1, V, V1);
    //mpfa_mul(t2, -2, t1);
}

void f_N (mpfa_ptr out, mpfa_srcptr N, mpfa_srcptr V) {
    // compute steady-state N
    mpfa_sub(t1, V, V3);

}

int main (int argc, char *argv[])
{
    mpfa_inits(V, M, N, I, C, gL, gCa, gK, VL, VCa, VK,
               V1, V2, V3, V4, phi, t1, t2, t3, NULL);





    /* // (nu) / (1 - nu)
       mpfr_mul_si(temp, &(zNew->u), (n - 1), MPFR_RNDU);
       mpfr_si_sub(error, 1, temp, MPFR_RNDD);
       mpfr_div(error, temp, error, MPFR_RNDU);
    */




    mpfa_clears(V, M, N, I, C, gL, gCa, gK, VL, VCa, VK,
                V1, V2, V3, V4, phi, t1, t2, t3, NULL);

    return 0;
}
