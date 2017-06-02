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

int main (int argc, char *argv[]) {
    mpfa_t V;   // Membrane potential (mV)
    mpfa_t M;   // Fraction of open Ca++ channels
    mpfa_t N;   // Fraction of open K+ channels
    mpfa_t I;   // Applied current (uA/cm^2)
    mpfa_t C;   // Membrane capacitance (uF/cm^2)
    mpfa_t gL;  // Maximum leak conductance (mmho/cm^2)
    mpfa_t gCa; // Maximum Ca++ conductance (mmho/cm^2)
    mpfa_t gK;  // Maximum K+ conductance (mmho/cm^2)
    mpfa_t VL;  // Equilibrium potential of leak conductance (mV)
    mpfa_t VCa; // Equilibrium potential of Ca++ conductance (mV)
    mpfa_t VK;  // Equilibrium potential of K+ conductance (mV)
    mpfa_t V1;  // Potential at which Mss(V) = 0.5 (mV)
    mpfa_t V2;  // Reciprocal of slope of voltage dependence of Mss(V) (mV)
    mpfa_t V3;  // Potential at which Nss(V) = 0.5 (mV)
    mpfa_t V4;  // Reciprocal of slope of voltage dependence of Nss(V) (mV)
    mpfa_t phi; // Reference frequency

    mpfa_inits(V, M, N, I, C, gL, gCa, gK, VL, VCa, VK, V1, V2, V3, V4, phi);


    /* // (nu) / (1 - nu)
       mpfr_mul_si(temp, &(zNew->u), (n - 1), MPFR_RNDU);
       mpfr_si_sub(error, 1, temp, MPFR_RNDD);
       mpfr_div(error, temp, error, MPFR_RNDU);
    */




    mpfa_clears(V, M, N, I, C, gL, gCa, gK, VL, VCa, VK, V1, V2, V3, V4, phi);
}
