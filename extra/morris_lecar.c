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

void debug (const arpra_mpfr x) {
    mpfr_out_str(stderr, 10, 80, &x, MPFR_RNDN);
    fputs("\n", stderr);
}

// General parameters
const static arpra_uint sim_steps = 100;
const static arpra_uint report_step = 20;
const static double step_size = 0.5;

const static arpra_uint grp1_size = 10;
const static arpra_uint grp2_size = 10;

const static arpra_uint reduce_step = 50;
const static double reduce_ratio = 0.3;

const static arpra_precision prec = 53;

// Neuron parameters
static arpra_range gL;    // Maximum leak conductance (mmho/cm^2)
static arpra_range VL;    // Equilibrium potential of leak conductance (mV)
static arpra_range gCa;   // Maximum Ca++ conductance (mmho/cm^2)
static arpra_range VCa;   // Equilibrium potential of Ca++ conductance (mV)
static arpra_range gK;    // Maximum K+ conductance (mmho/cm^2)
static arpra_range VK;    // Equilibrium potential of K+ conductance (mV)
static arpra_range V1;    // Potential at which M_ss(V) = 0.5 (mV)
static arpra_range V2;    // Reciprocal of voltage dependence slope of M_ss(V) (mV)
static arpra_range V3;    // Potential at which N_ss(V) = 0.5 (mV)
static arpra_range V4;    // Reciprocal of voltage dependence slope of N_ss(V) (mV)
static arpra_range phi;   // (s^-1)
static arpra_range C;     // Membrane capacitance (uF/cm^2)

// System state
// N: Fraction of open K+ channels
// V: Membrane potential (mV)
static const arpra_uint grp1_N_offset = 0;
static const arpra_uint grp2_N_offset = grp1_N_offset + grp1_size;
static const arpra_uint grp1_V_offset = grp2_N_offset + grp2_size;
static const arpra_uint grp2_V_offset = grp1_V_offset + grp1_size;
static const arpra_uint dimensions    = grp2_V_offset + grp2_size;
arpra_range *t, *dt, *x, *grp1_input, *grp2_input;

// Constants and temps
static arpra_range one, two, neg_two;
static arpra_range temp1, temp2, temp3;

static void file_init (char *grp, char *var, arpra_uint grp_size,
                       FILE **c, FILE **r, FILE **n, FILE **s, FILE **d)
{
    arpra_uint i;
    char fname[20];

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

static void file_clear (arpra_uint grp_size, FILE **c, FILE **r, FILE **n, FILE **s, FILE **d)
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

static void file_write (const arpra_range *A, arpra_uint grp_size,
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

static void dxdt (arpra_range *out,
                  const arpra_range *t, const arpra_range *x,
                  const arpra_uint x_i, const void *params)
{
    arpra_uint idx, N_offset, V_offset;
    arpra_range *input;

    // dN/dt
    if (x_i < (grp1_size + grp2_size)) {
        if (x_i < grp1_size) {
            idx = x_i;
            N_offset = grp1_N_offset;
            V_offset = grp1_V_offset;
            input = grp1_input;
        }
        else {
            idx = x_i % grp1_size;
            N_offset = grp2_N_offset;
            V_offset = grp2_V_offset;
            input = grp2_input;
        }

        // Compute K+ channel activation steady-state
        // N_ss = 1 / (1 + exp(-2 (V - V3) / V4))
        arpra_sub(&temp2, &(x[idx + V_offset]), &V3);
        arpra_mul(&temp1, &neg_two, &temp2);
        arpra_div(&temp1, &temp1, &V4);
        arpra_exp(&temp1, &temp1);
        arpra_add(&temp1, &one, &temp1);
        arpra_div(&temp1, &one, &temp1);

        // Compute tau of K+ channel activation
        // tau = 1 / (phi ((p + q) / 2))
        // p = exp(-(V - V3) / (2 V4))
        // q = exp( (V - V3) / (2 V4))
        arpra_mul(&temp3, &two, &V4);
        arpra_div(&temp3, &temp2, &temp3);
        arpra_neg(&temp2, &temp3);
        arpra_exp(&temp2, &temp2);
        arpra_exp(&temp3, &temp3);
        arpra_add(&temp2, &temp2, &temp3);
        arpra_div(&temp2, &temp2, &two);
        arpra_mul(&temp2, &phi, &temp2);
        arpra_div(&temp2, &one, &temp2);

        // Compute delta of K+ channel activation
        // dN / dt = (N_ss - N) / tau
        arpra_sub(&(out[x_i]), &temp1, &(x[idx + N_offset]));
        arpra_div(&(out[x_i]), &(out[x_i]), &temp2);
    }

    // dV/dt
    else {
        if (x_i < (grp1_size + grp2_size + grp1_size)) {
            idx = x_i % (grp1_size + grp2_size);
            N_offset = grp1_N_offset;
            V_offset = grp1_V_offset;
            input = grp1_input;
        }
        else {
            idx = x_i % (grp1_size + grp2_size + grp1_size);
            N_offset = grp2_N_offset;
            V_offset = grp2_V_offset;
            input = grp2_input;
        }

        // Compute Ca++ channel activation steady-state
        // M_ss = 1 / (1 + exp(-2 (V - V1) / V2))
        arpra_sub(&temp2, &(x[idx + V_offset]), &V1);
        arpra_mul(&temp2, &neg_two, &temp2);
        arpra_div(&temp2, &temp2, &V2);
        arpra_exp(&temp2, &temp2);
        arpra_add(&temp2, &one, &temp2);
        arpra_div(&temp2, &one, &temp2);

        // Compute leak current
        arpra_sub(&temp1, &(x[idx + V_offset]), &VL);
        arpra_mul(&temp1, &temp1, &gL);
        arpra_sub(&(out[x_i]), &(input[idx]), &temp1);

        // Compute Ca++ current
        arpra_sub(&temp1, &(x[idx + V_offset]), &VCa);
        arpra_mul(&temp1, &temp1, &gCa);
        arpra_mul(&temp1, &temp1, &temp2);
        arpra_sub(&(out[x_i]), &(out[x_i]), &temp1);

        // Compute K+ current
        arpra_sub(&temp1, &(x[idx + V_offset]), &VK);
        arpra_mul(&temp1, &temp1, &gK);
        arpra_mul(&temp1, &temp1, &(x[idx + N_offset]));
        arpra_sub(&(out[x_i]), &(out[x_i]), &temp1);

        // Compute delta of membrane potential
        // dV / dt = (I - gL (V - VL) - gCa M (V - VCa) - gK N (V - VK)) / C
        arpra_div(&(out[x_i]), &(out[x_i]), &C);
    }
}

int main (int argc, char *argv[])
{
    arpra_uint *reduce_epoch;
    clock_t run_time;
    arpra_uint i, j;

    // Allocate state memory
    t = malloc(sizeof(arpra_range));
    dt = malloc(sizeof(arpra_range));
    x = malloc(dimensions * sizeof(arpra_range));
    grp1_input = malloc(grp1_size * sizeof(arpra_range));
    grp2_input = malloc(grp2_size * sizeof(arpra_range));
    reduce_epoch = malloc(dimensions * sizeof(arpra_uint));

    // Initialise constants and variables
    for (i = 0; i < grp1_size; i++) {
        arpra_init2(&(x[grp1_N_offset + i]), prec);
        arpra_init2(&(x[grp1_V_offset + i]), prec);
        arpra_init2(&(grp1_input[i]), prec);
    }
    for (i = 0; i < grp2_size; i++) {
        arpra_init2(&(x[grp2_N_offset + i]), prec);
        arpra_init2(&(x[grp2_V_offset + i]), prec);
        arpra_init2(&(grp2_input[i]), prec);
    }
    arpra_init2(t, prec);
    arpra_init2(dt, prec);
    arpra_init2(&gL, prec);
    arpra_init2(&gCa, prec);
    arpra_init2(&gK, prec);
    arpra_init2(&VL, prec);
    arpra_init2(&VCa, prec);
    arpra_init2(&VK, prec);
    arpra_init2(&V1, prec);
    arpra_init2(&V2, prec);
    arpra_init2(&V3, prec);
    arpra_init2(&V4, prec);
    arpra_init2(&phi, prec);
    arpra_init2(&C, prec);
    arpra_init2(&one, prec);
    arpra_init2(&two, prec);
    arpra_init2(&neg_two, prec);
    arpra_init2(&temp1, prec);
    arpra_init2(&temp2, prec);
    arpra_init2(&temp3, prec);

    // Set constants and variables
    for (i = 0; i < grp1_size; i++) {
        arpra_set_d(&(x[grp1_N_offset + i]), 0.0);
        arpra_set_d(&(x[grp1_V_offset + i]), -60.0);
        arpra_set_d(&(grp1_input[i]), 80.0);
    }
    for (i = 0; i < grp2_size; i++) {
        arpra_set_d(&(x[grp2_N_offset + i]), 0.0);
        arpra_set_d(&(x[grp2_V_offset + i]), -60.0);
        arpra_set_d(&(grp2_input[i]), 80.0);
    }
    arpra_set_d(t, 0.0);
    arpra_set_d(dt, step_size);
    arpra_set_d(&gL, 2.0);
    arpra_set_d(&gCa, 4.0); // Class 1 excitability
    //arpra_set_d(&gCa, 4.4); // Class 2 excitability
    arpra_set_d(&gK, 8.0);
    arpra_set_d(&VL, -60.0);
    arpra_set_d(&VCa, 120.0);
    arpra_set_d(&VK, -80.0);
    arpra_set_d(&V1, -1.2);
    arpra_set_d(&V2, 18.0);
    arpra_set_d(&V3, 12.0); // Class 1 excitability
    //arpra_set_d(&V3, 2.0); // Class 2 excitability
    arpra_set_d(&V4, 17.4); // Class 1 excitability
    //arpra_set_d(&V4, 30.0); // Class 2 excitability
    arpra_set_d(&phi, 1.0 / 15.0); // Class 1 excitability
    //arpra_set_d(&phi, 1.0 / 25.0); // Class 2 excitability
    arpra_set_d(&C, 20.0);
    arpra_set_d(&one, 1.0);
    arpra_set_d(&two, 2.0);
    arpra_set_d(&neg_two, -2.0);

    // Initialise report files
    FILE **f_nrn1_N_c = malloc(grp1_size * sizeof(FILE *));
    FILE **f_nrn1_N_r = malloc(grp1_size * sizeof(FILE *));
    FILE **f_nrn1_N_n = malloc(grp1_size * sizeof(FILE *));
    FILE **f_nrn1_N_s = malloc(grp1_size * sizeof(FILE *));
    FILE **f_nrn1_N_d = malloc(grp1_size * sizeof(FILE *));
    file_init("nrn1", "N", grp1_size, f_nrn1_N_c, f_nrn1_N_r, f_nrn1_N_n, f_nrn1_N_s, f_nrn1_N_d);
    FILE **f_nrn1_V_c = malloc(grp1_size * sizeof(FILE *));
    FILE **f_nrn1_V_r = malloc(grp1_size * sizeof(FILE *));
    FILE **f_nrn1_V_n = malloc(grp1_size * sizeof(FILE *));
    FILE **f_nrn1_V_s = malloc(grp1_size * sizeof(FILE *));
    FILE **f_nrn1_V_d = malloc(grp1_size * sizeof(FILE *));
    file_init("nrn1", "V", grp1_size, f_nrn1_V_c, f_nrn1_V_r, f_nrn1_V_n, f_nrn1_V_s, f_nrn1_V_d);
    FILE **f_nrn2_N_c = malloc(grp2_size * sizeof(FILE *));
    FILE **f_nrn2_N_r = malloc(grp2_size * sizeof(FILE *));
    FILE **f_nrn2_N_n = malloc(grp2_size * sizeof(FILE *));
    FILE **f_nrn2_N_s = malloc(grp2_size * sizeof(FILE *));
    FILE **f_nrn2_N_d = malloc(grp2_size * sizeof(FILE *));
    file_init("nrn2", "N", grp2_size, f_nrn2_N_c, f_nrn2_N_r, f_nrn2_N_n, f_nrn2_N_s, f_nrn2_N_d);
    FILE **f_nrn2_V_c = malloc(grp2_size * sizeof(FILE *));
    FILE **f_nrn2_V_r = malloc(grp2_size * sizeof(FILE *));
    FILE **f_nrn2_V_n = malloc(grp2_size * sizeof(FILE *));
    FILE **f_nrn2_V_s = malloc(grp2_size * sizeof(FILE *));
    FILE **f_nrn2_V_d = malloc(grp2_size * sizeof(FILE *));
    file_init("nrn2", "V", grp2_size, f_nrn2_V_c, f_nrn2_V_r, f_nrn2_V_n, f_nrn2_V_s, f_nrn2_V_d);

    // ODE system
    arpra_ode_system ode_system;
    ode_system.f = dxdt;
    ode_system.t = t;
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

    for (i = 0; i < sim_steps; i++) {
        if (i % report_step == 0) printf("%u\n", i);

        for (j = 0; j < dimensions; j++) {
            reduce_epoch[j] = ode_system.x[j].nTerms;
        }

        arpra_ode_stepper_step(&ode_stepper, dt);

        for (j = 0; j < dimensions; j++) {
            arpra_reduce_last_n(&(ode_system.x[j]), (ode_system.x[j].nTerms - reduce_epoch[j]));
            if (i % reduce_step == 0) {
                arpra_reduce_small(&(ode_system.x[j]), reduce_ratio);
            }
        }

        file_write(x + grp1_N_offset, grp1_size, f_nrn1_N_c, f_nrn1_N_r, f_nrn1_N_n, f_nrn1_N_s, f_nrn1_N_d);
        file_write(x + grp1_V_offset, grp1_size, f_nrn1_V_c, f_nrn1_V_r, f_nrn1_V_n, f_nrn1_V_s, f_nrn1_V_d);
        file_write(x + grp2_N_offset, grp2_size, f_nrn2_N_c, f_nrn2_N_r, f_nrn2_N_n, f_nrn2_N_s, f_nrn2_N_d);
        file_write(x + grp2_V_offset, grp2_size, f_nrn2_V_c, f_nrn2_V_r, f_nrn2_V_n, f_nrn2_V_s, f_nrn2_V_d);
    }

    run_time = clock() - run_time;
    printf("Finished in %f seconds.\n", ((float) run_time) / CLOCKS_PER_SEC);

    // End simulation loop
    // ===================


    // Clear constants and variables
    for (i = 0; i < grp1_size; i++) {
        arpra_clear(&(x[grp1_N_offset + i]));
        arpra_clear(&(x[grp1_V_offset + i]));
        arpra_clear(&(grp1_input[i]));
    }
    for (i = 0; i < grp2_size; i++) {
        arpra_clear(&(x[grp2_N_offset + i]));
        arpra_clear(&(x[grp2_V_offset + i]));
        arpra_clear(&(grp2_input[i]));
    }
    arpra_clear(t);
    arpra_clear(dt);
    arpra_clear(&gL);
    arpra_clear(&gCa);
    arpra_clear(&gK);
    arpra_clear(&VL);
    arpra_clear(&VCa);
    arpra_clear(&VK);
    arpra_clear(&V1);
    arpra_clear(&V2);
    arpra_clear(&V3);
    arpra_clear(&V4);
    arpra_clear(&phi);
    arpra_clear(&C);
    arpra_clear(&one);
    arpra_clear(&two);
    arpra_clear(&neg_two);
    arpra_clear(&temp1);
    arpra_clear(&temp2);
    arpra_clear(&temp3);

    // Free state variables
    free(t);
    free(x);
    free(grp1_input);
    free(grp2_input);
    free(reduce_epoch);

    // Clear report files
    file_clear(grp1_size, f_nrn1_N_c, f_nrn1_N_r, f_nrn1_N_n, f_nrn1_N_s, f_nrn1_N_d);
    free(f_nrn1_N_c);
    free(f_nrn1_N_r);
    free(f_nrn1_N_n);
    free(f_nrn1_N_s);
    free(f_nrn1_N_d);
    file_clear(grp1_size, f_nrn1_V_c, f_nrn1_V_r, f_nrn1_V_n, f_nrn1_V_s, f_nrn1_V_d);
    free(f_nrn1_V_c);
    free(f_nrn1_V_r);
    free(f_nrn1_V_n);
    free(f_nrn1_V_s);
    free(f_nrn1_V_d);
    file_clear(grp2_size, f_nrn2_N_c, f_nrn2_N_r, f_nrn2_N_n, f_nrn2_N_s, f_nrn2_N_d);
    free(f_nrn2_N_c);
    free(f_nrn2_N_r);
    free(f_nrn2_N_n);
    free(f_nrn2_N_s);
    free(f_nrn2_N_d);
    file_clear(grp2_size, f_nrn2_V_c, f_nrn2_V_r, f_nrn2_V_n, f_nrn2_V_s, f_nrn2_V_d);
    free(f_nrn2_V_c);
    free(f_nrn2_V_r);
    free(f_nrn2_V_n);
    free(f_nrn2_V_s);
    free(f_nrn2_V_d);

    arpra_ode_stepper_clear(&ode_stepper);

    mpfr_free_cache();
    return 0;
}
