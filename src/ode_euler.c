/*
 * ode_euler.c -- Explicit Euler ODE stepper.
 *
 * Copyright 2018 James Paul Turner.
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

#include "arpra-impl.h"

#define euler_stages 1

typedef struct euler_scratch_struct
{
    arpra_range *k_0;
    arpra_range *x_new;
    arpra_range temp_x;
} euler_scratch;

static void euler_init (arpra_ode_stepper *stepper, arpra_ode_system *system)
{
    arpra_uint x_idx;
    arpra_precision prec_x, prec_internal;
    euler_scratch *scratch;

    // Allocate scratch memory.
    scratch = malloc(sizeof(euler_scratch));
    scratch->k_0 = malloc(system->dims * sizeof(arpra_range));
    scratch->x_new = malloc(system->dims * sizeof(arpra_range));

    // Initialise scratch memory.
    prec_internal = arpra_get_internal_precision();
    for (x_idx = 0; x_idx < system->dims; x_idx++) {
        prec_x = arpra_get_precision(&(system->x[x_idx]));
        arpra_init2(&(scratch->k_0[x_idx]), prec_x);
        arpra_init2(&(scratch->x_new[x_idx]), prec_x);
    }
    arpra_init2(&(scratch->temp_x), prec_internal);

    // Set stepper parameters.
    stepper->method = arpra_ode_euler;
    stepper->system = system;
    stepper->error = NULL;
    stepper->scratch = scratch;
}

static void euler_clear (arpra_ode_stepper *stepper)
{
    arpra_uint x_idx;
    arpra_ode_system *system;
    euler_scratch *scratch;

    system = stepper->system;
    scratch = (euler_scratch *) stepper->scratch;

    // Clear scratch memory.
    for (x_idx = 0; x_idx < system->dims; x_idx++) {
        arpra_clear(&(scratch->k_0[x_idx]));
        arpra_clear(&(scratch->x_new[x_idx]));
    }
    arpra_clear(&(scratch->temp_x));

    // Free scratch memory.
    free(scratch->k_0);
    free(scratch->x_new);
    free(scratch);
}

static void euler_step (arpra_ode_stepper *stepper, const arpra_range *h)
{
    arpra_uint x_idx;
    arpra_precision prec_x;
    arpra_ode_system *system;
    euler_scratch *scratch;

    system = stepper->system;
    scratch = (euler_scratch *) stepper->scratch;

    // Synchronise scratch precision and prepare step parameters.
    for (x_idx = 0; x_idx < system->dims; x_idx++) {
        prec_x = arpra_get_precision(&(system->x[x_idx]));
        arpra_set_precision(&(scratch->k_0[x_idx]), prec_x);
        arpra_set_precision(&(scratch->x_new[x_idx]), prec_x);
    }

    // k[0] = f(t, x(t))
    for (x_idx = 0; x_idx < system->dims; x_idx++) {
        system->f(scratch->k_0,
                  system->t, system->x,
                  x_idx, system->params);
    }

    // x(t + h) = x(t) + h k[0]
    for (x_idx = 0; x_idx < system->dims; x_idx++) {
        prec_x = arpra_get_precision(&(system->x[x_idx]));
        arpra_set_precision(&(scratch->temp_x), prec_x);
        arpra_mul(&(scratch->temp_x), h, &(scratch->k_0[x_idx]));
        arpra_add(&(scratch->x_new[x_idx]), &(system->x[x_idx]), &(scratch->temp_x));
    }

    // Advance system.
    arpra_add(system->t, system->t, h);
    for (x_idx = 0; x_idx < system->dims; x_idx++) {
        arpra_set(&(system->x[x_idx]), &(scratch->x_new[x_idx]));
    }
}

static const arpra_ode_method euler =
{
    .init = &euler_init,
    .clear = &euler_clear,
    .step = &euler_step,
    .stages = euler_stages,
};

const arpra_ode_method *arpra_ode_euler = &euler;
