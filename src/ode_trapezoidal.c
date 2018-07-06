/*
 * ode_trapezoidal.c -- Explicit Trapezoidal Rule ODE stepper.
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

#define trapezoidal_stages 2

typedef struct trapezoidal_scratch_struct
{ 
    arpra_range *k_0;
    arpra_range *k_1;
    arpra_range *x_new;
    arpra_range temp;
    arpra_range half_h;
    arpra_range half;
} trapezoidal_scratch;

static void trapezoidal_init (arpra_ode_stepper *stepper, arpra_ode_system *system)
{
    arpra_uint i;
    arpra_precision prec_x, prec_internal;
    trapezoidal_scratch *scratch;

    // Allocate scratch memory.
    scratch = malloc(sizeof(trapezoidal_scratch));
    scratch->k_0 = malloc(system->dims * sizeof(arpra_range));
    scratch->k_1 = malloc(system->dims * sizeof(arpra_range));
    scratch->x_new = malloc(system->dims * sizeof(arpra_range));

    // Initialise scratch memory.
    for (i = 0; i < system->dims; i++) {
        prec_x = arpra_get_precision(&(system->x[i]));
        arpra_init2(&(scratch->k_0[i]), prec_x);
        arpra_init2(&(scratch->k_1[i]), prec_x);
        arpra_init2(&(scratch->x_new[i]), prec_x);
    }
    prec_internal = arpra_get_internal_precision();
    arpra_init2(&(scratch->temp), prec_internal);
    arpra_init2(&(scratch->half_h), prec_internal);
    arpra_init2(&(scratch->half), 2);

    // Precompute constants.
    arpra_set_d(&(scratch->half), 0.5);

    // Set stepper parameters.
    stepper->method = arpra_ode_trapezoidal;
    stepper->system = system;
    stepper->error = NULL;
    stepper->scratch = scratch;
}

static void trapezoidal_clear (arpra_ode_stepper *stepper)
{
    arpra_uint i;
    arpra_ode_system *system;
    trapezoidal_scratch *scratch;

    system = stepper->system;
    scratch = (trapezoidal_scratch *) stepper->scratch;

    // Clear scratch memory.
    for (i = 0; i < system->dims; i++) {
        arpra_clear(&(scratch->k_0[i]));
        arpra_clear(&(scratch->k_1[i]));
        arpra_clear(&(scratch->x_new[i]));
    }
    arpra_clear(&(scratch->temp));
    arpra_clear(&(scratch->half_h));
    arpra_clear(&(scratch->half));

    // Free scratch memory.
    free(scratch->k_0);
    free(scratch->k_1);
    free(scratch->x_new);
    free(scratch);
}

static void trapezoidal_step (arpra_ode_stepper *stepper, const arpra_range *h)
{
    arpra_uint i;
    arpra_precision prec_t, prec_x;
    arpra_ode_system *system;
    trapezoidal_scratch *scratch;

    system = stepper->system;
    scratch = (trapezoidal_scratch *) stepper->scratch;

    // Synchronise scratch memory precision.
    for (i = 0; i < system->dims; i++) {
        prec_x = arpra_get_precision(&(system->x[i]));
        arpra_set_precision(&(scratch->k_0[i]), prec_x);
        arpra_set_precision(&(scratch->k_1[i]), prec_x);
        arpra_set_precision(&(scratch->x_new[i]), prec_x);
    }
    prec_t = arpra_get_precision(system->t);
    arpra_set_precision(&(scratch->half_h), prec_t);

    // k_0 = f(t, x(t))
    system->f(scratch->k_0,
              system->t, system->x,
              system->dims, system->params);

    // k_1 = f(t + h, x(t) + h k_0)
    for (i = 0; i < system->dims; i++) {
        prec_x = arpra_get_precision(&(system->x[i]));
        arpra_set_precision(&(scratch->temp), prec_x);
        arpra_mul(&(scratch->temp), h, &(scratch->k_0[i]));
        arpra_add(&(scratch->x_new[i]), &(system->x[i]), &(scratch->temp));
    }
    arpra_set_precision(&(scratch->temp), prec_t);
    arpra_add(&(scratch->temp), system->t, h);
    system->f(scratch->k_1,
              &(scratch->temp), scratch->x_new,
              system->dims, system->params);

    // x(t + h) = x(t) + h/2 k_0 + h/2 k_1
    arpra_mul(&(scratch->half_h), &(scratch->half), h);
    for (i = 0; i < system->dims; i++) {
        prec_x = arpra_get_precision(&(system->x[i]));
        arpra_set_precision(&(scratch->temp), prec_x);
        arpra_mul(&(scratch->temp), &(scratch->half_h), &(scratch->k_0[i]));
        arpra_add(&(scratch->x_new[i]), &(system->x[i]), &(scratch->temp));
        arpra_mul(&(scratch->temp), &(scratch->half_h), &(scratch->k_1[i]));
        arpra_add(&(scratch->x_new[i]), &(scratch->x_new[i]), &(scratch->temp));
    }

    // Advance system.
    arpra_add(system->t, system->t, h);
    arpra_range *x_temp = system->x;
    system->x = scratch->x_new;
    scratch->x_new = x_temp;
}

static const arpra_ode_method trapezoidal =
{
    .init = &trapezoidal_init,
    .clear = &trapezoidal_clear,
    .step = &trapezoidal_step,
    .stages = trapezoidal_stages,
};

const arpra_ode_method *arpra_ode_trapezoidal = &trapezoidal;
