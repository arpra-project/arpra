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
    arpra_range half;
    arpra_range half_h;
    arpra_range temp_t;
    arpra_range temp_x;
} trapezoidal_scratch;

static void trapezoidal_init (arpra_ode_stepper *stepper, arpra_ode_system *system)
{
    arpra_uint x_i;
    arpra_precision prec_x, prec_internal;
    trapezoidal_scratch *scratch;

    // Allocate scratch memory.
    scratch = malloc(sizeof(trapezoidal_scratch));
    scratch->k_0 = malloc(system->dims * sizeof(arpra_range));
    scratch->k_1 = malloc(system->dims * sizeof(arpra_range));
    scratch->x_new = malloc(system->dims * sizeof(arpra_range));

    // Initialise scratch memory.
    prec_internal = arpra_get_internal_precision();
    for (x_i = 0; x_i < system->dims; x_i++) {
        prec_x = arpra_get_precision(&(system->x[x_i]));
        arpra_init2(&(scratch->k_0[x_i]), prec_x);
        arpra_init2(&(scratch->k_1[x_i]), prec_x);
        arpra_init2(&(scratch->x_new[x_i]), prec_x);
    }
    arpra_init2(&(scratch->half), 2);
    arpra_init2(&(scratch->half_h), prec_internal);
    arpra_init2(&(scratch->temp_t), prec_internal);
    arpra_init2(&(scratch->temp_x), prec_internal);

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
    arpra_uint x_i;
    arpra_ode_system *system;
    trapezoidal_scratch *scratch;

    system = stepper->system;
    scratch = (trapezoidal_scratch *) stepper->scratch;

    // Clear scratch memory.
    for (x_i = 0; x_i < system->dims; x_i++) {
        arpra_clear(&(scratch->k_0[x_i]));
        arpra_clear(&(scratch->k_1[x_i]));
        arpra_clear(&(scratch->x_new[x_i]));
    }
    arpra_clear(&(scratch->half));
    arpra_clear(&(scratch->half_h));
    arpra_clear(&(scratch->temp_t));
    arpra_clear(&(scratch->temp_x));

    // Free scratch memory.
    free(scratch->k_0);
    free(scratch->k_1);
    free(scratch->x_new);
    free(scratch);
}

static void trapezoidal_step (arpra_ode_stepper *stepper, const arpra_range *h)
{
    arpra_uint x_i;
    arpra_precision prec_t, prec_x;
    arpra_range *x_new;
    arpra_ode_system *system;
    trapezoidal_scratch *scratch;

    system = stepper->system;
    scratch = (trapezoidal_scratch *) stepper->scratch;

    // Synchronise scratch precision and prepare step parameters.
    prec_t = arpra_get_precision(system->t);
    for (x_i = 0; x_i < system->dims; x_i++) {
        prec_x = arpra_get_precision(&(system->x[x_i]));
        arpra_set_precision(&(scratch->k_0[x_i]), prec_x);
        arpra_set_precision(&(scratch->k_1[x_i]), prec_x);
        arpra_set_precision(&(scratch->x_new[x_i]), prec_x);
    }
    arpra_set_precision(&(scratch->half_h), prec_t);
    arpra_mul(&(scratch->half_h), &(scratch->half), h);
    arpra_set_precision(&(scratch->temp_t), prec_t);
    arpra_add(&(scratch->temp_t), system->t, h);

    // k[0] = f(t, x(t))
    for (x_i = 0; x_i < system->dims; x_i++) {
        system->f(scratch->k_0,
                  system->t, system->x,
                  x_i, system->params);
    }

    // x(t + h) = x(t) + h k[0]
    for (x_i = 0; x_i < system->dims; x_i++) {
        prec_x = arpra_get_precision(&(system->x[x_i]));
        arpra_set_precision(&(scratch->temp_x), prec_x);
        arpra_mul(&(scratch->temp_x), h, &(scratch->k_0[x_i]));
        arpra_add(&(scratch->x_new[x_i]), &(system->x[x_i]), &(scratch->temp_x));
    }

    // k[1] = f(t + h, x(t) + h k[0])
    for (x_i = 0; x_i < system->dims; x_i++) {
        system->f(scratch->k_1,
                  &(scratch->temp_t), scratch->x_new,
                  x_i, system->params);
    }

    // x(t + h) = x(t) + 1/2 h k[0] + 1/2 h k[1]
    for (x_i = 0; x_i < system->dims; x_i++) {
        prec_x = arpra_get_precision(&(system->x[x_i]));
        arpra_set_precision(&(scratch->temp_x), prec_x);
        arpra_mul(&(scratch->temp_x), &(scratch->half_h), &(scratch->k_0[x_i]));
        arpra_add(&(scratch->x_new[x_i]), &(system->x[x_i]), &(scratch->temp_x));
        arpra_mul(&(scratch->temp_x), &(scratch->half_h), &(scratch->k_1[x_i]));
        arpra_add(&(scratch->x_new[x_i]), &(scratch->x_new[x_i]), &(scratch->temp_x));
    }

    // Advance system.
    arpra_add(system->t, system->t, h);
    x_new = scratch->x_new;
    scratch->x_new = system->x;
    system->x = x_new;
}

static const arpra_ode_method trapezoidal =
{
    .init = &trapezoidal_init,
    .clear = &trapezoidal_clear,
    .step = &trapezoidal_step,
    .stages = trapezoidal_stages,
};

const arpra_ode_method *arpra_ode_trapezoidal = &trapezoidal;
