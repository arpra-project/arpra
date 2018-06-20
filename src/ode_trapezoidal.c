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

typedef struct trapezoidal_scratch_struct
{
    arpra_range *k1;
    arpra_range *k2;
    arpra_range *temp_x;
    arpra_range temp;
} trapezoidal_scratch;

static void trapezoidal_init (arpra_ode_stepper *stepper, arpra_ode_system *system)
{
    arpra_uint i;
    arpra_precision prec;
    trapezoidal_scratch *scratch;

    scratch = malloc(sizeof(trapezoidal_scratch));
    scratch->k1 = malloc(system->dims * sizeof(arpra_range));
    scratch->k2 = malloc(system->dims * sizeof(arpra_range));
    scratch->temp_x = malloc(system->dims * sizeof(arpra_range));
    for (i = 0; i < system->dims; i++) {
        prec = arpra_get_precision(&(system->x[i]));
        arpra_init2(&(scratch->k1[i]), prec);
        arpra_init2(&(scratch->k2[i]), prec);
        arpra_init2(&(scratch->temp_x[i]), prec);
    }
    prec = arpra_get_default_precision();
    arpra_init2(&scratch->temp, prec);
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
    for (i = 0; i < system->dims; i++) {
        arpra_clear(&(scratch->k1[i]));
        arpra_clear(&(scratch->k2[i]));
        arpra_clear(&(scratch->temp_x[i]));
    }
    arpra_clear(&scratch->temp);
    free(scratch);
}

static void trapezoidal_reset (arpra_ode_stepper *stepper)
{
    arpra_uint i;
    arpra_ode_system *system;
    trapezoidal_scratch *scratch;

    system = stepper->system;
    scratch = (trapezoidal_scratch *) stepper->scratch;
    for (i = 0; i < system->dims; i++) {
        arpra_set_zero(&(scratch->k1[i]));
        arpra_set_zero(&(scratch->k2[i]));
        arpra_set_zero(&(scratch->temp_x[i]));
    }
    arpra_set_zero(&scratch->temp);
}

static void trapezoidal_step (arpra_ode_stepper *stepper, const arpra_range *h)
{
    arpra_uint i;
    arpra_precision prec;
    arpra_ode_system *system;
    trapezoidal_scratch *scratch;

    system = stepper->system;
    scratch = (trapezoidal_scratch *) stepper->scratch;

    // Compute k1.
    system->f(scratch->k1,
              system->x, system->t,
              system->dims, system->params);

    // Compute k2.
    arpra_add(system->t, system->t, h);
    arpra_mul(&(scratch->temp_x[i]), h, &(scratch->k1[i]));
    arpra_add(&(scratch->temp_x[i]), &(system->x[i]), &(scratch->temp_x[i]));
    system->f(scratch->k2,
              scratch->temp_x, system->t,
              system->dims, system->params);

    // Step x by h.
    for (i = 0; i < system->dims; i++) {
        prec = arpra_get_precision(&(system->x[i]));
        arpra_set_precision(&scratch->temp, prec);

        // 1/2 k1 + 1/2 k2
        arpra_set_d(&scratch->temp, 0.5);
        arpra_mul(&(scratch->temp_x[i]), &(scratch->k1[i]), &scratch->temp);
        arpra_mul(&scratch->temp, &(scratch->k2[i]), &scratch->temp);
        arpra_add(&(scratch->temp_x[i]), &(scratch->temp_x[i]), &scratch->temp);

        // Scale by h and step.
        arpra_mul(&(scratch->temp_x[i]), h, &(scratch->temp_x[i]));
        arpra_add(&(system->x[i]), &(system->x[i]), &(scratch->temp_x[i]));
    }
}

static const arpra_ode_method trapezoidal =
{
    .init = &trapezoidal_init,
    .clear = &trapezoidal_clear,
    .reset = &trapezoidal_reset,
    .step = &trapezoidal_step
};

const arpra_ode_method *arpra_ode_trapezoidal = &trapezoidal;
