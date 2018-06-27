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
    arpra_range *k_1;
    arpra_range *k_2;
    arpra_range *k_weights;
    arpra_range *next_t;
    arpra_range *next_x;
    arpra_range *temp;
} trapezoidal_scratch;

static void trapezoidal_init (arpra_ode_stepper *stepper, arpra_ode_system *system)
{
    arpra_uint i;
    arpra_precision prec;
    trapezoidal_scratch *scratch;

    scratch = malloc(sizeof(trapezoidal_scratch));
    scratch->k_1 = malloc(system->dims * sizeof(arpra_range));
    scratch->k_2 = malloc(system->dims * sizeof(arpra_range));
    scratch->k_weights = malloc(2 * sizeof(arpra_range));
    scratch->next_t = malloc(sizeof(arpra_range));
    scratch->next_x = malloc(system->dims * sizeof(arpra_range));
    scratch->temp = malloc(sizeof(arpra_range));
    for (i = 0; i < system->dims; i++) {
        prec = arpra_get_precision(&(system->x[i]));
        arpra_init2(&(scratch->k_1[i]), prec);
        arpra_init2(&(scratch->k_2[i]), prec);
        arpra_init2(&(scratch->next_x[i]), prec);
    }
    prec = arpra_get_default_precision();
    for (i = 0; i < 2; i++) {
        arpra_init2(&(scratch->k_weights[i]), prec);
    }
    arpra_init2(scratch->next_t, prec);
    arpra_init2(scratch->temp, prec);
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
        arpra_clear(&(scratch->k_1[i]));
        arpra_clear(&(scratch->k_2[i]));
        arpra_clear(&(scratch->next_x[i]));
    }
    for (i = 0; i < 2; i++) {
        arpra_clear(&(scratch->k_weights[i]));
    }
    arpra_clear(scratch->next_t);
    arpra_clear(scratch->temp);
    free(scratch->k_1);
    free(scratch->k_2);
    free(scratch->k_weights);
    free(scratch->next_t);
    free(scratch->next_x);
    free(scratch->temp);
    free(scratch);
}

static void trapezoidal_step (arpra_ode_stepper *stepper, const arpra_range *h)
{
    arpra_uint i;
    arpra_precision prec;
    arpra_ode_system *system;
    trapezoidal_scratch *scratch;

    system = stepper->system;
    scratch = (trapezoidal_scratch *) stepper->scratch;

    // Synchronise scratch memory precision.
    for (i = 0; i < system->dims; i++) {
        prec = arpra_get_precision(&(system->x[i]));
        arpra_set_precision(&(scratch->k_1[i]), prec);
        arpra_set_precision(&(scratch->k_2[i]), prec);
        arpra_set_precision(&(scratch->next_x[i]), prec);
    }
    prec = arpra_get_precision(system->t);
    for (i = 0; i < 2; i++) {
        arpra_set_precision(&(scratch->k_weights[i]), prec);
    }
    arpra_set_precision(scratch->next_t, prec);
    arpra_set_precision(scratch->temp, prec);

    // k_1 = f([t], [x(t)])
    system->f(scratch->k_1,
              system->t, system->x,
              system->dims, system->params);

    // k_2 = f([t + h], [x(t) + h k_1])
    for (i = 0; i < system->dims; i++) {
        arpra_mul(&(scratch->next_x[i]), h, &(scratch->k_1[i]));
        arpra_add(&(scratch->next_x[i]), &(system->x[i]), &(scratch->next_x[i]));
    }
    arpra_add(scratch->next_t, system->t, h);
    system->f(scratch->k_2,
              scratch->next_t, scratch->next_x,
              system->dims, system->params);

    // x(t + h) = x(t) + h/2 k_1 + h/2 k_2
    arpra_set_d(scratch->temp, 2.0);
    arpra_div(&(scratch->k_weights[0]), h, scratch->temp);
    arpra_set(&(scratch->k_weights[1]), &(scratch->k_weights[0]));
    for (i = 0; i < system->dims; i++) {
        prec = arpra_get_precision(&(system->x[i]));
        arpra_set_precision(scratch->temp, prec);
        arpra_mul(&(scratch->next_x[i]), &(scratch->k_weights[0]), &(scratch->k_1[i]));
        arpra_mul(scratch->temp, &(scratch->k_weights[0]), &(scratch->k_2[i]));
        arpra_add(&(scratch->next_x[i]), &(scratch->next_x[i]), scratch->temp);
        arpra_add(&(system->x[i]), &(system->x[i]), &(scratch->next_x[i]));
    }

    // Advance t.
    arpra_add(system->t, system->t, h);
}

static const arpra_ode_method trapezoidal =
{
    .init = &trapezoidal_init,
    .clear = &trapezoidal_clear,
    .step = &trapezoidal_step,
    .stages = 2,
};

const arpra_ode_method *arpra_ode_trapezoidal = &trapezoidal;
