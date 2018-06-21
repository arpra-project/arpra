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

typedef struct euler_scratch_struct
{
    arpra_range *k_1;
    arpra_range *temp;
} euler_scratch;

static void euler_init (arpra_ode_stepper *stepper, arpra_ode_system *system)
{
    arpra_uint i;
    arpra_precision prec;
    euler_scratch *scratch;

    scratch = malloc(sizeof(euler_scratch));
    scratch->k_1 = malloc(system->dims * sizeof(arpra_range));
    scratch->temp = malloc(sizeof(arpra_range));
    for (i = 0; i < system->dims; i++) {
        prec = arpra_get_precision(&(system->x[i]));
        arpra_init2(&(scratch->k_1[i]), prec);
    }
    prec = arpra_get_default_precision();
    arpra_init2(scratch->temp, prec);
    stepper->method = arpra_ode_euler;
    stepper->system = system;
    stepper->error = NULL;
    stepper->scratch = scratch;
}

static void euler_clear (arpra_ode_stepper *stepper)
{
    arpra_uint i;
    arpra_ode_system *system;
    euler_scratch *scratch;

    system = stepper->system;
    scratch = (euler_scratch *) stepper->scratch;
    for (i = 0; i < system->dims; i++) {
        arpra_clear(&(scratch->k_1[i]));
    }
    arpra_clear(scratch->temp);
    free(scratch->k_1);
    free(scratch->temp);
    free(scratch);
}

static void euler_step (arpra_ode_stepper *stepper, const arpra_range *h)
{
    arpra_uint i;
    arpra_precision prec;
    arpra_ode_system *system;
    euler_scratch *scratch;

    system = stepper->system;
    scratch = (euler_scratch *) stepper->scratch;

    // Synchronise scratch memory precision.
    for (i = 0; i < system->dims; i++) {
        prec = arpra_get_precision(&(system->x[i]));
        arpra_set_precision(&(scratch->k_1[i]), prec);
    }

    // k_1 = f([t], [x(t)])
    system->f(scratch->k_1,
              system->t, system->x,
              system->dims, system->params);

    // x(t + h) = x(t) + h k_1
    arpra_add(system->t, system->t, h);
    for (i = 0; i < system->dims; i++) {
        prec = arpra_get_precision(&(system->x[i]));
        arpra_set_precision(scratch->temp, prec);
        arpra_mul(scratch->temp, h, &(scratch->k_1[i]));
        arpra_add(&(system->x[i]), &(system->x[i]), scratch->temp);
    }
}

static const arpra_ode_method euler =
{
    .init = &euler_init,
    .clear = &euler_clear,
    .step = &euler_step
};

const arpra_ode_method *arpra_ode_euler = &euler;
