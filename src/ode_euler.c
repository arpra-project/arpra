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
    arpra_range k1;
} euler_scratch;

static void euler_init (arpra_ode_stepper *stepper, arpra_ode_system *system)
{
    arpra_uint i;
    arpra_precision prec;
    euler_scratch *scratch;

    stepper->method = arpra_ode_euler;
    stepper->system = system;
    stepper->scratch = malloc(system->dims * sizeof(euler_scratch));
    scratch = (euler_scratch *) stepper->scratch;
    for (i = 0; i < system->dims; i++) {
        prec = arpra_get_precision(&(system->x[i]));
        arpra_init2(&(scratch[i].k1), prec);
    }
}

static void euler_clear (arpra_ode_stepper *stepper)
{
    arpra_uint i;
    euler_scratch *scratch;

    scratch = (euler_scratch *) stepper->scratch;
    for (i = 0; i < stepper->system->dims; i++) {
        arpra_clear(&(scratch[i].k1));
    }
    free(scratch);
}

static void euler_reset (arpra_ode_stepper *stepper)
{
    arpra_uint i;
    euler_scratch *scratch;

    scratch = (euler_scratch *) stepper->scratch;
    for (i = 0; i < stepper->system->dims; i++) {
        arpra_set_zero(&(scratch[i].k1));
    }
}

static void euler_step (arpra_ode_stepper *stepper, const arpra_range *h)
{
    arpra_uint i;
    euler_scratch *scratch;

}

static const arpra_ode_method euler =
{
    .init = &euler_init,
    .clear = &euler_clear,
    .reset = &euler_reset,
    .step = &euler_step
};

const arpra_ode_method *arpra_ode_euler = &euler;
