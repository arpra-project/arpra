/*
 * ode_rk2.c -- Runge-Kutta 2nd order ODE stepper.
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

typedef struct rk2_scratch_struct
{
    arpra_range *k1;
    arpra_range *k2;
    arpra_range *k3;
} rk2_scratch;

static void rk2_init (arpra_ode_stepper *stepper, arpra_ode_system *system)
{
    arpra_uint i;
    arpra_precision prec;
    rk2_scratch *scratch;

    scratch = malloc(sizeof(rk2_scratch));
    scratch->k1 = malloc(system->dims * sizeof(arpra_range));
    scratch->k2 = malloc(system->dims * sizeof(arpra_range));
    scratch->k3 = malloc(system->dims * sizeof(arpra_range));
    for (i = 0; i < system->dims; i++) {
        prec = arpra_get_precision(&(system->x[i]));
        arpra_init2(&(scratch->k1[i]), prec);
        arpra_init2(&(scratch->k2[i]), prec);
        arpra_init2(&(scratch->k3[i]), prec);
    }
    stepper->method = arpra_ode_rk2;
    stepper->system = system;
    stepper->scratch = scratch;
}

static void rk2_clear (arpra_ode_stepper *stepper)
{
    arpra_uint i;
    rk2_scratch *scratch;

    scratch = (rk2_scratch *) stepper->scratch;
    for (i = 0; i < stepper->system->dims; i++) {
        arpra_clear(&(scratch->k1[i]));
        arpra_clear(&(scratch->k2[i]));
        arpra_clear(&(scratch->k3[i]));
    }
    free(scratch->k1);
    free(scratch->k2);
    free(scratch->k3);
    free(scratch);
}

static void rk2_reset (arpra_ode_stepper *stepper)
{
    arpra_uint i;
    rk2_scratch *scratch;

    scratch = (rk2_scratch *) stepper->scratch;
    for (i = 0; i < stepper->system->dims; i++) {
        arpra_set_zero(&(scratch->k1[i]));
        arpra_set_zero(&(scratch->k2[i]));
        arpra_set_zero(&(scratch->k3[i]));
    }
}

static void rk2_step (arpra_ode_stepper *stepper, const arpra_range *h)
{
    arpra_uint i;
    rk2_scratch *scratch;

    scratch = (rk2_scratch *) stepper->scratch;

    // ======================== SCRATCH SPACE SHOULD HOLD ERROR, IF NEEDED

    // Compute k1.
    stepper->system->f(scratch->k1,
                       stepper->system->x, stepper->system->t,
                       stepper->system->dims, stepper->system->params);

    // Compute k2.


    // Compute k3.

    // Step x by h.
}

static const arpra_ode_method rk2 =
{
    .init = &rk2_init,
    .clear = &rk2_clear,
    .reset = &rk2_reset,
    .step = &rk2_step
};

const arpra_ode_method *arpra_ode_rk2 = &rk2;
