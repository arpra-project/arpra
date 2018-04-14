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
    arpra_range k1;
    arpra_range k2;
    arpra_range k3;
    arpra_range temp;
} rk2_scratch;

static void rk2_init (arpra_ode_stepper *stepper, arpra_ode_system *system)
{
    arpra_uint i;
    arpra_precision prec;
    rk2_scratch *scratch;

    stepper->method = arpra_ode_rk2;
    stepper->system = system;
    stepper->scratch = malloc(system->dims * sizeof(rk2_scratch));
    scratch = (rk2_scratch *) stepper->scratch;
    for (i = 0; i < system->dims; i++) {
        prec = arpra_get_precision(&(system->x[i]));
        arpra_init2(&(scratch[i].k1), prec);
        arpra_init2(&(scratch[i].k2), prec);
        arpra_init2(&(scratch[i].k3), prec);
        arpra_init2(&(scratch[i].temp), prec);
    }
}

static void rk2_clear (arpra_ode_stepper *stepper)
{
    arpra_uint i;
    rk2_scratch *scratch;

    scratch = (rk2_scratch *) stepper->scratch;
    for (i = 0; i < stepper->system->dims; i++) {
        arpra_clear(&(scratch[i].k1));
        arpra_clear(&(scratch[i].k2));
        arpra_clear(&(scratch[i].k3));
        arpra_clear(&(scratch[i].temp));
    }
    free(scratch);
}

static void rk2_reset (arpra_ode_stepper *stepper)
{
    arpra_uint i;
    rk2_scratch *scratch;

    scratch = (rk2_scratch *) stepper->scratch;
    for (i = 0; i < stepper->system->dims; i++) {
        arpra_set_zero(&(scratch[i].k1));
        arpra_set_zero(&(scratch[i].k2));
        arpra_set_zero(&(scratch[i].k3));
        arpra_set_zero(&(scratch[i].temp));
    }
}

static void rk2_step (arpra_ode_stepper *stepper, const arpra_range *h)
{
    arpra_uint i;
    rk2_scratch *scratch;

    // EVALUATE K1 = DYDT|T
    //stepper->system->f(t, x, &(scratch->k1), stepper->params);
}

static const arpra_ode_method rk2 =
{
    .init = &rk2_init,
    .clear = &rk2_clear,
    .reset = &rk2_reset,
    .step = &rk2_step
};

const arpra_ode_method *arpra_ode_rk2 = &rk2;
