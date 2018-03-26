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

struct rk2_state
{
    struct arpra_range k1;
    struct arpra_range k2;
    struct arpra_range k3;
    struct arpra_range temp;
};


static void rk2_init (struct arpra_ode_stepper *stepper, const struct arpra_ode_system *system)
{
    struct rk2_state *state;

    stepper->system = system;
    state = (struct rk2_state *) stepper->state;
    state = malloc(sizeof(struct rk2_state));
    arpra_init(&(state->k1));
    arpra_init(&(state->k2));
    arpra_init(&(state->k3));
    arpra_init(&(state->temp));
}


static void rk2_init2 (struct arpra_ode_stepper *stepper, const struct arpra_ode_system *system)
{
    struct rk2_state state;

}


static void rk2_clear (struct arpra_ode_stepper *stepper)
{
    struct rk2_state *state;

    state = (struct rk2_state *) stepper->state;
    arpra_clear(&(state->k1));
    arpra_clear(&(state->k2));
    arpra_clear(&(state->k3));
    arpra_clear(&(state->temp));
    free(state);
}


static void rk2_reset (struct arpra_ode_stepper *stepper)
{
    struct rk2_state *state;

    state = (struct rk2_state *) stepper->state;
    arpra_set_zero(&(state->k1));
    arpra_set_zero(&(state->k2));
    arpra_set_zero(&(state->k3));
    arpra_set_zero(&(state->temp));
}


static void rk2_step (const struct arpra_ode_stepper *stepper,
                      struct arpra_range *dxdt, struct arpra_range *t, struct arpra_range *error,
                      struct arpra_range *x, const struct arpra_range *h)
{
    struct rk2_state *state;

    state = (struct rk2_state *) stepper->state;

    // EVALUATE K1 = DYDT|T
    stepper->system->function (t, x, &(state->k1), stepper->system->parameters);

}


static const struct arpra_ode_stepper rk2 = {};

const struct arpra_ode_stepper *arpra_ode_rk2 = &rk2;
