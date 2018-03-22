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

struct __rk2_state_struct
{
    struct __arpra_struct k1;
    struct __arpra_struct k2;
    struct __arpra_struct k3;
    struct __arpra_struct temp;
};

typedef struct __rk2_state_struct rk2_state_t;
typedef struct __rk2_state_struct *rk2_state_ptr;
typedef const struct __rk2_state_struct *rk2_state_srcptr;


static void rk2_init (arpra_ode_stepper_ptr stepper, arpra_ode_system_srcptr system)
{
    rk2_state_ptr state_rk2;

    stepper->system = system;
    state_rk2 = (rk2_state_ptr) stepper->state;
    state_rk2 = malloc(sizeof(rk2_state_t));
    arpra_init(&(state_rk2->k1));
    arpra_init(&(state_rk2->k2));
    arpra_init(&(state_rk2->k3));
    arpra_init(&(state_rk2->temp));
}


static void rk2_init2 (arpra_ode_stepper_ptr stepper, arpra_ode_system_srcptr system)
{
    rk2_state_ptr state_rk2;

}


static void rk2_clear (arpra_ode_stepper_ptr stepper)
{
    rk2_state_ptr state_rk2;

    state_rk2 = (rk2_state_ptr) stepper->state;
    arpra_clear(&(state_rk2->k1));
    arpra_clear(&(state_rk2->k2));
    arpra_clear(&(state_rk2->k3));
    arpra_clear(&(state_rk2->temp));
    free(state_rk2);
}


static void rk2_reset (arpra_ode_stepper_ptr stepper)
{
    rk2_state_ptr state_rk2;

    state_rk2 = (rk2_state_ptr) stepper->state;
    arpra_set_zero(&(state_rk2->k1));
    arpra_set_zero(&(state_rk2->k2));
    arpra_set_zero(&(state_rk2->k3));
    arpra_set_zero(&(state_rk2->temp));
}


static void rk2_step (arpra_ode_stepper_srcptr stepper,
                      arpra_ptr dxdt, arpra_ptr t, arpra_ptr error,
                      arpra_srcptr x, arpra_srcptr h)
{
    rk2_state_ptr state_rk2;

    state_rk2 = (rk2_state_ptr) stepper->state;

    // EVALUATE K1 = DYDT|T
    stepper->system->function (t, x, &(state_rk2->k1), stepper->system->parameters);

}


static const arpra_ode_stepper_t rk2 = {};

arpra_ode_stepper_srcptr arpra_ode_rk2 = &rk2;
