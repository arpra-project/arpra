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

typedef struct
{
    arpra_range k1;
    arpra_range k2;
    arpra_range k3;
    arpra_range temp;
} rk2_workspace;


static void rk2_init (arpra_ode_stepper *stepper)
{
    rk2_workspace *workspace;

    workspace = (rk2_workspace *) stepper->workspace;
    workspace = malloc(sizeof(rk2_workspace));
    arpra_init(&(workspace->k1));
    arpra_init(&(workspace->k2));
    arpra_init(&(workspace->k3));
    arpra_init(&(workspace->temp));
}


static void rk2_init2 (arpra_ode_stepper *stepper, const arpra_precision prec)
{
    rk2_workspace *workspace;

}


static void rk2_clear (arpra_ode_stepper *stepper)
{
    rk2_workspace *workspace;

    workspace = (rk2_workspace *) stepper->workspace;
    arpra_clear(&(workspace->k1));
    arpra_clear(&(workspace->k2));
    arpra_clear(&(workspace->k3));
    arpra_clear(&(workspace->temp));
    free(workspace);
}


static void rk2_reset (arpra_ode_stepper *stepper)
{
    rk2_workspace *workspace;

    workspace = (rk2_workspace *) stepper->workspace;
    arpra_set_zero(&(workspace->k1));
    arpra_set_zero(&(workspace->k2));
    arpra_set_zero(&(workspace->k3));
    arpra_set_zero(&(workspace->temp));
}


static void rk2_step (arpra_ode_stepper *stepper, const arpra_range *h)
{
    rk2_workspace *workspace;

    workspace = (rk2_workspace *) stepper->workspace;

    // EVALUATE K1 = DYDT|T
    //stepper->function (t, x, &(workspace->k1), stepper->parameters);
}

// stepper doesnt need to save old state - can be done by the step control object

static const arpra_ode_method rk2 = {};

const arpra_ode_method *arpra_ode_rk2 = &rk2;
