/*
 * ode_euler.c -- Explicit Euler ODE stepper.
 *
 * Copyright 2018-2020 James Paul Turner.
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

#define euler_stages 1

typedef struct euler_scratch_struct
{
    arpra_range *_k_0;
    arpra_range **k_0;
    arpra_range *_x_new;
    arpra_range **x_new;
    arpra_range temp_x;
} euler_scratch;

static void euler_init (arpra_ode_stepper *stepper, arpra_ode_system *system)
{
    arpra_uint x_grp, x_dim, state_size;
    arpra_prec prec_x, prec_internal;
    euler_scratch *scratch;

    // Allocate scratch memory.
    scratch = malloc(sizeof(euler_scratch));
    for (x_grp = 0, state_size = 0; x_grp < system->grps; x_grp++) {
        state_size += system->dims[x_grp];
    }
    scratch->_k_0 = malloc(state_size * sizeof(arpra_range));
    scratch->k_0 = malloc(system->grps * sizeof(arpra_range *));
    scratch->_x_new = malloc(state_size * sizeof(arpra_range));
    scratch->x_new = malloc(system->grps * sizeof(arpra_range *));

    // Initialise scratch memory.
    prec_internal = arpra_get_internal_precision();
    scratch->k_0[0] = scratch->_k_0;
    scratch->x_new[0] = scratch->_x_new;
    for (x_grp = 1; x_grp < system->grps; x_grp++) {
        scratch->k_0[x_grp] = scratch->k_0[x_grp - 1] + system->dims[x_grp - 1];
        scratch->x_new[x_grp] = scratch->x_new[x_grp - 1] + system->dims[x_grp - 1];
    }
    for (x_grp = 0; x_grp < system->grps; x_grp++) {
        for (x_dim = 0; x_dim < system->dims[x_grp]; x_dim++) {
            prec_x = arpra_get_precision(&(system->x[x_grp][x_dim]));
            arpra_init2(&(scratch->k_0[x_grp][x_dim]), prec_x);
            arpra_init2(&(scratch->x_new[x_grp][x_dim]), prec_x);
        }
    }
    arpra_init2(&(scratch->temp_x), prec_internal);

    // Set stepper parameters.
    stepper->method = arpra_ode_euler;
    stepper->system = system;
    stepper->error = NULL;
    stepper->scratch = scratch;
}

static void euler_clear (arpra_ode_stepper *stepper)
{
    arpra_uint x_grp, x_dim;
    arpra_ode_system *system;
    euler_scratch *scratch;

    system = stepper->system;
    scratch = (euler_scratch *) stepper->scratch;

    // Clear scratch memory.
    for (x_grp = 0; x_grp < system->grps; x_grp++) {
        for (x_dim = 0; x_dim < system->dims[x_grp]; x_dim++) {
            arpra_clear(&(scratch->k_0[x_grp][x_dim]));
            arpra_clear(&(scratch->x_new[x_grp][x_dim]));
        }
    }
    arpra_clear(&(scratch->temp_x));

    // Free scratch memory.
    free(scratch->_k_0);
    free(scratch->k_0);
    free(scratch->_x_new);
    free(scratch->x_new);
    free(scratch);
}

static void euler_step (arpra_ode_stepper *stepper, const arpra_range *h)
{
    arpra_uint x_grp, x_dim;
    arpra_prec prec_x;
    arpra_ode_system *system;
    euler_scratch *scratch;

    system = stepper->system;
    scratch = (euler_scratch *) stepper->scratch;

    // Synchronise scratch precision and prepare step parameters.
    for (x_grp = 0; x_grp < system->grps; x_grp++) {
        for (x_dim = 0; x_dim < system->dims[x_grp]; x_dim++) {
            prec_x = arpra_get_precision(&(system->x[x_grp][x_dim]));
            arpra_set_precision(&(scratch->k_0[x_grp][x_dim]), prec_x);
            arpra_set_precision(&(scratch->x_new[x_grp][x_dim]), prec_x);
        }
    }

    // k[0] = f(t, x(t))
    for (x_grp = 0; x_grp < system->grps; x_grp++) {
        for (x_dim = 0; x_dim < system->dims[x_grp]; x_dim++) {
            system->f[x_grp](&(scratch->k_0[x_grp][x_dim]), system->params[x_grp],
                             system->t, (const arpra_range **) system->x, x_grp, x_dim);
        }
    }

    // x(t + h) = x(t) + h k[0]
    for (x_grp = 0; x_grp < system->grps; x_grp++) {
        for (x_dim = 0; x_dim < system->dims[x_grp]; x_dim++) {
            prec_x = arpra_get_precision(&(system->x[x_grp][x_dim]));
            arpra_set_precision(&(scratch->temp_x), prec_x);
            arpra_mul(&(scratch->temp_x), h, &(scratch->k_0[x_grp][x_dim]));
            arpra_add(&(scratch->x_new[x_grp][x_dim]), &(system->x[x_grp][x_dim]), &(scratch->temp_x));
        }
    }

    // Advance system.
    arpra_add(system->t, system->t, h);
    for (x_grp = 0; x_grp < system->grps; x_grp++) {
        for (x_dim = 0; x_dim < system->dims[x_grp]; x_dim++) {
            arpra_set(&(system->x[x_grp][x_dim]), &(scratch->x_new[x_grp][x_dim]));
        }
    }
}

static const arpra_ode_method euler =
{
    .init = &euler_init,
    .clear = &euler_clear,
    .step = &euler_step,
    .stages = euler_stages,
};

const arpra_ode_method *arpra_ode_euler = &euler;
