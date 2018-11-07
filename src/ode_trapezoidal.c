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

#define trapezoidal_stages 2

typedef struct trapezoidal_scratch_struct
{
    arpra_range *_k_0;
    arpra_range **k_0;
    arpra_range *_k_1;
    arpra_range **k_1;
    arpra_range *_x_new;
    arpra_range **x_new;
    arpra_range half;
    arpra_range half_h;
    arpra_range temp_t;
    arpra_range temp_x;
} trapezoidal_scratch;

static void trapezoidal_init (arpra_ode_stepper *stepper, arpra_ode_system *system)
{
    arpra_uint x_grp, x_dim, state_size;
    arpra_precision prec_x, prec_internal;
    trapezoidal_scratch *scratch;

    // Allocate scratch memory.
    scratch = malloc(sizeof(trapezoidal_scratch));
    for (x_grp = 0, state_size = 0; x_grp < system->grps; x_grp++) {
        state_size += system->dims[x_grp];
    }
    scratch->_k_0 = malloc(state_size * sizeof(arpra_range));
    scratch->k_0 = malloc(system->grps * sizeof(arpra_range *));
    scratch->_k_1 = malloc(state_size * sizeof(arpra_range));
    scratch->k_1 = malloc(system->grps * sizeof(arpra_range *));
    scratch->_x_new = malloc(state_size * sizeof(arpra_range));
    scratch->x_new = malloc(system->grps * sizeof(arpra_range *));

    // Initialise scratch memory.
    prec_internal = arpra_get_internal_precision();
    scratch->k_0[0] = scratch->_k_0;
    scratch->k_1[0] = scratch->_k_1;
    scratch->x_new[0] = scratch->_x_new;
    for (x_grp = 1; x_grp < system->grps; x_grp++) {
        scratch->k_0[x_grp] = scratch->k_0[x_grp - 1] + system->dims[x_grp - 1];
        scratch->k_1[x_grp] = scratch->k_1[x_grp - 1] + system->dims[x_grp - 1];
        scratch->x_new[x_grp] = scratch->x_new[x_grp - 1] + system->dims[x_grp - 1];
    }
    for (x_grp = 0; x_grp < system->grps; x_grp++) {
        for (x_dim = 0; x_dim < system->dims[x_grp]; x_dim++) {
            prec_x = arpra_get_precision(&(system->x[x_grp][x_dim]));
            arpra_init2(&(scratch->k_0[x_grp][x_dim]), prec_x);
            arpra_init2(&(scratch->k_1[x_grp][x_dim]), prec_x);
            arpra_init2(&(scratch->x_new[x_grp][x_dim]), prec_x);
        }
    }
    arpra_init2(&(scratch->half), 2);
    arpra_init2(&(scratch->half_h), prec_internal);
    arpra_init2(&(scratch->temp_t), prec_internal);
    arpra_init2(&(scratch->temp_x), prec_internal);

    // Set stepper parameters.
    stepper->method = arpra_ode_trapezoidal;
    stepper->system = system;
    stepper->error = NULL;
    stepper->scratch = scratch;

    // Precompute constants.
    arpra_set_d(&(scratch->half), 0.5);
}

static void trapezoidal_clear (arpra_ode_stepper *stepper)
{
    arpra_uint x_grp, x_dim;
    arpra_ode_system *system;
    trapezoidal_scratch *scratch;

    system = stepper->system;
    scratch = (trapezoidal_scratch *) stepper->scratch;

    // Clear scratch memory.
    for (x_grp = 0; x_grp < system->grps; x_grp++) {
        for (x_dim = 0; x_dim < system->dims[x_grp]; x_dim++) {
            arpra_clear(&(scratch->k_0[x_grp][x_dim]));
            arpra_clear(&(scratch->k_1[x_grp][x_dim]));
            arpra_clear(&(scratch->x_new[x_grp][x_dim]));
        }
    }
    arpra_clear(&(scratch->half));
    arpra_clear(&(scratch->half_h));
    arpra_clear(&(scratch->temp_t));
    arpra_clear(&(scratch->temp_x));

    // Free scratch memory.
    free(scratch->_k_0);
    free(scratch->k_0);
    free(scratch->_k_1);
    free(scratch->k_1);
    free(scratch->_x_new);
    free(scratch->x_new);
    free(scratch);
}

static void trapezoidal_step (arpra_ode_stepper *stepper, const arpra_range *h)
{
    arpra_uint x_grp, x_dim;
    arpra_precision prec_t, prec_x;
    arpra_ode_system *system;
    trapezoidal_scratch *scratch;

    system = stepper->system;
    scratch = (trapezoidal_scratch *) stepper->scratch;

    // Synchronise scratch precision and prepare step parameters.
    prec_t = arpra_get_precision(system->t);
    for (x_grp = 0; x_grp < system->grps; x_grp++) {
        for (x_dim = 0; x_dim < system->dims[x_grp]; x_dim++) {
            prec_x = arpra_get_precision(&(system->x[x_grp][x_dim]));
            arpra_set_precision(&(scratch->k_0[x_grp][x_dim]), prec_x);
            arpra_set_precision(&(scratch->k_1[x_grp][x_dim]), prec_x);
            arpra_set_precision(&(scratch->x_new[x_grp][x_dim]), prec_x);
        }
    }
    arpra_set_precision(&(scratch->half_h), prec_t);
    arpra_mul(&(scratch->half_h), &(scratch->half), h);
    arpra_set_precision(&(scratch->temp_t), prec_t);
    arpra_add(&(scratch->temp_t), system->t, h);

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

    // k[1] = f(t + h, x(t) + h k[0])
    for (x_grp = 0; x_grp < system->grps; x_grp++) {
        for (x_dim = 0; x_dim < system->dims[x_grp]; x_dim++) {
            system->f[x_grp](&(scratch->k_1[x_grp][x_dim]), system->params[x_grp],
                             &(scratch->temp_t), (const arpra_range **) scratch->x_new, x_grp, x_dim);
        }
    }

    // x(t + h) = x(t) + 1/2 h k[0] + 1/2 h k[1]
    for (x_grp = 0; x_grp < system->grps; x_grp++) {
        for (x_dim = 0; x_dim < system->dims[x_grp]; x_dim++) {
            prec_x = arpra_get_precision(&(system->x[x_grp][x_dim]));
            arpra_set_precision(&(scratch->temp_x), prec_x);
            arpra_mul(&(scratch->temp_x), &(scratch->half_h), &(scratch->k_0[x_grp][x_dim]));
            arpra_add(&(scratch->x_new[x_grp][x_dim]), &(system->x[x_grp][x_dim]), &(scratch->temp_x));
            arpra_mul(&(scratch->temp_x), &(scratch->half_h), &(scratch->k_1[x_grp][x_dim]));
            arpra_add(&(scratch->x_new[x_grp][x_dim]), &(scratch->x_new[x_grp][x_dim]), &(scratch->temp_x));
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

static const arpra_ode_method trapezoidal =
{
    .init = &trapezoidal_init,
    .clear = &trapezoidal_clear,
    .step = &trapezoidal_step,
    .stages = trapezoidal_stages,
};

const arpra_ode_method *arpra_ode_trapezoidal = &trapezoidal;
