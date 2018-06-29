/*
 * ode_dopri54.c -- Dormand-Prince 5(4) ODE stepper.
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

#define dopri54_stages 7

typedef struct dopri54_scratch_struct
{
    arpra_range *k[dopri54_stages];


    arpra_range a[dopri54_stages][dopri54_stages];
    arpra_range b_5[dopri54_stages];
    arpra_range b_4[dopri54_stages];
    arpra_range c[dopri54_stages];


    arpra_range *k_weights_5;
    arpra_range *k_weights_4;


    arpra_range *next_t;
    arpra_range *next_x_5;
    arpra_range *next_x_4;
    arpra_range temp;
} dopri54_scratch;

static void dopri54_init (arpra_ode_stepper *stepper, arpra_ode_system *system)
{
    arpra_uint i;
    arpra_precision prec;
    dopri54_scratch *scratch;

    scratch = malloc(sizeof(dopri54_scratch));
    scratch->k[0] = malloc(system->dims * sizeof(arpra_range));
    scratch->k[1] = malloc(system->dims * sizeof(arpra_range));
    scratch->k[2] = malloc(system->dims * sizeof(arpra_range));
    scratch->k[3] = malloc(system->dims * sizeof(arpra_range));
    scratch->k[4] = malloc(system->dims * sizeof(arpra_range));
    scratch->k[5] = malloc(system->dims * sizeof(arpra_range));
    scratch->k[6] = malloc(system->dims * sizeof(arpra_range));
    scratch->k_weights_5 = malloc(dopri54_stages * sizeof(arpra_range));
    scratch->k_weights_4 = malloc(dopri54_stages * sizeof(arpra_range));
    scratch->next_t = malloc(sizeof(arpra_range));
    scratch->next_x_5 = malloc(system->dims * sizeof(arpra_range));
    scratch->next_x_4 = malloc(system->dims * sizeof(arpra_range));
    for (i = 0; i < system->dims; i++) {
        prec = arpra_get_precision(&(system->x[i]));
        arpra_init2(&(scratch->k[0][i]), prec);
        arpra_init2(&(scratch->k[1][i]), prec);
        arpra_init2(&(scratch->k[2][i]), prec);
        arpra_init2(&(scratch->k[3][i]), prec);
        arpra_init2(&(scratch->k[4][i]), prec);
        arpra_init2(&(scratch->k[5][i]), prec);
        arpra_init2(&(scratch->k[6][i]), prec);
        arpra_init2(&(scratch->next_x_5[i]), prec);
        arpra_init2(&(scratch->next_x_4[i]), prec);
    }
    prec = arpra_get_default_precision();
    for (i = 0; i < dopri54_stages; i++) {
        arpra_init2(&(scratch->k_weights_5[i]), prec);
        arpra_init2(&(scratch->k_weights_4[i]), prec);
    }
    arpra_init2(scratch->next_t, prec);
    arpra_init2(&(scratch->temp), prec);
    stepper->method = arpra_ode_dopri54;
    stepper->system = system;
    stepper->error = NULL;
    stepper->scratch = scratch;
}

static void dopri54_clear (arpra_ode_stepper *stepper)
{
    arpra_uint i;
    arpra_ode_system *system;
    dopri54_scratch *scratch;

    system = stepper->system;
    scratch = (dopri54_scratch *) stepper->scratch;
    for (i = 0; i < system->dims; i++) {
        arpra_clear(&(scratch->k[0][i]));
        arpra_clear(&(scratch->k[1][i]));
        arpra_clear(&(scratch->k[2][i]));
        arpra_clear(&(scratch->k[3][i]));
        arpra_clear(&(scratch->k[4][i]));
        arpra_clear(&(scratch->k[5][i]));
        arpra_clear(&(scratch->k[6][i]));
        arpra_clear(&(scratch->next_x_5[i]));
        arpra_clear(&(scratch->next_x_4[i]));
    }
    for (i = 0; i < dopri54_stages; i++) {
        arpra_clear(&(scratch->k_weights_5[i]));
        arpra_clear(&(scratch->k_weights_4[i]));
    }
    arpra_clear(scratch->next_t);
    arpra_clear(&(scratch->temp));
    free(scratch->k[0]);
    free(scratch->k[1]);
    free(scratch->k[2]);
    free(scratch->k[3]);
    free(scratch->k[4]);
    free(scratch->k[5]);
    free(scratch->k[6]);
    free(scratch->k_weights_5);
    free(scratch->k_weights_4);
    free(scratch->next_t);
    free(scratch->next_x_5);
    free(scratch->next_x_4);
    free(scratch);
}

static void dopri54_step (arpra_ode_stepper *stepper, const arpra_range *h)
{
    arpra_uint i;
    arpra_precision prec;
    arpra_ode_system *system;
    dopri54_scratch *scratch;

    system = stepper->system;
    scratch = (dopri54_scratch *) stepper->scratch;

    // Synchronise scratch memory precision.
    for (i = 0; i < system->dims; i++) {
        prec = arpra_get_precision(&(system->x[i]));
        arpra_set_precision(&(scratch->k[0][i]), prec);
        arpra_set_precision(&(scratch->k[1][i]), prec);
        arpra_set_precision(&(scratch->k[2][i]), prec);
        arpra_set_precision(&(scratch->k[3][i]), prec);
        arpra_set_precision(&(scratch->k[4][i]), prec);
        arpra_set_precision(&(scratch->k[5][i]), prec);
        arpra_set_precision(&(scratch->k[6][i]), prec);
        arpra_set_precision(&(scratch->next_x_5[i]), prec);
        arpra_set_precision(&(scratch->next_x_4[i]), prec);
    }
    prec = arpra_get_precision(system->t);
    for (i = 0; i < dopri54_stages; i++) {
        arpra_set_precision(&(scratch->k_weights_5[i]), prec);
        arpra_set_precision(&(scratch->k_weights_4[i]), prec);
    }
    arpra_set_precision(scratch->next_t, prec);
    arpra_set_precision(&(scratch->temp), prec);

    // k[0] = f([t], [x(t)])
    system->f(scratch->k[0],
              system->t, system->x,
              system->dims, system->params);

    // k[1] = f([t + h/5], [x(t) + h/5 k[0]])
    arpra_set_d(&(scratch->temp), 2.0);
    arpra_div(&(scratch->k_weights_5[0]), h, &(scratch->temp));
    for (i = 0; i < system->dims; i++) {
        arpra_mul(&(scratch->next_x_5[i]), &(scratch->k_weights_5[0]), &(scratch->k[0][i]));
        arpra_add(&(scratch->next_x_5[i]), &(system->x[i]), &(scratch->next_x_5[i]));
    }
    arpra_add(scratch->next_t, system->t, &(scratch->k_weights_5[0]));
    system->f(scratch->k[1],
              scratch->next_t, scratch->next_x_5,
              system->dims, system->params);

    // k[2] = f([t + 12h/40], [x(t) + 3h/40 k[0] + 9h/40 k[1]])
    arpra_set_d(&(scratch->temp), 4.0);
    arpra_div(&(scratch->k_weights_5[0]), h, &(scratch->temp));
    arpra_set_d(&(scratch->temp), 3.0);
    arpra_mul(&(scratch->k_weights_5[0]), &(scratch->temp), &(scratch->k_weights_5[0]));
    for (i = 0; i < system->dims; i++) {
        arpra_mul(&(scratch->next_x_5[i]), &(scratch->k_weights_5[0]), &(scratch->k[1][i]));
        arpra_add(&(scratch->next_x_5[i]), &(system->x[i]), &(scratch->next_x_5[i]));
    }
    arpra_add(scratch->next_t, system->t, &(scratch->k_weights_5[0]));
    system->f(scratch->k[2],
              scratch->next_t, scratch->next_x_5,
              system->dims, system->params);

    // k[3] = f([t + 36h/45], [x(t) + 44h/45 k[0] - 168h/45 k[1] + 160h/45])
    arpra_set_d(&(scratch->temp), 9.0);
    arpra_div(&(scratch->k_weights_5[3]), h, &(scratch->temp));
    arpra_set_d(&(scratch->temp), 2.0);
    arpra_mul(&(scratch->k_weights_5[0]), &(scratch->temp), &(scratch->k_weights_5[3]));
    arpra_set_d(&(scratch->temp), 3.0);
    arpra_mul(&(scratch->k_weights_5[1]), &(scratch->temp), &(scratch->k_weights_5[3]));
    arpra_set_d(&(scratch->temp), 4.0);
    arpra_mul(&(scratch->k_weights_5[2]), &(scratch->temp), &(scratch->k_weights_5[3]));
    arpra_set_zero(&(scratch->k_weights_5[3]));
    for (i = 0; i < system->dims; i++) {
        prec = arpra_get_precision(&(system->x[i]));
        arpra_set_precision(&(scratch->temp), prec);
        arpra_mul(&(scratch->temp), &(scratch->k_weights_5[0]), &(scratch->k[0][i]));
        arpra_add(&(scratch->next_x_5[i]), &(system->x[i]), &(scratch->temp));
        arpra_mul(&(scratch->temp), &(scratch->k_weights_5[1]), &(scratch->k[1][i]));
        arpra_add(&(scratch->next_x_5[i]), &(scratch->next_x_5[i]), &(scratch->temp));
        arpra_mul(&(scratch->temp), &(scratch->k_weights_5[2]), &(scratch->k[2][i]));
        arpra_add(&(scratch->next_x_5[i]), &(scratch->next_x_5[i]), &(scratch->temp));
    }
    arpra_add(scratch->next_t, system->t, h);
    system->f(scratch->k[3],
              scratch->next_t, scratch->next_x_5,
              system->dims, system->params);

    // k[4] = f([t + 5832h/6561],
    //          [x(t) + 19372h/6561 k[0] - 76080h/6561 k[1]
    //                + 64448h/6561 k[2] - 1908h/6561 k[3]])


    // k[5] = f([t + h],
    //          [x(t) + ])


    // k[6] = f([t + h],
    //          [x(t) + ])


    // x_5(t + h) = x(t) +
    // This has already been computed above.

    // x_4(t + h) = x(t) +
    arpra_set_d(&(scratch->temp), 24.0);
    arpra_div(&(scratch->k_weights_4[3]), h, &(scratch->temp));
    arpra_set_d(&(scratch->temp), 7.0);
    arpra_mul(&(scratch->k_weights_4[0]), &(scratch->temp), &(scratch->k_weights_4[3]));
    arpra_set_d(&(scratch->temp), 6.0);
    arpra_mul(&(scratch->k_weights_4[1]), &(scratch->temp), &(scratch->k_weights_4[3]));
    arpra_set_d(&(scratch->temp), 8.0);
    arpra_mul(&(scratch->k_weights_4[2]), &(scratch->temp), &(scratch->k_weights_4[3]));
    arpra_set_d(&(scratch->temp), 3.0);
    arpra_mul(&(scratch->k_weights_4[3]), &(scratch->temp), &(scratch->k_weights_4[3]));
    for (i = 0; i < system->dims; i++) {
        prec = arpra_get_precision(&(system->x[i]));
        arpra_set_precision(&(scratch->temp), prec);
        arpra_mul(&(scratch->temp), &(scratch->k_weights_4[0]), &(scratch->k[0][i]));
        arpra_add(&(scratch->next_x_4[i]), &(system->x[i]), &(scratch->temp));
        arpra_mul(&(scratch->temp), &(scratch->k_weights_4[1]), &(scratch->k[1][i]));
        arpra_add(&(scratch->next_x_4[i]), &(scratch->next_x_4[i]), &(scratch->temp));
        arpra_mul(&(scratch->temp), &(scratch->k_weights_4[2]), &(scratch->k[2][i]));
        arpra_add(&(scratch->next_x_4[i]), &(scratch->next_x_4[i]), &(scratch->temp));
        arpra_mul(&(scratch->temp), &(scratch->k_weights_4[3]), &(scratch->k[3][i]));
        arpra_add(&(scratch->next_x_4[i]), &(scratch->next_x_4[i]), &(scratch->temp));
    }

    // Advance system.
    arpra_add(system->t, system->t, h);
    arpra_range *temp_x = system->x;
    system->x = scratch->next_x_5;
    scratch->next_x_5 = temp_x;
}

static const arpra_ode_method dopri54 =
{
    .init = &dopri54_init,
    .clear = &dopri54_clear,
    .step = &dopri54_step,
    .stages = dopri54_stages,
};

const arpra_ode_method *arpra_ode_dopri54 = &dopri54;
