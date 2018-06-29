/*
 * ode_bogsham32.c -- Bogacki-Shampine 3(2) ODE stepper.
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

#define bogsham32_stages 4

typedef struct bogsham32_scratch_struct
{
    arpra_range a[bogsham32_stages][bogsham32_stages];
    arpra_range b_3[bogsham32_stages];
    arpra_range b_2[bogsham32_stages];
    arpra_range c[bogsham32_stages];


    arpra_range *k_weights_3;
    arpra_range *k_weights_2;


    arpra_range *k[bogsham32_stages];
    arpra_range *x_new_3;
    arpra_range *x_new_2;
    arpra_range temp;
} bogsham32_scratch;

static void bogsham32_init (arpra_ode_stepper *stepper, arpra_ode_system *system)
{
    arpra_uint i;
    arpra_precision prec;
    bogsham32_scratch *scratch;

    scratch = malloc(sizeof(bogsham32_scratch));
    scratch->k[0] = malloc(system->dims * sizeof(arpra_range));
    scratch->k[1] = malloc(system->dims * sizeof(arpra_range));
    scratch->k[2] = malloc(system->dims * sizeof(arpra_range));
    scratch->k[3] = malloc(system->dims * sizeof(arpra_range));
    scratch->k_weights_3 = malloc(bogsham32_stages * sizeof(arpra_range));
    scratch->k_weights_2 = malloc(bogsham32_stages * sizeof(arpra_range));
    scratch->x_new_3 = malloc(system->dims * sizeof(arpra_range));
    scratch->x_new_2 = malloc(system->dims * sizeof(arpra_range));
    for (i = 0; i < system->dims; i++) {
        prec = arpra_get_precision(&(system->x[i]));
        arpra_init2(&(scratch->k[0][i]), prec);
        arpra_init2(&(scratch->k[1][i]), prec);
        arpra_init2(&(scratch->k[2][i]), prec);
        arpra_init2(&(scratch->k[3][i]), prec);
        arpra_init2(&(scratch->x_new_3[i]), prec);
        arpra_init2(&(scratch->x_new_2[i]), prec);
    }
    prec = arpra_get_default_precision();
    for (i = 0; i < bogsham32_stages; i++) {
        arpra_init2(&(scratch->k_weights_3[i]), prec);
        arpra_init2(&(scratch->k_weights_2[i]), prec);
    }
    arpra_init2(&(scratch->temp), prec);
    stepper->method = arpra_ode_bogsham32;
    stepper->system = system;
    stepper->error = NULL;
    stepper->scratch = scratch;
}

static void bogsham32_clear (arpra_ode_stepper *stepper)
{
    arpra_uint i;
    arpra_ode_system *system;
    bogsham32_scratch *scratch;

    system = stepper->system;
    scratch = (bogsham32_scratch *) stepper->scratch;
    for (i = 0; i < system->dims; i++) {
        arpra_clear(&(scratch->k[0][i]));
        arpra_clear(&(scratch->k[1][i]));
        arpra_clear(&(scratch->k[2][i]));
        arpra_clear(&(scratch->k[3][i]));
        arpra_clear(&(scratch->x_new_3[i]));
        arpra_clear(&(scratch->x_new_2[i]));
    }
    for (i = 0; i < bogsham32_stages; i++) {
        arpra_clear(&(scratch->k_weights_3[i]));
        arpra_clear(&(scratch->k_weights_2[i]));
    }
    arpra_clear(&(scratch->temp));
    free(scratch->k[0]);
    free(scratch->k[1]);
    free(scratch->k[2]);
    free(scratch->k[3]);
    free(scratch->k_weights_3);
    free(scratch->k_weights_2);
    free(scratch->x_new_3);
    free(scratch->x_new_2);
    free(scratch);
}

static void bogsham32_step (arpra_ode_stepper *stepper, const arpra_range *h)
{
    arpra_uint i;
    arpra_precision prec;
    arpra_ode_system *system;
    bogsham32_scratch *scratch;

    system = stepper->system;
    scratch = (bogsham32_scratch *) stepper->scratch;

    // Synchronise scratch memory precision.
    for (i = 0; i < system->dims; i++) {
        prec = arpra_get_precision(&(system->x[i]));
        arpra_set_precision(&(scratch->k[0][i]), prec);
        arpra_set_precision(&(scratch->k[1][i]), prec);
        arpra_set_precision(&(scratch->k[2][i]), prec);
        arpra_set_precision(&(scratch->k[3][i]), prec);
        arpra_set_precision(&(scratch->x_new_3[i]), prec);
        arpra_set_precision(&(scratch->x_new_2[i]), prec);
    }
    prec = arpra_get_precision(system->t);
    for (i = 0; i < bogsham32_stages; i++) {
        arpra_set_precision(&(scratch->k_weights_3[i]), prec);
        arpra_set_precision(&(scratch->k_weights_2[i]), prec);
    }
    arpra_set_precision(&(scratch->temp), prec);

    // k[0] = f([t], [x(t)])
    system->f(scratch->k[0],
              system->t, system->x,
              system->dims, system->params);

    // k[1] = f([t + h/2], [x(t) + h/2 k[0]])
    arpra_set_d(&(scratch->temp), 2.0);
    arpra_div(&(scratch->k_weights_3[0]), h, &(scratch->temp));
    for (i = 0; i < system->dims; i++) {
        arpra_mul(&(scratch->x_new_3[i]), &(scratch->k_weights_3[0]), &(scratch->k[0][i]));
        arpra_add(&(scratch->x_new_3[i]), &(system->x[i]), &(scratch->x_new_3[i]));
    }
    arpra_add(&(scratch->temp), system->t, &(scratch->k_weights_3[0]));
    system->f(scratch->k[1],
              &(scratch->temp), scratch->x_new_3,
              system->dims, system->params);

    // k[2] = f([t + 3h/4], [x(t) + 3h/4 k[1]])
    arpra_set_d(&(scratch->temp), 4.0);
    arpra_div(&(scratch->k_weights_3[0]), h, &(scratch->temp));
    arpra_set_d(&(scratch->temp), 3.0);
    arpra_mul(&(scratch->k_weights_3[0]), &(scratch->temp), &(scratch->k_weights_3[0]));
    for (i = 0; i < system->dims; i++) {
        arpra_mul(&(scratch->x_new_3[i]), &(scratch->k_weights_3[0]), &(scratch->k[1][i]));
        arpra_add(&(scratch->x_new_3[i]), &(system->x[i]), &(scratch->x_new_3[i]));
    }
    arpra_add(&(scratch->temp), system->t, &(scratch->k_weights_3[0]));
    system->f(scratch->k[2],
              &(scratch->temp), scratch->x_new_3,
              system->dims, system->params);

    // k[3] = f([t + h], [x(t) + 2h/9 k[0] + 3h/9 k[1] + 4h/9 k[2]])
    arpra_set_d(&(scratch->temp), 9.0);
    arpra_div(&(scratch->k_weights_3[3]), h, &(scratch->temp));
    arpra_set_d(&(scratch->temp), 2.0);
    arpra_mul(&(scratch->k_weights_3[0]), &(scratch->temp), &(scratch->k_weights_3[3]));
    arpra_set_d(&(scratch->temp), 3.0);
    arpra_mul(&(scratch->k_weights_3[1]), &(scratch->temp), &(scratch->k_weights_3[3]));
    arpra_set_d(&(scratch->temp), 4.0);
    arpra_mul(&(scratch->k_weights_3[2]), &(scratch->temp), &(scratch->k_weights_3[3]));
    arpra_set_zero(&(scratch->k_weights_3[3]));
    for (i = 0; i < system->dims; i++) {
        prec = arpra_get_precision(&(system->x[i]));
        arpra_set_precision(&(scratch->temp), prec);
        arpra_mul(&(scratch->temp), &(scratch->k_weights_3[0]), &(scratch->k[0][i]));
        arpra_add(&(scratch->x_new_3[i]), &(system->x[i]), &(scratch->temp));
        arpra_mul(&(scratch->temp), &(scratch->k_weights_3[1]), &(scratch->k[1][i]));
        arpra_add(&(scratch->x_new_3[i]), &(scratch->x_new_3[i]), &(scratch->temp));
        arpra_mul(&(scratch->temp), &(scratch->k_weights_3[2]), &(scratch->k[2][i]));
        arpra_add(&(scratch->x_new_3[i]), &(scratch->x_new_3[i]), &(scratch->temp));
    }
    arpra_add(&(scratch->temp), system->t, h);
    system->f(scratch->k[3],
              &(scratch->temp), scratch->x_new_3,
              system->dims, system->params);

    // x_3(t + h) = x(t) + 2h/9 k[0] + 3h/9 k[1] + 4h/9 k[2]
    // This has already been computed above.

    // x_2(t + h) = x(t) + 7h/24 k[0] + 6h/24 k[1] + 8h/24 k[2] + 3h/24 k[3]
    arpra_set_d(&(scratch->temp), 24.0);
    arpra_div(&(scratch->k_weights_2[3]), h, &(scratch->temp));
    arpra_set_d(&(scratch->temp), 7.0);
    arpra_mul(&(scratch->k_weights_2[0]), &(scratch->temp), &(scratch->k_weights_2[3]));
    arpra_set_d(&(scratch->temp), 6.0);
    arpra_mul(&(scratch->k_weights_2[1]), &(scratch->temp), &(scratch->k_weights_2[3]));
    arpra_set_d(&(scratch->temp), 8.0);
    arpra_mul(&(scratch->k_weights_2[2]), &(scratch->temp), &(scratch->k_weights_2[3]));
    arpra_set_d(&(scratch->temp), 3.0);
    arpra_mul(&(scratch->k_weights_2[3]), &(scratch->temp), &(scratch->k_weights_2[3]));
    for (i = 0; i < system->dims; i++) {
        prec = arpra_get_precision(&(system->x[i]));
        arpra_set_precision(&(scratch->temp), prec);
        arpra_mul(&(scratch->temp), &(scratch->k_weights_2[0]), &(scratch->k[0][i]));
        arpra_add(&(scratch->x_new_2[i]), &(system->x[i]), &(scratch->temp));
        arpra_mul(&(scratch->temp), &(scratch->k_weights_2[1]), &(scratch->k[1][i]));
        arpra_add(&(scratch->x_new_2[i]), &(scratch->x_new_2[i]), &(scratch->temp));
        arpra_mul(&(scratch->temp), &(scratch->k_weights_2[2]), &(scratch->k[2][i]));
        arpra_add(&(scratch->x_new_2[i]), &(scratch->x_new_2[i]), &(scratch->temp));
        arpra_mul(&(scratch->temp), &(scratch->k_weights_2[3]), &(scratch->k[3][i]));
        arpra_add(&(scratch->x_new_2[i]), &(scratch->x_new_2[i]), &(scratch->temp));
    }

    // Advance system.
    arpra_add(system->t, system->t, h);
    arpra_range *x_temp = system->x;
    system->x = scratch->x_new_3;
    scratch->x_new_3 = x_temp;
}

static const arpra_ode_method bogsham32 =
{
    .init = &bogsham32_init,
    .clear = &bogsham32_clear,
    .step = &bogsham32_step,
    .stages = bogsham32_stages,
};

const arpra_ode_method *arpra_ode_bogsham32 = &bogsham32;
