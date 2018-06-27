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

typedef struct bogsham32_scratch_struct
{
    arpra_range *k_1;
    arpra_range *k_2;
    arpra_range *k_3;
    arpra_range *k_4;
    arpra_range *k_weights_3;
    arpra_range *k_weights_2;
    arpra_range *temp;
    arpra_range *temp_t;
    arpra_range *temp_x;
} bogsham32_scratch;

static void bogsham32_init (arpra_ode_stepper *stepper, arpra_ode_system *system)
{
    arpra_uint i;
    arpra_precision prec;
    bogsham32_scratch *scratch;

    scratch = malloc(sizeof(bogsham32_scratch));
    scratch->k_1 = malloc(system->dims * sizeof(arpra_range));
    scratch->k_2 = malloc(system->dims * sizeof(arpra_range));
    scratch->k_3 = malloc(system->dims * sizeof(arpra_range));
    scratch->k_4 = malloc(system->dims * sizeof(arpra_range));
    scratch->k_weights_3 = malloc(4 * sizeof(arpra_range));
    scratch->k_weights_2 = malloc(4 * sizeof(arpra_range));
    scratch->temp = malloc(sizeof(arpra_range));
    scratch->temp_t = malloc(sizeof(arpra_range));
    scratch->temp_x = malloc(system->dims * sizeof(arpra_range));
    for (i = 0; i < system->dims; i++) {
        prec = arpra_get_precision(&(system->x[i]));
        arpra_init2(&(scratch->k_1[i]), prec);
        arpra_init2(&(scratch->k_2[i]), prec);
        arpra_init2(&(scratch->k_3[i]), prec);
        arpra_init2(&(scratch->k_4[i]), prec);
        arpra_init2(&(scratch->temp_x[i]), prec);
    }
    prec = arpra_get_default_precision();
    for (i = 0; i < 4; i++) {
        arpra_init2(&(scratch->k_weights_3[i]), prec);
        arpra_init2(&(scratch->k_weights_2[i]), prec);
    }
    arpra_init2(scratch->temp, prec);
    arpra_init2(scratch->temp_t, prec);
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
        arpra_clear(&(scratch->k_1[i]));
        arpra_clear(&(scratch->k_2[i]));
        arpra_clear(&(scratch->k_3[i]));
        arpra_clear(&(scratch->k_4[i]));
        arpra_clear(&(scratch->temp_x[i]));
    }
    for (i = 0; i < 4; i++) {
        arpra_clear(&(scratch->k_weights_3[i]));
        arpra_clear(&(scratch->k_weights_2[i]));
    }
    arpra_clear(scratch->temp);
    arpra_clear(scratch->temp_t);
    free(scratch->k_1);
    free(scratch->k_2);
    free(scratch->k_3);
    free(scratch->k_4);
    free(scratch->k_weights_3);
    free(scratch->k_weights_2);
    free(scratch->temp);
    free(scratch->temp_t);
    free(scratch->temp_x);
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
        arpra_set_precision(&(scratch->k_1[i]), prec);
        arpra_set_precision(&(scratch->k_2[i]), prec);
        arpra_set_precision(&(scratch->k_3[i]), prec);
        arpra_set_precision(&(scratch->k_4[i]), prec);
        arpra_set_precision(&(scratch->temp_x[i]), prec);
    }
    prec = arpra_get_precision(system->t);
    for (i = 0; i < 4; i++) {
        arpra_set_precision(&(scratch->k_weights_3[i]), prec);
        arpra_set_precision(&(scratch->k_weights_2[i]), prec);
    }
    arpra_set_precision(scratch->temp, prec);
    arpra_set_precision(scratch->temp_t, prec);

    // k_1 = f([t], [x(t)])
    system->f(scratch->k_1,
              system->t, system->x,
              system->dims, system->params);

    // k_2 = f([t + h/2], [x(t) + h/2 k_1])
    arpra_set_d(scratch->temp, 2.0);
    arpra_div(&(scratch->k_weights_3[0]), h, scratch->temp);
    for (i = 0; i < system->dims; i++) {
        arpra_mul(&(scratch->temp_x[i]), &(scratch->k_weights_3[0]), &(scratch->k_1[i]));
        arpra_add(&(scratch->temp_x[i]), &(system->x[i]), &(scratch->temp_x[i]));
    }
    arpra_add(scratch->temp_t, system->t, &(scratch->k_weights_3[0]));
    system->f(scratch->k_2,
              scratch->temp_t, scratch->temp_x,
              system->dims, system->params);

    // k_3 = f([t + 3h/4], [x(t) + 3h/4 k_2])
    arpra_set_d(scratch->temp, 4.0);
    arpra_div(&(scratch->k_weights_3[0]), h, scratch->temp);
    arpra_set_d(scratch->temp, 3.0);
    arpra_mul(&(scratch->k_weights_3[0]), scratch->temp, &(scratch->k_weights_3[0]));
    for (i = 0; i < system->dims; i++) {
        arpra_mul(&(scratch->temp_x[i]), &(scratch->k_weights_3[0]), &(scratch->k_2[i]));
        arpra_add(&(scratch->temp_x[i]), &(system->x[i]), &(scratch->temp_x[i]));
    }
    arpra_add(scratch->temp_t, system->t, &(scratch->k_weights_3[0]));
    system->f(scratch->k_3,
              scratch->temp_t, scratch->temp_x,
              system->dims, system->params);

    // k_4 = f([t + h], [x(t) + 2h/9 k_1 + 3h/9 k_2 + 4h/9 k_3])
    arpra_set_d(scratch->temp, 9.0);
    arpra_div(&(scratch->k_weights_3[2]), h, scratch->temp);
    arpra_set_d(scratch->temp, 2.0);
    arpra_mul(&(scratch->k_weights_3[0]), scratch->temp, &(scratch->k_weights_3[2]));
    arpra_set_d(scratch->temp, 3.0);
    arpra_mul(&(scratch->k_weights_3[1]), scratch->temp, &(scratch->k_weights_3[2]));
    arpra_set_d(scratch->temp, 4.0);
    arpra_mul(&(scratch->k_weights_3[2]), scratch->temp, &(scratch->k_weights_3[2]));
    for (i = 0; i < system->dims; i++) {
        prec = arpra_get_precision(&(system->x[i]));
        arpra_set_precision(scratch->temp, prec);
        arpra_mul(scratch->temp, &(scratch->k_weights_3[0]), &(scratch->k_1[i]));
        arpra_add(&(scratch->temp_x[i]), &(system->x[i]), scratch->temp);
        arpra_mul(scratch->temp, &(scratch->k_weights_3[1]), &(scratch->k_2[i]));
        arpra_add(&(scratch->temp_x[i]), &(scratch->temp_x[i]), scratch->temp);
        arpra_mul(scratch->temp, &(scratch->k_weights_3[2]), &(scratch->k_3[i]));
        arpra_add(&(scratch->temp_x[i]), &(scratch->temp_x[i]), scratch->temp);
    }
    arpra_add(scratch->temp_t, system->t, h);
    system->f(scratch->k_4,
              scratch->temp_t, scratch->temp_x,
              system->dims, system->params);

    // x(t + h) = x(t) + 2h/9 k_1 + 3h/9 k_2 + 4h/9 k_3

    // THIS IS ALREADY IN TEMP_X HERE

    for (i = 0; i < system->dims; i++) {
        prec = arpra_get_precision(&(system->x[i]));
        arpra_set_precision(scratch->temp, prec);


        //arpra_add(&(system->x[i]), &(system->x[i]), &(scratch->temp_x[i]));
    }

    // e(t + h) = x(t) + 7h/24 k_1 + 6h/24 k_2 + 8h/24 k_3 + 3h/24 k_4

    for (i = 0; i < system->dims; i++) {
        prec = arpra_get_precision(&(system->x[i]));
        arpra_set_precision(scratch->temp, prec);


        //arpra_add(&(system->x[i]), &(system->x[i]), &(scratch->temp_x[i]));
    }

    // Advance t.
    arpra_add(system->t, system->t, h);
}

static const arpra_ode_method bogsham32 =
{
    .init = &bogsham32_init,
    .clear = &bogsham32_clear,
    .step = &bogsham32_step,
    .stages = 4,
};

const arpra_ode_method *arpra_ode_bogsham32 = &bogsham32;
