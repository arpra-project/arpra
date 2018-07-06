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
    arpra_range *k[bogsham32_stages];
    arpra_range *x_new_3;
    arpra_range *x_new_2;
    arpra_range temp;
    arpra_range weighted_h;
    arpra_range a[bogsham32_stages][bogsham32_stages];
    arpra_range b_3[bogsham32_stages];
    arpra_range b_2[bogsham32_stages];
    arpra_range c[bogsham32_stages];
} bogsham32_scratch;

static void bogsham32_compute_constants (arpra_ode_stepper *stepper)
{
    arpra_uint i, j;
    arpra_precision prec_internal;
    arpra_range numerator, denominator;
    bogsham32_scratch *scratch;

    scratch = (bogsham32_scratch *) stepper->scratch;

    // Update constant memory to internal precision.
    prec_internal = arpra_get_internal_precision();
    arpra_init2(&numerator, prec_internal);
    arpra_init2(&denominator, prec_internal);
    for (i = 0; i < bogsham32_stages; i++) {
        for (j = 0; j < bogsham32_stages; j++) {
            arpra_set_precision(&(scratch->a[i][j]), prec_internal);
        }
        arpra_set_precision(&(scratch->b_3[i]), prec_internal);
        arpra_set_precision(&(scratch->b_2[i]), prec_internal);
        arpra_set_precision(&(scratch->c[i]), prec_internal);
    }

    arpra_set_zero(&(scratch->c[0]));
    arpra_set_zero(&(scratch->a[0][0]));
    arpra_set_zero(&(scratch->a[0][1]));
    arpra_set_zero(&(scratch->a[0][2]));
    arpra_set_zero(&(scratch->a[0][3]));

    arpra_set_d(&numerator, 1.);
    arpra_set_d(&denominator, 2.);
    arpra_div(&(scratch->c[1]), &numerator, &denominator);
    arpra_set(&(scratch->a[1][0]), &(scratch->c[1]));
    arpra_set_zero(&(scratch->a[1][1]));
    arpra_set_zero(&(scratch->a[1][2]));
    arpra_set_zero(&(scratch->a[1][3]));

    arpra_set_d(&numerator, 3.);
    arpra_set_d(&denominator, 4.);
    arpra_div(&(scratch->c[2]), &numerator, &denominator);
    arpra_set_zero(&(scratch->a[2][0]));
    arpra_set(&(scratch->a[2][1]), &(scratch->c[2]));
    arpra_set_zero(&(scratch->a[2][2]));
    arpra_set_zero(&(scratch->a[2][3]));

    arpra_set_d(&(scratch->c[3]), 1.);
    arpra_set_d(&numerator, 2.);
    arpra_set_d(&denominator, 9.);
    arpra_div(&(scratch->a[3][0]), &numerator, &denominator);
    arpra_set(&(scratch->b_3[0]), &(scratch->a[3][0]));
    arpra_set_d(&numerator, 1.);
    arpra_set_d(&denominator, 3.);
    arpra_div(&(scratch->a[3][1]), &numerator, &denominator);
    arpra_set(&(scratch->b_3[1]), &(scratch->a[3][1]));
    arpra_set_d(&numerator, 4.);
    arpra_set_d(&denominator, 9.);
    arpra_div(&(scratch->a[3][2]), &numerator, &denominator);
    arpra_set(&(scratch->b_3[2]), &(scratch->a[3][2]));
    arpra_set_zero(&(scratch->a[3][3]));
    arpra_set_zero(&(scratch->b_3[3]));

    arpra_set_d(&numerator, 7.);
    arpra_set_d(&denominator, 24.);
    arpra_div(&(scratch->b_2[0]), &numerator, &denominator);
    arpra_set_d(&numerator, 1.);
    arpra_set_d(&denominator, 4.);
    arpra_div(&(scratch->b_2[1]), &numerator, &denominator);
    arpra_set_d(&numerator, 1.);
    arpra_set_d(&denominator, 3.);
    arpra_div(&(scratch->b_2[2]), &numerator, &denominator);
    arpra_set_d(&numerator, 1.);
    arpra_set_d(&denominator, 8.);
    arpra_div(&(scratch->b_2[3]), &numerator, &denominator);

    arpra_clear(&numerator);
    arpra_clear(&denominator);
}

static void bogsham32_init (arpra_ode_stepper *stepper, arpra_ode_system *system)
{
    arpra_uint i, j;
    arpra_precision prec_x, prec_internal;
    bogsham32_scratch *scratch;

    // Allocate scratch memory.
    scratch = malloc(sizeof(bogsham32_scratch));
    for (j = 0; j < bogsham32_stages; j++) {
        scratch->k[j] = malloc(system->dims * sizeof(arpra_range));
    }
    scratch->x_new_3 = malloc(system->dims * sizeof(arpra_range));
    scratch->x_new_2 = malloc(system->dims * sizeof(arpra_range));

    // Initialise scratch memory.
    for (i = 0; i < system->dims; i++) {
        prec_x = arpra_get_precision(&(system->x[i]));
        for (j = 0; j < bogsham32_stages; j++) {
            arpra_init2(&(scratch->k[j][i]), prec_x);
        }
        arpra_init2(&(scratch->x_new_3[i]), prec_x);
        arpra_init2(&(scratch->x_new_2[i]), prec_x);
    }
    prec_internal = arpra_get_internal_precision();
    arpra_init2(&(scratch->temp), prec_internal);
    arpra_init2(&(scratch->weighted_h), prec_internal);
    for (i = 0; i < bogsham32_stages; i++) {
        for (j = 0; j < bogsham32_stages; j++) {
            arpra_init2(&(scratch->a[i][j]), prec_internal);
        }
        arpra_init2(&(scratch->b_3[i]), prec_internal);
        arpra_init2(&(scratch->b_2[i]), prec_internal);
        arpra_init2(&(scratch->c[i]), prec_internal);
    }

    // Precompute constants.
    bogsham32_compute_constants(stepper);

    // Set stepper parameters.
    stepper->method = arpra_ode_bogsham32;
    stepper->system = system;
    stepper->error = NULL;
    stepper->scratch = scratch;
}

static void bogsham32_clear (arpra_ode_stepper *stepper)
{
    arpra_uint i, j;
    arpra_ode_system *system;
    bogsham32_scratch *scratch;

    system = stepper->system;
    scratch = (bogsham32_scratch *) stepper->scratch;

    // Clear scratch memory.
    for (i = 0; i < system->dims; i++) {
        for (j = 0; j < bogsham32_stages; j++) {
            arpra_clear(&(scratch->k[j][i]));
        }
        arpra_clear(&(scratch->x_new_3[i]));
        arpra_clear(&(scratch->x_new_2[i]));
    }
    arpra_clear(&(scratch->temp));
    arpra_clear(&(scratch->weighted_h));
    for (i = 0; i < bogsham32_stages; i++) {
        for (j = 0; j < bogsham32_stages; j++) {
            arpra_clear(&(scratch->a[i][j]));
        }
        arpra_clear(&(scratch->b_3[i]));
        arpra_clear(&(scratch->b_2[i]));
        arpra_clear(&(scratch->c[i]));
    }

    // Free scratch memory.
    for (j = 0; j < bogsham32_stages; j++) {
        free(scratch->k[j]);
    }
    free(scratch->x_new_3);
    free(scratch->x_new_2);
    free(scratch);
}

static void bogsham32_step (arpra_ode_stepper *stepper, const arpra_range *h)
{
    arpra_uint i, j;
    arpra_precision prec_t, prec_x;
    arpra_ode_system *system;
    bogsham32_scratch *scratch;

    system = stepper->system;
    scratch = (bogsham32_scratch *) stepper->scratch;

    // Synchronise scratch memory precision.
    for (i = 0; i < system->dims; i++) {
        prec_x = arpra_get_precision(&(system->x[i]));
        for (j = 0; j < bogsham32_stages; j++) {
            arpra_set_precision(&(scratch->k[j][i]), prec_x);
        }
        arpra_set_precision(&(scratch->x_new_3[i]), prec_x);
        arpra_set_precision(&(scratch->x_new_2[i]), prec_x);
    }
    prec_t = arpra_get_precision(system->t);
    arpra_set_precision(&(scratch->weighted_h), prec_t);

    // k[0] = f(t, x(t))
    system->f(scratch->k[0],
              system->t, system->x,
              system->dims, system->params);

    // k[1] = f(t + h/2, x(t) + h/2 k[0])
    arpra_mul(&(scratch->weighted_h), &(scratch->c[1]), h);
    for (i = 0; i < system->dims; i++) {
        prec_x = arpra_get_precision(&(system->x[i]));
        arpra_set_precision(&(scratch->temp), prec_x);
        arpra_mul(&(scratch->temp), &(scratch->weighted_h), &(scratch->k[0][i]));
        arpra_add(&(scratch->x_new_3[i]), &(system->x[i]), &(scratch->temp));
    }
    arpra_set_precision(&(scratch->temp), prec_t);
    arpra_add(&(scratch->temp), system->t, &(scratch->weighted_h));
    system->f(scratch->k[1],
              &(scratch->temp), scratch->x_new_3,
              system->dims, system->params);

    // k[2] = f(t + 3h/4, x(t) + 3h/4 k[1])
    arpra_mul(&(scratch->weighted_h), &(scratch->c[2]), h);
    for (i = 0; i < system->dims; i++) {
        prec_x = arpra_get_precision(&(system->x[i]));
        arpra_set_precision(&(scratch->temp), prec_x);
        arpra_mul(&(scratch->temp), &(scratch->weighted_h), &(scratch->k[1][i]));
        arpra_add(&(scratch->x_new_3[i]), &(system->x[i]), &(scratch->temp));
    }
    arpra_set_precision(&(scratch->temp), prec_t);
    arpra_add(&(scratch->temp), system->t, &(scratch->weighted_h));
    system->f(scratch->k[2],
              &(scratch->temp), scratch->x_new_3,
              system->dims, system->params);

    // k[3] = f(t + h, x(t) + 2h/9 k[0] + h/3 k[1] + 4h/9 k[2])
    for (i = 0; i < system->dims; i++) {
        prec_x = arpra_get_precision(&(system->x[i]));
        arpra_set_precision(&(scratch->temp), prec_x);
        arpra_mul(&(scratch->weighted_h), &(scratch->a[3][0]), h);
        arpra_mul(&(scratch->temp), &(scratch->weighted_h), &(scratch->k[0][i]));
        arpra_add(&(scratch->x_new_3[i]), &(system->x[i]), &(scratch->temp));
        arpra_mul(&(scratch->weighted_h), &(scratch->a[3][1]), h);
        arpra_mul(&(scratch->temp), &(scratch->weighted_h), &(scratch->k[1][i]));
        arpra_add(&(scratch->x_new_3[i]), &(scratch->x_new_3[i]), &(scratch->temp));
        arpra_mul(&(scratch->weighted_h), &(scratch->a[3][2]), h);
        arpra_mul(&(scratch->temp), &(scratch->weighted_h), &(scratch->k[2][i]));
        arpra_add(&(scratch->x_new_3[i]), &(scratch->x_new_3[i]), &(scratch->temp));
    }
    arpra_set_precision(&(scratch->temp), prec_t);
    arpra_add(&(scratch->temp), system->t, h);
    system->f(scratch->k[3],
              &(scratch->temp), scratch->x_new_3,
              system->dims, system->params);

    // Already been computed in x_new_3.
    // x_3(t + h) = x(t) + 2h/9 k[0] + h/3 k[1] + 4h/9 k[2]

    // x_2(t + h) = x(t) + 7h/24 k[0] + h/4 k[1] + h/3 k[2] + h/8 k[3]
    for (i = 0; i < system->dims; i++) {
        prec_x = arpra_get_precision(&(system->x[i]));
        arpra_set_precision(&(scratch->temp), prec_x);
        arpra_mul(&(scratch->weighted_h), &(scratch->b_2[0]), h);
        arpra_mul(&(scratch->temp), &(scratch->weighted_h), &(scratch->k[0][i]));
        arpra_add(&(scratch->x_new_2[i]), &(system->x[i]), &(scratch->temp));
        arpra_mul(&(scratch->weighted_h), &(scratch->b_2[1]), h);
        arpra_mul(&(scratch->temp), &(scratch->weighted_h), &(scratch->k[1][i]));
        arpra_add(&(scratch->x_new_2[i]), &(scratch->x_new_2[i]), &(scratch->temp));
        arpra_mul(&(scratch->weighted_h), &(scratch->b_2[2]), h);
        arpra_mul(&(scratch->temp), &(scratch->weighted_h), &(scratch->k[2][i]));
        arpra_add(&(scratch->x_new_2[i]), &(scratch->x_new_2[i]), &(scratch->temp));
        arpra_mul(&(scratch->weighted_h), &(scratch->b_2[3]), h);
        arpra_mul(&(scratch->temp), &(scratch->weighted_h), &(scratch->k[3][i]));
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
