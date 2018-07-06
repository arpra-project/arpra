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
    arpra_range *x_new_5;
    arpra_range *x_new_4;
    arpra_range temp;
    arpra_range weighted_h;
    arpra_range a[dopri54_stages][dopri54_stages];
    arpra_range b_5[dopri54_stages];
    arpra_range b_4[dopri54_stages];
    arpra_range c[dopri54_stages];
} dopri54_scratch;

static void dopri54_compute_constants (arpra_ode_stepper *stepper)
{
    arpra_uint i, j;
    arpra_precision prec_internal;
    arpra_range numerator, denominator;
    dopri54_scratch *scratch;

    scratch = (dopri54_scratch *) stepper->scratch;

    // Update constant memory to internal precision.
    prec_internal = arpra_get_internal_precision();
    arpra_init2(&numerator, prec_internal);
    arpra_init2(&denominator, prec_internal);
    for (i = 0; i < dopri54_stages; i++) {
        for (j = 0; j < dopri54_stages; j++) {
            arpra_set_precision(&(scratch->a[i][j]), prec_internal);
        }
        arpra_set_precision(&(scratch->b_5[i]), prec_internal);
        arpra_set_precision(&(scratch->b_4[i]), prec_internal);
        arpra_set_precision(&(scratch->c[i]), prec_internal);
    }

    arpra_set_zero(&(scratch->c[0]));
    arpra_set_zero(&(scratch->a[0][0]));
    arpra_set_zero(&(scratch->a[0][1]));
    arpra_set_zero(&(scratch->a[0][2]));
    arpra_set_zero(&(scratch->a[0][3]));
    arpra_set_zero(&(scratch->a[0][4]));
    arpra_set_zero(&(scratch->a[0][5]));
    arpra_set_zero(&(scratch->a[0][6]));

    arpra_set_d(&numerator, 1.);
    arpra_set_d(&denominator, 5.);
    arpra_div(&(scratch->c[1]), &numerator, &denominator);
    arpra_set(&(scratch->a[1][0]), &(scratch->c[1]));
    arpra_set_zero(&(scratch->a[1][1]));
    arpra_set_zero(&(scratch->a[1][2]));
    arpra_set_zero(&(scratch->a[1][3]));
    arpra_set_zero(&(scratch->a[1][4]));
    arpra_set_zero(&(scratch->a[1][5]));
    arpra_set_zero(&(scratch->a[1][6]));

    arpra_set_d(&numerator, 3.);
    arpra_set_d(&denominator, 10.);
    arpra_div(&(scratch->c[2]), &numerator, &denominator);
    arpra_set_d(&numerator, 3.);
    arpra_set_d(&denominator, 40.);
    arpra_div(&(scratch->a[2][0]), &numerator, &denominator);
    arpra_set_d(&numerator, 9.);
    arpra_set_d(&denominator, 40.);
    arpra_div(&(scratch->a[2][1]), &numerator, &denominator);
    arpra_set_zero(&(scratch->a[2][2]));
    arpra_set_zero(&(scratch->a[2][3]));
    arpra_set_zero(&(scratch->a[2][4]));
    arpra_set_zero(&(scratch->a[2][5]));
    arpra_set_zero(&(scratch->a[2][6]));

    arpra_set_d(&numerator, 4.);
    arpra_set_d(&denominator, 5.);
    arpra_div(&(scratch->c[3]), &numerator, &denominator);
    arpra_set_d(&numerator, 44.);
    arpra_set_d(&denominator, 45.);
    arpra_div(&(scratch->a[3][0]), &numerator, &denominator);
    arpra_set_d(&numerator, -56.);
    arpra_set_d(&denominator, 15.);
    arpra_div(&(scratch->a[3][1]), &numerator, &denominator);
    arpra_set_d(&numerator, 32.);
    arpra_set_d(&denominator, 9.);
    arpra_div(&(scratch->a[3][2]), &numerator, &denominator);
    arpra_set_zero(&(scratch->a[3][3]));
    arpra_set_zero(&(scratch->a[3][4]));
    arpra_set_zero(&(scratch->a[3][5]));
    arpra_set_zero(&(scratch->a[3][6]));

    arpra_set_d(&numerator, 8.);
    arpra_set_d(&denominator, 9.);
    arpra_div(&(scratch->c[4]), &numerator, &denominator);
    arpra_set_d(&numerator, 19372.);
    arpra_set_d(&denominator, 6561.);
    arpra_div(&(scratch->a[4][0]), &numerator, &denominator);
    arpra_set_d(&numerator, -25360.);
    arpra_set_d(&denominator, 2187.);
    arpra_div(&(scratch->a[4][1]), &numerator, &denominator);
    arpra_set_d(&numerator, 64448.);
    arpra_set_d(&denominator, 6561.);
    arpra_div(&(scratch->a[4][2]), &numerator, &denominator);
    arpra_set_d(&numerator, -212.);
    arpra_set_d(&denominator, 729.);
    arpra_div(&(scratch->a[4][3]), &numerator, &denominator);
    arpra_set_zero(&(scratch->a[4][4]));
    arpra_set_zero(&(scratch->a[4][5]));
    arpra_set_zero(&(scratch->a[4][6]));

    arpra_set_d(&(scratch->c[5]), 1.);
    arpra_set_d(&numerator, 9017.);
    arpra_set_d(&denominator, 3168.);
    arpra_div(&(scratch->a[5][0]), &numerator, &denominator);
    arpra_set_d(&numerator, -355.);
    arpra_set_d(&denominator, 33.);
    arpra_div(&(scratch->a[5][1]), &numerator, &denominator);
    arpra_set_d(&numerator, 46732.);
    arpra_set_d(&denominator, 5247.);
    arpra_div(&(scratch->a[5][2]), &numerator, &denominator);
    arpra_set_d(&numerator, 49.);
    arpra_set_d(&denominator, 176.);
    arpra_div(&(scratch->a[5][3]), &numerator, &denominator);
    arpra_set_d(&numerator, -5103.);
    arpra_set_d(&denominator, 18656.);
    arpra_div(&(scratch->a[5][4]), &numerator, &denominator);
    arpra_set_zero(&(scratch->a[5][5]));
    arpra_set_zero(&(scratch->a[5][6]));

    arpra_set_d(&(scratch->c[6]), 1.);
    arpra_set_d(&numerator, 35.);
    arpra_set_d(&denominator, 384.);
    arpra_div(&(scratch->a[6][0]), &numerator, &denominator);
    arpra_set(&(scratch->b_5[0]), &(scratch->a[6][0]));
    arpra_set_zero(&(scratch->a[6][1]));
    arpra_set_zero(&(scratch->b_5[1]));
    arpra_set_d(&numerator, 500.);
    arpra_set_d(&denominator, 1113.);
    arpra_div(&(scratch->a[6][2]), &numerator, &denominator);
    arpra_set(&(scratch->b_5[2]), &(scratch->a[6][2]));
    arpra_set_d(&numerator, 125.);
    arpra_set_d(&denominator, 192.);
    arpra_div(&(scratch->a[6][3]), &numerator, &denominator);
    arpra_set(&(scratch->b_5[3]), &(scratch->a[6][3]));
    arpra_set_d(&numerator, -2187.);
    arpra_set_d(&denominator, 6784.);
    arpra_div(&(scratch->a[6][4]), &numerator, &denominator);
    arpra_set(&(scratch->b_5[4]), &(scratch->a[6][4]));
    arpra_set_d(&numerator, 11.);
    arpra_set_d(&denominator, 84.);
    arpra_div(&(scratch->a[6][5]), &numerator, &denominator);
    arpra_set(&(scratch->b_5[5]), &(scratch->a[6][5]));
    arpra_set_zero(&(scratch->a[6][6]));
    arpra_set_zero(&(scratch->b_5[6]));

    arpra_set_d(&numerator, 5179.);
    arpra_set_d(&denominator, 57600.);
    arpra_div(&(scratch->b_4[0]), &numerator, &denominator);
    arpra_set_zero(&(scratch->b_4[1]));
    arpra_set_d(&numerator, 7571.);
    arpra_set_d(&denominator, 16695.);
    arpra_div(&(scratch->b_4[2]), &numerator, &denominator);
    arpra_set_d(&numerator, 393.);
    arpra_set_d(&denominator, 640.);
    arpra_div(&(scratch->b_4[3]), &numerator, &denominator);
    arpra_set_d(&numerator, -92097.);
    arpra_set_d(&denominator, 339200.);
    arpra_div(&(scratch->b_4[4]), &numerator, &denominator);
    arpra_set_d(&numerator, 187.);
    arpra_set_d(&denominator, 2100.);
    arpra_div(&(scratch->b_4[5]), &numerator, &denominator);
    arpra_set_d(&numerator, 1.);
    arpra_set_d(&denominator, 40.);
    arpra_div(&(scratch->b_4[6]), &numerator, &denominator);

    arpra_clear(&numerator);
    arpra_clear(&denominator);
}

static void dopri54_init (arpra_ode_stepper *stepper, arpra_ode_system *system)
{
    arpra_uint i, j;
    arpra_precision prec_x, prec_internal;
    dopri54_scratch *scratch;

    // Allocate scratch memory.
    scratch = malloc(sizeof(dopri54_scratch));
    for (j = 0; j < dopri54_stages; j++) {
        scratch->k[j] = malloc(system->dims * sizeof(arpra_range));
    }
    scratch->x_new_5 = malloc(system->dims * sizeof(arpra_range));
    scratch->x_new_4 = malloc(system->dims * sizeof(arpra_range));

    // Initialise scratch memory.
    for (i = 0; i < system->dims; i++) {
        prec_x = arpra_get_precision(&(system->x[i]));
        for (j = 0; j < dopri54_stages; j++) {
            arpra_init2(&(scratch->k[j][i]), prec_x);
        }
        arpra_init2(&(scratch->x_new_5[i]), prec_x);
        arpra_init2(&(scratch->x_new_4[i]), prec_x);
    }
    prec_internal = arpra_get_internal_precision();
    arpra_init2(&(scratch->temp), prec_internal);
    arpra_init2(&(scratch->weighted_h), prec_internal);
    for (i = 0; i < dopri54_stages; i++) {
        for (j = 0; j < dopri54_stages; j++) {
            arpra_init2(&(scratch->a[i][j]), prec_internal);
        }
        arpra_init2(&(scratch->b_5[i]), prec_internal);
        arpra_init2(&(scratch->b_4[i]), prec_internal);
        arpra_init2(&(scratch->c[i]), prec_internal);
    }

    // Precompute constants.
    dopri54_compute_constants(stepper);

    // Set stepper parameters.
    stepper->method = arpra_ode_dopri54;
    stepper->system = system;
    stepper->error = NULL;
    stepper->scratch = scratch;
}

static void dopri54_clear (arpra_ode_stepper *stepper)
{
    arpra_uint i, j;
    arpra_ode_system *system;
    dopri54_scratch *scratch;

    system = stepper->system;
    scratch = (dopri54_scratch *) stepper->scratch;

    // Clear scratch memory.
    for (i = 0; i < system->dims; i++) {
        for (j = 0; j < dopri54_stages; j++) {
            arpra_clear(&(scratch->k[j][i]));
        }
        arpra_clear(&(scratch->x_new_5[i]));
        arpra_clear(&(scratch->x_new_4[i]));
    }
    arpra_clear(&(scratch->temp));
    arpra_clear(&(scratch->weighted_h));
    for (i = 0; i < dopri54_stages; i++) {
        for (j = 0; j < dopri54_stages; j++) {
            arpra_clear(&(scratch->a[i][j]));
        }
        arpra_clear(&(scratch->b_5[i]));
        arpra_clear(&(scratch->b_4[i]));
        arpra_clear(&(scratch->c[i]));
    }

    // Free scratch memory.
    for (j = 0; j < dopri54_stages; j++) {
        free(scratch->k[j]);
    }
    free(scratch->x_new_5);
    free(scratch->x_new_4);
    free(scratch);
}

static void dopri54_step (arpra_ode_stepper *stepper, const arpra_range *h)
{
    arpra_uint i, j;
    arpra_precision prec_t, prec_x;
    arpra_ode_system *system;
    dopri54_scratch *scratch;

    system = stepper->system;
    scratch = (dopri54_scratch *) stepper->scratch;

    // Synchronise scratch memory precision.
    for (i = 0; i < system->dims; i++) {
        prec_x = arpra_get_precision(&(system->x[i]));
        for (j = 0; j < dopri54_stages; j++) {
            arpra_set_precision(&(scratch->k[j][i]), prec_x);
        }
        arpra_set_precision(&(scratch->x_new_5[i]), prec_x);
        arpra_set_precision(&(scratch->x_new_4[i]), prec_x);
    }
    prec_t = arpra_get_precision(system->t);
    arpra_set_precision(&(scratch->weighted_h), prec_t);

    // k[0] = f(t, x(t))
    system->f(scratch->k[0],
              system->t, system->x,
              system->dims, system->params);

    // k[1] = f(t + h/5, x(t) + h/5 k[0])
    arpra_mul(&(scratch->weighted_h), &(scratch->c[1]), h);
    for (i = 0; i < system->dims; i++) {
        prec_x = arpra_get_precision(&(system->x[i]));
        arpra_set_precision(&(scratch->temp), prec_x);
        arpra_mul(&(scratch->temp), &(scratch->weighted_h), &(scratch->k[0][i]));
        arpra_add(&(scratch->x_new_5[i]), &(system->x[i]), &(scratch->temp));
    }
    arpra_set_precision(&(scratch->temp), prec_t);
    arpra_add(&(scratch->temp), system->t, &(scratch->weighted_h));
    system->f(scratch->k[1],
              &(scratch->temp), scratch->x_new_5,
              system->dims, system->params);

    // k[2] = f(t + 3h/10, x(t) + 3h/40 k[0] + 9h/40 k[1])
    for (i = 0; i < system->dims; i++) {
        prec_x = arpra_get_precision(&(system->x[i]));
        arpra_set_precision(&(scratch->temp), prec_x);
        arpra_mul(&(scratch->weighted_h), &(scratch->a[2][0]), h);
        arpra_mul(&(scratch->temp), &(scratch->weighted_h), &(scratch->k[0][i]));
        arpra_add(&(scratch->x_new_5[i]), &(system->x[i]), &(scratch->temp));
        arpra_mul(&(scratch->weighted_h), &(scratch->a[2][1]), h);
        arpra_mul(&(scratch->temp), &(scratch->weighted_h), &(scratch->k[1][i]));
        arpra_add(&(scratch->x_new_5[i]), &(scratch->x_new_5[i]), &(scratch->temp));
    }
    arpra_set_precision(&(scratch->temp), prec_t);
    arpra_mul(&(scratch->weighted_h), &(scratch->c[2]), h);
    arpra_add(&(scratch->temp), system->t, &(scratch->weighted_h));
    system->f(scratch->k[2],
              &(scratch->temp), scratch->x_new_5,
              system->dims, system->params);

    // k[3] = f(t + 4h/5, x(t) + 44h/45 k[0] - 56h/15 k[1] + 32h/9 k[2])
    for (i = 0; i < system->dims; i++) {
        prec_x = arpra_get_precision(&(system->x[i]));
        arpra_set_precision(&(scratch->temp), prec_x);
        arpra_mul(&(scratch->weighted_h), &(scratch->a[3][0]), h);
        arpra_mul(&(scratch->temp), &(scratch->weighted_h), &(scratch->k[0][i]));
        arpra_add(&(scratch->x_new_5[i]), &(system->x[i]), &(scratch->temp));
        arpra_mul(&(scratch->weighted_h), &(scratch->a[3][1]), h);
        arpra_mul(&(scratch->temp), &(scratch->weighted_h), &(scratch->k[1][i]));
        arpra_add(&(scratch->x_new_5[i]), &(scratch->x_new_5[i]), &(scratch->temp));
        arpra_mul(&(scratch->weighted_h), &(scratch->a[3][2]), h);
        arpra_mul(&(scratch->temp), &(scratch->weighted_h), &(scratch->k[2][i]));
        arpra_add(&(scratch->x_new_5[i]), &(scratch->x_new_5[i]), &(scratch->temp));
    }
    arpra_set_precision(&(scratch->temp), prec_t);
    arpra_mul(&(scratch->weighted_h), &(scratch->c[3]), h);
    arpra_add(&(scratch->temp), system->t, &(scratch->weighted_h));
    system->f(scratch->k[3],
              &(scratch->temp), scratch->x_new_5,
              system->dims, system->params);

    // k[4] = f(t + 8h/9,
    //          x(t) + 19372h/6561 k[0]
    //               - 25360h/2187 k[1]
    //               + 64448h/6561 k[2]
    //               - 212h/729 k[3])


    // k[5] = f(t + h,
    //          x(t) + 9017h/3168 k[0]
    //               - 355h/33 k[1]
    //               + 46732h/5247 k[2]
    //               + 49h/176 k[3]
    //               - 5103h/18656 k[4])


    // k[6] = f(t + h,
    //          x(t) + 35h/384 k[0]
    //               + 500h/1113 k[2]
    //               + 125h/192 k[3]
    //               − 2187h/6784 k[4]
    //               + 11h/84 k[5])


    // Already been computed in x_new_5.
    // x_5(t + h) = x(t) + 35h/384 k[0]
    //                   + 500h/1113 k[2]
    //                   + 125h/192 k[3]
    //                   − 2187h/6784 k[4]
    //                   + 11h/84 k[5])


    // x_4(t + h) = x(t) + 5179h/57600 k[0]
    //                   + 7571h/16695 k[2]
    //                   + 393h/640 k[3]
    //                   - 92097h/339200 k[4]
    //                   + 187h/2100 k[5]
    //                   + 1h/40 k[6])


    for (i = 0; i < system->dims; i++) {
        prec_x = arpra_get_precision(&(system->x[i]));
        arpra_set_precision(&(scratch->temp), prec_x);

        arpra_mul(&(scratch->temp), &(scratch->weighted_h), &(scratch->k[0][i]));
        arpra_add(&(scratch->x_new_4[i]), &(system->x[i]), &(scratch->temp));

        arpra_mul(&(scratch->temp), &(scratch->weighted_h), &(scratch->k[1][i]));
        arpra_add(&(scratch->x_new_4[i]), &(scratch->x_new_4[i]), &(scratch->temp));

        arpra_mul(&(scratch->temp), &(scratch->weighted_h), &(scratch->k[2][i]));
        arpra_add(&(scratch->x_new_4[i]), &(scratch->x_new_4[i]), &(scratch->temp));

        arpra_mul(&(scratch->temp), &(scratch->weighted_h), &(scratch->k[3][i]));
        arpra_add(&(scratch->x_new_4[i]), &(scratch->x_new_4[i]), &(scratch->temp));
    }

    // Advance system.
    arpra_add(system->t, system->t, h);
    arpra_range *x_temp = system->x;
    system->x = scratch->x_new_5;
    scratch->x_new_5 = x_temp;
}

static const arpra_ode_method dopri54 =
{
    .init = &dopri54_init,
    .clear = &dopri54_clear,
    .step = &dopri54_step,
    .stages = dopri54_stages,
};

const arpra_ode_method *arpra_ode_dopri54 = &dopri54;
