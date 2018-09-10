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
    arpra_range a_[(bogsham32_stages * (bogsham32_stages - 1)) / 2];
    arpra_range *a[bogsham32_stages];
    arpra_range b_3[bogsham32_stages];
    arpra_range b_2[bogsham32_stages];
    arpra_range c[bogsham32_stages];
    arpra_range ah_[(bogsham32_stages * (bogsham32_stages - 1)) / 2];
    arpra_range *ah[bogsham32_stages];
    arpra_range bh_3[bogsham32_stages];
    arpra_range bh_2[bogsham32_stages];
    arpra_range ch[bogsham32_stages];
    arpra_range temp_t[bogsham32_stages];
    arpra_range temp_x;
} bogsham32_scratch;

static void bogsham32_compute_constants (arpra_ode_stepper *stepper, const arpra_precision prec)
{
    arpra_uint k_i, k_j;
    arpra_range numerator, denominator;
    bogsham32_scratch *scratch;

    scratch = (bogsham32_scratch *) stepper->scratch;

    // Init temp vars.
    arpra_init2(&numerator, prec);
    arpra_init2(&denominator, prec);

    // Update constant memory to internal precision.
    for (k_i = 0; k_i < bogsham32_stages; k_i++) {
        for (k_j = 0; k_j < k_i; k_j++) {
            arpra_set_precision(&(scratch->a[k_i][k_j]), prec);
        }
        arpra_set_precision(&(scratch->b_3[k_i]), prec);
        arpra_set_precision(&(scratch->b_2[k_i]), prec);
        arpra_set_precision(&(scratch->c[k_i]), prec);
    }

    // k[0] = f(t, x(t))
    arpra_set_zero(&(scratch->c[0]));

    // k[1] = f(t + 1/2 h,
    //          x(t) + 1/2 h k[0])
    arpra_set_d(&numerator, 1.);
    arpra_set_d(&denominator, 2.);
    arpra_div(&(scratch->c[1]), &numerator, &denominator);
    arpra_set(&(scratch->a[1][0]), &(scratch->c[1]));

    // k[2] = f(t + 3/4 h,
    //          x(t) + 0   h k[0]
    //               + 3/4 h k[1])
    arpra_set_d(&numerator, 3.);
    arpra_set_d(&denominator, 4.);
    arpra_div(&(scratch->c[2]), &numerator, &denominator);
    arpra_set_zero(&(scratch->a[2][0]));
    arpra_set(&(scratch->a[2][1]), &(scratch->c[2]));

    // k[3] = f(t + h,
    //          x(t) + 2/9 h k[0]
    //               + 1/3 h k[1]
    //               + 4/9 h k[2])
    arpra_set_d(&(scratch->c[3]), 1.);
    arpra_set_d(&numerator, 2.);
    arpra_set_d(&denominator, 9.);
    arpra_div(&(scratch->a[3][0]), &numerator, &denominator);
    arpra_set_d(&numerator, 1.);
    arpra_set_d(&denominator, 3.);
    arpra_div(&(scratch->a[3][1]), &numerator, &denominator);
    arpra_set_d(&numerator, 4.);
    arpra_set_d(&denominator, 9.);
    arpra_div(&(scratch->a[3][2]), &numerator, &denominator);

    // Already been computed in x_new_3.
    // x_3(t + h) = x(t) + 2/9 h k[0]
    //                   + 1/3 h k[1]
    //                   + 4/9 h k[2]
    //                   + 0   h k[3]
    arpra_set(&(scratch->b_3[0]), &(scratch->a[3][0]));
    arpra_set(&(scratch->b_3[1]), &(scratch->a[3][1]));
    arpra_set(&(scratch->b_3[2]), &(scratch->a[3][2]));
    arpra_set_zero(&(scratch->b_3[3]));

    // x_2(t + h) = x(t) + 7/24 h k[0]
    //                   + 1/4  h k[1]
    //                   + 1/3  h k[2]
    //                   + 1/8  h k[3]
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

    // Clear temp vars.
    arpra_clear(&numerator);
    arpra_clear(&denominator);
}

static void bogsham32_init (arpra_ode_stepper *stepper, arpra_ode_system *system)
{
    arpra_uint x_idx, k_i, k_j;
    arpra_precision prec_x, prec_internal;
    bogsham32_scratch *scratch;

    // Allocate scratch memory.
    scratch = malloc(sizeof(bogsham32_scratch));
    for (k_i = 0; k_i < bogsham32_stages; k_i++) {
        scratch->k[k_i] = malloc(system->dims * sizeof(arpra_range));
    }
    scratch->x_new_3 = malloc(system->dims * sizeof(arpra_range));
    scratch->x_new_2 = malloc(system->dims * sizeof(arpra_range));

    // Initialise scratch memory.
    prec_internal = arpra_get_internal_precision();
    for (x_idx = 0; x_idx < system->dims; x_idx++) {
        prec_x = arpra_get_precision(&(system->x[x_idx]));
        for (k_i = 0; k_i < bogsham32_stages; k_i++) {
            arpra_init2(&(scratch->k[k_i][x_idx]), prec_x);
        }
        arpra_init2(&(scratch->x_new_3[x_idx]), prec_x);
        arpra_init2(&(scratch->x_new_2[x_idx]), prec_x);
    }
    for (k_i = 0; k_i < bogsham32_stages; k_i++) {
        scratch->a[k_i] = &(scratch->a_[(k_i * (k_i - 1)) / 2]);
        scratch->ah[k_i] = &(scratch->ah_[(k_i * (k_i - 1)) / 2]);
        for (k_j = 0; k_j < k_i; k_j++) {
            arpra_init2(&(scratch->a[k_i][k_j]), prec_internal);
            arpra_init2(&(scratch->ah[k_i][k_j]), prec_internal);
        }
        arpra_init2(&(scratch->b_3[k_i]), prec_internal);
        arpra_init2(&(scratch->bh_3[k_i]), prec_internal);
        arpra_init2(&(scratch->b_2[k_i]), prec_internal);
        arpra_init2(&(scratch->bh_2[k_i]), prec_internal);
        arpra_init2(&(scratch->c[k_i]), prec_internal);
        arpra_init2(&(scratch->ch[k_i]), prec_internal);
        arpra_init2(&(scratch->temp_t[k_i]), prec_internal);
    }
    arpra_init2(&(scratch->temp_x), prec_internal);

    // Set stepper parameters.
    stepper->method = arpra_ode_bogsham32;
    stepper->system = system;
    stepper->error = NULL;
    stepper->scratch = scratch;

    // Precompute constants.
    bogsham32_compute_constants(stepper, prec_internal);
}

static void bogsham32_clear (arpra_ode_stepper *stepper)
{
    arpra_uint x_idx, k_i, k_j;
    arpra_ode_system *system;
    bogsham32_scratch *scratch;

    system = stepper->system;
    scratch = (bogsham32_scratch *) stepper->scratch;

    // Clear scratch memory.
    for (x_idx = 0; x_idx < system->dims; x_idx++) {
        for (k_i = 0; k_i < bogsham32_stages; k_i++) {
            arpra_clear(&(scratch->k[k_i][x_idx]));
        }
        arpra_clear(&(scratch->x_new_3[x_idx]));
        arpra_clear(&(scratch->x_new_2[x_idx]));
    }
    for (k_i = 0; k_i < bogsham32_stages; k_i++) {
        for (k_j = 0; k_j < k_i; k_j++) {
            arpra_clear(&(scratch->a[k_i][k_j]));
            arpra_clear(&(scratch->ah[k_i][k_j]));
        }
        arpra_clear(&(scratch->b_3[k_i]));
        arpra_clear(&(scratch->bh_3[k_i]));
        arpra_clear(&(scratch->b_2[k_i]));
        arpra_clear(&(scratch->bh_2[k_i]));
        arpra_clear(&(scratch->c[k_i]));
        arpra_clear(&(scratch->ch[k_i]));
        arpra_clear(&(scratch->temp_t[k_i]));
    }
    arpra_clear(&(scratch->temp_x));

    // Free scratch memory.
    for (k_i = 0; k_i < bogsham32_stages; k_i++) {
        free(scratch->k[k_i]);
    }
    free(scratch->x_new_3);
    free(scratch->x_new_2);
    free(scratch);
}

static void bogsham32_step (arpra_ode_stepper *stepper, const arpra_range *h)
{
    arpra_uint x_idx, k_i, k_j;
    arpra_precision prec_t, prec_x;
    arpra_range *x_old, *x_sum;
    arpra_ode_system *system;
    bogsham32_scratch *scratch;

    system = stepper->system;
    scratch = (bogsham32_scratch *) stepper->scratch;

    // Synchronise scratch precision and prepare step parameters.
    prec_t = arpra_get_precision(system->t);
    for (x_idx = 0; x_idx < system->dims; x_idx++) {
        prec_x = arpra_get_precision(&(system->x[x_idx]));
        for (k_i = 0; k_i < bogsham32_stages; k_i++) {
            arpra_set_precision(&(scratch->k[k_i][x_idx]), prec_x);
        }
        arpra_set_precision(&(scratch->x_new_3[x_idx]), prec_x);
        arpra_set_precision(&(scratch->x_new_2[x_idx]), prec_x);
    }
    for (k_i = 0; k_i < bogsham32_stages; k_i++) {
        for (k_j = 0; k_j < k_i; k_j++) {
            arpra_set_precision(&(scratch->ah[k_i][k_j]), prec_t);
            arpra_mul(&(scratch->ah[k_i][k_j]), &(scratch->a[k_i][k_j]), h);
        }
        arpra_set_precision(&(scratch->bh_3[k_i]), prec_t);
        arpra_mul(&(scratch->bh_3[k_i]), &(scratch->b_3[k_i]), h);
        arpra_set_precision(&(scratch->bh_2[k_i]), prec_t);
        arpra_mul(&(scratch->bh_2[k_i]), &(scratch->b_2[k_i]), h);
        arpra_set_precision(&(scratch->ch[k_i]), prec_t);
        arpra_mul(&(scratch->ch[k_i]), &(scratch->c[k_i]), h);
        arpra_set_precision(&(scratch->temp_t[k_i]), prec_t);
        arpra_add(&(scratch->temp_t[k_i]), system->t, &(scratch->ch[k_i]));
    }

    // Compute k stages and third-order approximation.
    for (k_i = 0; k_i < bogsham32_stages; k_i++) {
        x_old = (k_i == 0) ? system->x : scratch->x_new_3;

        // x(t + c_i h) = x(t) + a_i0 h k[0] + ... + a_is h k[s]
        for (x_idx = 0; x_idx < system->dims; x_idx++) {
            prec_x = arpra_get_precision(&(system->x[x_idx]));
            arpra_set_precision(&(scratch->temp_x), prec_x);
            for (k_j = 0; k_j < k_i; k_j++) {
                x_sum = (k_j == 0) ? system->x : scratch->x_new_3;
                arpra_mul(&(scratch->temp_x), &(scratch->ah[k_i][k_j]), &(scratch->k[k_j][x_idx]));
                arpra_add(&(scratch->x_new_3[x_idx]), &(x_sum[x_idx]), &(scratch->temp_x));
            }
        }

        // k[i] = f(t + c_i h, x(t) + a_i0 h k[0] + ... + a_is h k[s])
        for (x_idx = 0; x_idx < system->dims; x_idx++) {
            system->f(scratch->k[k_i],
                      &(scratch->temp_t[k_i]), x_old,
                      x_idx, system->params);
        }
    }

    // Compute second-order approximation.
    for (x_idx = 0; x_idx < system->dims; x_idx++) {
        prec_x = arpra_get_precision(&(system->x[x_idx]));
        arpra_set_precision(&(scratch->temp_x), prec_x);
        for (k_j = 0; k_j < bogsham32_stages; k_j++) {
            x_sum = (k_j == 0) ? system->x : scratch->x_new_2;
            arpra_mul(&(scratch->temp_x), &(scratch->bh_2[k_j]), &(scratch->k[k_j][x_idx]));
            arpra_add(&(scratch->x_new_2[x_idx]), &(x_sum[x_idx]), &(scratch->temp_x));
        }
    }

    // Advance system.
    arpra_add(system->t, system->t, h);
    for (x_idx = 0; x_idx < system->dims; x_idx++) {
        arpra_set(&(system->x[x_idx]), &(scratch->x_new_3[x_idx]));
    }
}

static const arpra_ode_method bogsham32 =
{
    .init = &bogsham32_init,
    .clear = &bogsham32_clear,
    .step = &bogsham32_step,
    .stages = bogsham32_stages,
};

const arpra_ode_method *arpra_ode_bogsham32 = &bogsham32;
