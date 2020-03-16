/*
 * ode_dopri54.c -- Dormand-Prince 5(4) ODE stepper.
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

#define dopri54_stages 7

typedef struct dopri54_scratch_struct
{
    arpra_range *_k[dopri54_stages];
    arpra_range **k[dopri54_stages];
    arpra_range *_x_new_5;
    arpra_range **x_new_5;
    arpra_range *_x_new_4;
    arpra_range **x_new_4;
    arpra_range _a[(dopri54_stages * (dopri54_stages - 1)) / 2];
    arpra_range *a[dopri54_stages];
    arpra_range b_5[dopri54_stages];
    arpra_range b_4[dopri54_stages];
    arpra_range c[dopri54_stages];
    arpra_range _ah[(dopri54_stages * (dopri54_stages - 1)) / 2];
    arpra_range *ah[dopri54_stages];
    arpra_range bh_5[dopri54_stages];
    arpra_range bh_4[dopri54_stages];
    arpra_range ch[dopri54_stages];
    arpra_range temp_t[dopri54_stages];
    arpra_range temp_x;
} dopri54_scratch;

static void dopri54_compute_constants (arpra_ode_stepper *stepper, const arpra_prec prec)
{
    arpra_uint k_i, k_j;
    arpra_range numerator, denominator;
    dopri54_scratch *scratch;

    scratch = (dopri54_scratch *) stepper->scratch;

    // Init temp vars.
    arpra_init2(&numerator, prec);
    arpra_init2(&denominator, prec);

    // Update constant memory to internal precision.
    for (k_i = 0; k_i < dopri54_stages; k_i++) {
        for (k_j = 0; k_j < k_i; k_j++) {
            arpra_set_precision(&(scratch->a[k_i][k_j]), prec);
        }
        arpra_set_precision(&(scratch->b_5[k_i]), prec);
        arpra_set_precision(&(scratch->b_4[k_i]), prec);
        arpra_set_precision(&(scratch->c[k_i]), prec);
    }

    // k[0] = f(t, x(t))
    arpra_set_zero(&(scratch->c[0]));

    // k[1] = f(t + 1/5 h,
    //          x(t) + 1/5 h k[0])
    arpra_set_d(&numerator, 1.);
    arpra_set_d(&denominator, 5.);
    arpra_div(&(scratch->c[1]), &numerator, &denominator);
    arpra_set(&(scratch->a[1][0]), &(scratch->c[1]));

    // k[2] = f(t + 3/10 h,
    //          x(t) + 3/40 h k[0]
    //               + 9/40 h k[1])
    arpra_set_d(&numerator, 3.);
    arpra_set_d(&denominator, 10.);
    arpra_div(&(scratch->c[2]), &numerator, &denominator);
    arpra_set_d(&numerator, 3.);
    arpra_set_d(&denominator, 40.);
    arpra_div(&(scratch->a[2][0]), &numerator, &denominator);
    arpra_set_d(&numerator, 9.);
    arpra_set_d(&denominator, 40.);
    arpra_div(&(scratch->a[2][1]), &numerator, &denominator);

    // k[3] = f(t + 4/5 h,
    //          x(t) + 44/45 h k[0]
    //               - 56/15 h k[1]
    //               + 32/9  h k[2])
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

    // k[4] = f(t + 8/9 h,
    //          x(t) + 19372/6561 h k[0]
    //               - 25360/2187 h k[1]
    //               + 64448/6561 h k[2]
    //               - 212/729    h k[3])
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

    // k[5] = f(t + h,
    //          x(t) + 9017/3168  h k[0]
    //               - 355/33     h k[1]
    //               + 46732/5247 h k[2]
    //               + 49/176     h k[3]
    //               - 5103/18656 h k[4])
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

    // k[6] = f(t + h,
    //          x(t) + 35/384    h k[0]
    //               + 0         h k[1]
    //               + 500/1113  h k[2]
    //               + 125/192   h k[3]
    //               - 2187/6784 h k[4]
    //               + 11/84     h k[5])
    arpra_set_d(&(scratch->c[6]), 1.);
    arpra_set_d(&numerator, 35.);
    arpra_set_d(&denominator, 384.);
    arpra_div(&(scratch->a[6][0]), &numerator, &denominator);
    arpra_set_zero(&(scratch->a[6][1]));
    arpra_set_d(&numerator, 500.);
    arpra_set_d(&denominator, 1113.);
    arpra_div(&(scratch->a[6][2]), &numerator, &denominator);
    arpra_set_d(&numerator, 125.);
    arpra_set_d(&denominator, 192.);
    arpra_div(&(scratch->a[6][3]), &numerator, &denominator);
    arpra_set_d(&numerator, -2187.);
    arpra_set_d(&denominator, 6784.);
    arpra_div(&(scratch->a[6][4]), &numerator, &denominator);
    arpra_set_d(&numerator, 11.);
    arpra_set_d(&denominator, 84.);
    arpra_div(&(scratch->a[6][5]), &numerator, &denominator);

    // Already been computed in x_new_5.
    // x_5(t + h) = x(t) + 35/384    h k[0]
    //                   + 0         h k[1]
    //                   + 500/1113  h k[2]
    //                   + 125/192   h k[3]
    //                   - 2187/6784 h k[4]
    //                   + 11/84     h k[5]
    //                   + 0         h k[6]
    arpra_set(&(scratch->b_5[0]), &(scratch->a[6][0]));
    arpra_set(&(scratch->b_5[1]), &(scratch->a[6][1]));
    arpra_set(&(scratch->b_5[2]), &(scratch->a[6][2]));
    arpra_set(&(scratch->b_5[3]), &(scratch->a[6][3]));
    arpra_set(&(scratch->b_5[4]), &(scratch->a[6][4]));
    arpra_set(&(scratch->b_5[5]), &(scratch->a[6][5]));
    arpra_set_zero(&(scratch->b_5[6]));

    // x_4(t + h) = x(t) + 5179/57600   h k[0]
    //                   + 0            h k[1]
    //                   + 7571/16695   h k[2]
    //                   + 393/640      h k[3]
    //                   - 92097/339200 h k[4]
    //                   + 187/2100     h k[5]
    //                   + 1/40         h k[6]
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

    // Clear temp vars.
    arpra_clear(&numerator);
    arpra_clear(&denominator);
}

static void dopri54_init (arpra_ode_stepper *stepper, arpra_ode_system *system)
{
    arpra_uint x_grp, x_dim, k_i, k_j, state_size;
    arpra_prec prec_x, prec_internal;
    dopri54_scratch *scratch;

    // Allocate scratch memory.
    scratch = malloc(sizeof(dopri54_scratch));
    for (x_grp = 0, state_size = 0; x_grp < system->grps; x_grp++) {
        state_size += system->dims[x_grp];
    }
    for (k_i = 0; k_i < dopri54_stages; k_i++) {
        scratch->_k[k_i] = malloc(state_size * sizeof(arpra_range));
        scratch->k[k_i] = malloc(system->grps * sizeof(arpra_range *));
    }
    scratch->_x_new_5 = malloc(state_size * sizeof(arpra_range));
    scratch->x_new_5 = malloc(system->grps * sizeof(arpra_range *));
    scratch->_x_new_4 = malloc(state_size * sizeof(arpra_range));
    scratch->x_new_4 = malloc(system->grps * sizeof(arpra_range *));

    // Initialise scratch memory.
    prec_internal = arpra_get_internal_precision();
    for (k_i = 0; k_i < dopri54_stages; k_i++) {
        scratch->k[k_i][0] = scratch->_k[k_i];
    }
    scratch->x_new_5[0] = scratch->_x_new_5;
    scratch->x_new_4[0] = scratch->_x_new_4;
    for (x_grp = 1; x_grp < system->grps; x_grp++) {
        for (k_i = 0; k_i < dopri54_stages; k_i++) {
            scratch->k[k_i][x_grp] = scratch->k[k_i][x_grp - 1] + system->dims[x_grp - 1];
        }
        scratch->x_new_5[x_grp] = scratch->x_new_5[x_grp - 1] + system->dims[x_grp - 1];
        scratch->x_new_4[x_grp] = scratch->x_new_4[x_grp - 1] + system->dims[x_grp - 1];
    }
    for (x_grp = 0; x_grp < system->grps; x_grp++) {
        for (x_dim = 0; x_dim < system->dims[x_grp]; x_dim++) {
            prec_x = arpra_get_precision(&(system->x[x_grp][x_dim]));
            for (k_i = 0; k_i < dopri54_stages; k_i++) {
                arpra_init2(&(scratch->k[k_i][x_grp][x_dim]), prec_x);
            }
            arpra_init2(&(scratch->x_new_5[x_grp][x_dim]), prec_x);
            arpra_init2(&(scratch->x_new_4[x_grp][x_dim]), prec_x);
        }
    }
    for (k_i = 0; k_i < dopri54_stages; k_i++) {
        scratch->a[k_i] = &(scratch->_a[(k_i * (k_i - 1)) / 2]);
        scratch->ah[k_i] = &(scratch->_ah[(k_i * (k_i - 1)) / 2]);
        for (k_j = 0; k_j < k_i; k_j++) {
            arpra_init2(&(scratch->a[k_i][k_j]), prec_internal);
            arpra_init2(&(scratch->ah[k_i][k_j]), prec_internal);
        }
        arpra_init2(&(scratch->b_5[k_i]), prec_internal);
        arpra_init2(&(scratch->bh_5[k_i]), prec_internal);
        arpra_init2(&(scratch->b_4[k_i]), prec_internal);
        arpra_init2(&(scratch->bh_4[k_i]), prec_internal);
        arpra_init2(&(scratch->c[k_i]), prec_internal);
        arpra_init2(&(scratch->ch[k_i]), prec_internal);
        arpra_init2(&(scratch->temp_t[k_i]), prec_internal);
    }
    arpra_init2(&(scratch->temp_x), prec_internal);

    // Set stepper parameters.
    stepper->method = arpra_ode_dopri54;
    stepper->system = system;
    stepper->error = NULL;
    stepper->scratch = scratch;

    // Precompute constants.
    dopri54_compute_constants(stepper, prec_internal);
}

static void dopri54_clear (arpra_ode_stepper *stepper)
{
    arpra_uint x_grp, x_dim, k_i, k_j;
    arpra_ode_system *system;
    dopri54_scratch *scratch;

    system = stepper->system;
    scratch = (dopri54_scratch *) stepper->scratch;

    // Clear scratch memory.
    for (x_grp = 0; x_grp < system->grps; x_grp++) {
        for (x_dim = 0; x_dim < system->dims[x_grp]; x_dim++) {
            for (k_i = 0; k_i < dopri54_stages; k_i++) {
                arpra_clear(&(scratch->k[k_i][x_grp][x_dim]));
            }
            arpra_clear(&(scratch->x_new_5[x_grp][x_dim]));
            arpra_clear(&(scratch->x_new_4[x_grp][x_dim]));
        }
    }
    for (k_i = 0; k_i < dopri54_stages; k_i++) {
        for (k_j = 0; k_j < k_i; k_j++) {
            arpra_clear(&(scratch->a[k_i][k_j]));
            arpra_clear(&(scratch->ah[k_i][k_j]));
        }
        arpra_clear(&(scratch->b_5[k_i]));
        arpra_clear(&(scratch->bh_5[k_i]));
        arpra_clear(&(scratch->b_4[k_i]));
        arpra_clear(&(scratch->bh_4[k_i]));
        arpra_clear(&(scratch->c[k_i]));
        arpra_clear(&(scratch->ch[k_i]));
        arpra_clear(&(scratch->temp_t[k_i]));
    }
    arpra_clear(&(scratch->temp_x));

    // Free scratch memory.
    for (k_i = 0; k_i < dopri54_stages; k_i++) {
        free(scratch->_k[k_i]);
        free(scratch->k[k_i]);
    }
    free(scratch->_x_new_5);
    free(scratch->x_new_5);
    free(scratch->_x_new_4);
    free(scratch->x_new_4);
    free(scratch);
}

static void dopri54_step (arpra_ode_stepper *stepper, const arpra_range *h)
{
    arpra_uint x_grp, x_dim, k_i, k_j;
    arpra_prec prec_t, prec_x;
    arpra_range **x_old, **x_sum;
    arpra_ode_system *system;
    dopri54_scratch *scratch;

    system = stepper->system;
    scratch = (dopri54_scratch *) stepper->scratch;

    // Synchronise scratch precision and prepare step parameters.
    prec_t = arpra_get_precision(system->t);
    for (x_grp = 0; x_grp < system->grps; x_grp++) {
        for (x_dim = 0; x_dim < system->dims[x_grp]; x_dim++) {
            prec_x = arpra_get_precision(&(system->x[x_grp][x_dim]));
            for (k_i = 0; k_i < dopri54_stages; k_i++) {
                arpra_set_precision(&(scratch->k[k_i][x_grp][x_dim]), prec_x);
            }
            arpra_set_precision(&(scratch->x_new_5[x_grp][x_dim]), prec_x);
            arpra_set_precision(&(scratch->x_new_4[x_grp][x_dim]), prec_x);
        }
    }
    for (k_i = 0; k_i < dopri54_stages; k_i++) {
        for (k_j = 0; k_j < k_i; k_j++) {
            arpra_set_precision(&(scratch->ah[k_i][k_j]), prec_t);
            arpra_mul(&(scratch->ah[k_i][k_j]), &(scratch->a[k_i][k_j]), h);
        }
        arpra_set_precision(&(scratch->bh_5[k_i]), prec_t);
        arpra_mul(&(scratch->bh_5[k_i]), &(scratch->b_5[k_i]), h);
        arpra_set_precision(&(scratch->bh_4[k_i]), prec_t);
        arpra_mul(&(scratch->bh_4[k_i]), &(scratch->b_4[k_i]), h);
        arpra_set_precision(&(scratch->ch[k_i]), prec_t);
        arpra_mul(&(scratch->ch[k_i]), &(scratch->c[k_i]), h);
        arpra_set_precision(&(scratch->temp_t[k_i]), prec_t);
        arpra_add(&(scratch->temp_t[k_i]), system->t, &(scratch->ch[k_i]));
    }

    // Compute k stages and fifth-order approximation.
    for (k_i = 0; k_i < dopri54_stages; k_i++) {
        x_old = (k_i == 0) ? system->x : scratch->x_new_5;

        // x(t + c_i h) = x(t) + a_i0 h k[0] + ... + a_is h k[s]
        for (x_grp = 0; x_grp < system->grps; x_grp++) {
            for (x_dim = 0; x_dim < system->dims[x_grp]; x_dim++) {
                prec_x = arpra_get_precision(&(system->x[x_grp][x_dim]));
                arpra_set_precision(&(scratch->temp_x), prec_x);
                for (k_j = 0; k_j < k_i; k_j++) {
                    x_sum = (k_j == 0) ? system->x : scratch->x_new_5;
                    arpra_mul(&(scratch->temp_x), &(scratch->ah[k_i][k_j]), &(scratch->k[k_j][x_grp][x_dim]));
                    arpra_add(&(scratch->x_new_5[x_grp][x_dim]), &(x_sum[x_grp][x_dim]), &(scratch->temp_x));
                }
            }
        }

        // k[i] = f(t + c_i h, x(t) + a_i0 h k[0] + ... + a_is h k[s])
        for (x_grp = 0; x_grp < system->grps; x_grp++) {
            for (x_dim = 0; x_dim < system->dims[x_grp]; x_dim++) {
                system->f[x_grp](&(scratch->k[k_i][x_grp][x_dim]), system->params[x_grp],
                                 &(scratch->temp_t[k_i]), (const arpra_range **) x_old, x_grp, x_dim);
            }
        }
    }

    // Compute fourth-order approximation.
    for (x_grp = 0; x_grp < system->grps; x_grp++) {
        for (x_dim = 0; x_dim < system->dims[x_grp]; x_dim++) {
            prec_x = arpra_get_precision(&(system->x[x_grp][x_dim]));
            arpra_set_precision(&(scratch->temp_x), prec_x);
            for (k_j = 0; k_j < dopri54_stages; k_j++) {
                x_sum = (k_j == 0) ? system->x : scratch->x_new_4;
                arpra_mul(&(scratch->temp_x), &(scratch->bh_4[k_j]), &(scratch->k[k_j][x_grp][x_dim]));
                arpra_add(&(scratch->x_new_4[x_grp][x_dim]), &(x_sum[x_grp][x_dim]), &(scratch->temp_x));
            }
        }
    }

    // Advance system.
    arpra_add(system->t, system->t, h);
    for (x_grp = 0; x_grp < system->grps; x_grp++) {
        for (x_dim = 0; x_dim < system->dims[x_grp]; x_dim++) {
            arpra_set(&(system->x[x_grp][x_dim]), &(scratch->x_new_5[x_grp][x_dim]));
        }
    }
}

static const arpra_ode_method dopri54 =
{
    .init = &dopri54_init,
    .clear = &dopri54_clear,
    .step = &dopri54_step,
    .stages = dopri54_stages,
};

const arpra_ode_method *arpra_ode_dopri54 = &dopri54;
