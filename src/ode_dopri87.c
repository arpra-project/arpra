/*
 * ode_dopri87.c -- Dormand-Prince 8(7) ODE stepper.
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

#define dopri87_stages 13

typedef struct dopri87_scratch_struct
{
    arpra_range *k[dopri87_stages];
    arpra_range *x_new_8;
    arpra_range *x_new_7;
    arpra_range a_[(dopri87_stages * (dopri87_stages - 1)) / 2];
    arpra_range *a[dopri87_stages];
    arpra_range b_8[dopri87_stages];
    arpra_range b_7[dopri87_stages];
    arpra_range c[dopri87_stages];
    arpra_range ah_[(dopri87_stages * (dopri87_stages - 1)) / 2];
    arpra_range *ah[dopri87_stages];
    arpra_range bh_8[dopri87_stages];
    arpra_range bh_7[dopri87_stages];
    arpra_range ch[dopri87_stages];
    arpra_range temp_t[dopri87_stages];
    arpra_range temp_x;
} dopri87_scratch;

static void dopri87_compute_constants (arpra_ode_stepper *stepper, const arpra_precision prec)
{
    arpra_uint k_i, k_j;
    arpra_range numerator, denominator;
    dopri87_scratch *scratch;

    scratch = (dopri87_scratch *) stepper->scratch;

    // Init temp vars.
    arpra_init2(&numerator, prec);
    arpra_init2(&denominator, prec);

    // Update constant memory to internal precision.
    for (k_i = 0; k_i < dopri87_stages; k_i++) {
        for (k_j = 0; k_j < k_i; k_j++) {
            arpra_set_precision(&(scratch->a[k_i][k_j]), prec);
        }
        arpra_set_precision(&(scratch->b_8[k_i]), prec);
        arpra_set_precision(&(scratch->b_7[k_i]), prec);
        arpra_set_precision(&(scratch->c[k_i]), prec);
    }

    // k[0] = f(t, x(t))
    arpra_set_zero(&(scratch->c[0]));

    // k[1] = f(t + 1/18 h,
    //          x(t) + 1/18 h k[0])
    arpra_set_d(&numerator, 1.);
    arpra_set_d(&denominator, 18.);
    arpra_div(&(scratch->c[1]), &numerator, &denominator);
    arpra_set(&(scratch->a[1][0]), &(scratch->c[1]));

    // k[2] = f(t + 1/12 h,
    //          x(t) + 1/48 h k[0]
    //               + 1/16 h k[1])
    arpra_set_d(&numerator, 1.);
    arpra_set_d(&denominator, 12.);
    arpra_div(&(scratch->c[2]), &numerator, &denominator);
    arpra_set_d(&denominator, 48.);
    arpra_div(&(scratch->a[2][0]), &numerator, &denominator);
    arpra_set_d(&denominator, 16.);
    arpra_div(&(scratch->a[2][1]), &numerator, &denominator);

    // k[3] = f(t + 1/8 h,
    //          x(t) + 1/32 h k[0]
    //               + 0    h k[1]
    //               + 3/32 h k[2])
    arpra_set_d(&numerator, 1.);
    arpra_set_d(&denominator, 8.);
    arpra_div(&(scratch->c[3]), &numerator, &denominator);
    arpra_set_d(&denominator, 32.);
    arpra_div(&(scratch->a[3][0]), &numerator, &denominator);
    arpra_set_zero(&(scratch->a[3][1]));
    arpra_set_d(&numerator, 3.);
    arpra_div(&(scratch->a[3][2]), &numerator, &denominator);

    // k[4] = f(t + 5/16 h,
    //          x(t) + 5/16  h k[0]
    //               + 0     h k[1]
    //               - 75/64 h k[2]
    //               + 75/64 h k[3])
    arpra_set_d(&numerator, 5.);
    arpra_set_d(&denominator, 16.);
    arpra_div(&(scratch->c[4]), &numerator, &denominator);
    arpra_set(&(scratch->a[4][0]), &(scratch->c[4]));
    arpra_set_zero(&(scratch->a[4][1]));
    arpra_set_d(&numerator, -75.);
    arpra_set_d(&denominator, 64.);
    arpra_div(&(scratch->a[4][2]), &numerator, &denominator);
    arpra_set_d(&numerator, 75.);
    arpra_div(&(scratch->a[4][3]), &numerator, &denominator);

    // k[5] = f(t + 3/8 h,
    //          x(t) + 3/80 h k[0]
    //               + 0    h k[1]
    //               + 0    h k[2]
    //               + 3/16 h k[3]
    //               + 3/20 h k[4])
    arpra_set_d(&numerator, 3.);
    arpra_set_d(&denominator, 8.);
    arpra_div(&(scratch->c[5]), &numerator, &denominator);
    arpra_set_d(&denominator, 80.);
    arpra_div(&(scratch->a[5][0]), &numerator, &denominator);
    arpra_set_zero(&(scratch->a[5][1]));
    arpra_set_zero(&(scratch->a[5][2]));
    arpra_set_d(&denominator, 16.);
    arpra_div(&(scratch->a[5][3]), &numerator, &denominator);
    arpra_set_d(&denominator, 20.);
    arpra_div(&(scratch->a[5][4]), &numerator, &denominator);

    // k[6] = f(t + 59/400 h,
    //          x(t) + 29443841/614563906  h k[0]
    //               + 0                   h k[1]
    //               + 0                   h k[2]
    //               + 77736538/692538347  h k[3]
    //               - 28693883/1125000000 h k[4]
    //               + 23124283/1800000000 h k[5])
    arpra_set_d(&numerator, 59.);
    arpra_set_d(&denominator, 400.);
    arpra_div(&(scratch->c[6]), &numerator, &denominator);
    arpra_set_d(&numerator, 29443841.);
    arpra_set_d(&denominator, 614563906.);
    arpra_div(&(scratch->a[6][0]), &numerator, &denominator);
    arpra_set_zero(&(scratch->a[6][1]));
    arpra_set_zero(&(scratch->a[6][2]));
    arpra_set_d(&numerator, 77736538.);
    arpra_set_d(&denominator, 692538347.);
    arpra_div(&(scratch->a[6][3]), &numerator, &denominator);
    arpra_set_d(&numerator, -28693883.);
    arpra_set_d(&denominator, 1125000000.);
    arpra_div(&(scratch->a[6][4]), &numerator, &denominator);
    arpra_set_d(&numerator, 23124283.);
    arpra_set_d(&denominator, 1800000000.);
    arpra_div(&(scratch->a[6][5]), &numerator, &denominator);

    // k[7] = f(t + 93/200 h,
    //          x(t) + 16016141/946692911   h k[0]
    //               + 0                    h k[1]
    //               + 0                    h k[2]
    //               + 61564180/158732637   h k[3]
    //               + 22789713/633445777   h k[4]
    //               + 545815736/2771057229 h k[5]
    //               - 180193667/1043307555 h k[6])
    arpra_set_d(&numerator, 93.);
    arpra_set_d(&denominator, 200.);
    arpra_div(&(scratch->c[7]), &numerator, &denominator);
    arpra_set_d(&numerator, 16016141.);
    arpra_set_d(&denominator, 946692911.);
    arpra_div(&(scratch->a[7][0]), &numerator, &denominator);
    arpra_set_zero(&(scratch->a[7][1]));
    arpra_set_zero(&(scratch->a[7][2]));
    arpra_set_d(&numerator, 61564180.);
    arpra_set_d(&denominator, 158732637.);
    arpra_div(&(scratch->a[7][3]), &numerator, &denominator);
    arpra_set_d(&numerator, 22789713.);
    arpra_set_d(&denominator, 633445777.);
    arpra_div(&(scratch->a[7][4]), &numerator, &denominator);
    arpra_set_d(&numerator, 545815736.);
    arpra_set_d(&denominator, 2771057229.);
    arpra_div(&(scratch->a[7][5]), &numerator, &denominator);
    arpra_set_d(&numerator, -180193667.);
    arpra_set_d(&denominator, 1043307555.);
    arpra_div(&(scratch->a[7][6]), &numerator, &denominator);

    // k[8] = f(t + 5490023248/9719169821 h,
    //          x(t) + 39632708/573591083   h k[0]
    //               + 0                    h k[1]
    //               + 0                    h k[2]
    //               - 433636366/683701615  h k[3]
    //               - 421739975/2616292301 h k[4]
    //               + 100302831/723423059  h k[5]
    //               + 790204164/839813087  h k[6]
    //               + 800635310/3783071287 h k[7])
    arpra_set_d(&numerator, 5490023248.);
    arpra_set_d(&denominator, 9719169821.);
    arpra_div(&(scratch->c[8]), &numerator, &denominator);
    arpra_set_d(&numerator, 39632708.);
    arpra_set_d(&denominator, 573591083.);
    arpra_div(&(scratch->a[8][0]), &numerator, &denominator);
    arpra_set_zero(&(scratch->a[8][1]));
    arpra_set_zero(&(scratch->a[8][2]));
    arpra_set_d(&numerator, -433636366.);
    arpra_set_d(&denominator, 683701615.);
    arpra_div(&(scratch->a[8][3]), &numerator, &denominator);
    arpra_set_d(&numerator, -421739975.);
    arpra_set_d(&denominator, 2616292301.);
    arpra_div(&(scratch->a[8][4]), &numerator, &denominator);
    arpra_set_d(&numerator, 100302831.);
    arpra_set_d(&denominator, 723423059.);
    arpra_div(&(scratch->a[8][5]), &numerator, &denominator);
    arpra_set_d(&numerator, 790204164.);
    arpra_set_d(&denominator, 839813087.);
    arpra_div(&(scratch->a[8][6]), &numerator, &denominator);
    arpra_set_d(&numerator, 800635310.);
    arpra_set_d(&denominator, 3783071287.);
    arpra_div(&(scratch->a[8][7]), &numerator, &denominator);

    // k[9] = f(t + 13/20 h,
    //          x(t) + 246121993/1340847787    h k[0]
    //               + 0                       h k[1]
    //               + 0                       h k[2]
    //               - 37695042795/15268766246 h k[3]
    //               - 309121744/1061227803    h k[4]
    //               - 12992083/490766935      h k[5]
    //               + 6005943493/2108947869   h k[6]
    //               + 393006217/1396673457    h k[7]
    //               + 123872331/1001029789    h k[8])
    arpra_set_d(&numerator, 13.);
    arpra_set_d(&denominator, 20.);
    arpra_div(&(scratch->c[9]), &numerator, &denominator);
    arpra_set_d(&numerator, 246121993.);
    arpra_set_d(&denominator, 1340847787.);
    arpra_div(&(scratch->a[9][0]), &numerator, &denominator);
    arpra_set_zero(&(scratch->a[9][1]));
    arpra_set_zero(&(scratch->a[9][2]));
    arpra_set_d(&numerator, -37695042795.);
    arpra_set_d(&denominator, 15268766246.);
    arpra_div(&(scratch->a[9][3]), &numerator, &denominator);
    arpra_set_d(&numerator, -309121744.);
    arpra_set_d(&denominator, 1061227803.);
    arpra_div(&(scratch->a[9][4]), &numerator, &denominator);
    arpra_set_d(&numerator, -12992083.);
    arpra_set_d(&denominator, 490766935.);
    arpra_div(&(scratch->a[9][5]), &numerator, &denominator);
    arpra_set_d(&numerator, 6005943493.);
    arpra_set_d(&denominator, 2108947869.);
    arpra_div(&(scratch->a[9][6]), &numerator, &denominator);
    arpra_set_d(&numerator, 393006217.);
    arpra_set_d(&denominator, 1396673457.);
    arpra_div(&(scratch->a[9][7]), &numerator, &denominator);
    arpra_set_d(&numerator, 123872331.);
    arpra_set_d(&denominator, 1001029789.);
    arpra_div(&(scratch->a[9][8]), &numerator, &denominator);

    // k[10] = f(t + 1201146811/1299019798 h,
    //           x(t) - 1028468189/846180014   h k[0]
    //                + 0                      h k[1]
    //                + 0                      h k[2]
    //                + 8478235783/508512852   h k[3]
    //                + 1311729495/1432422823  h k[4]
    //                - 10304129995/1701304382 h k[5]
    //                - 48777925059/3047939560 h k[6]
    //                + 15336726248/1032824649 h k[7]
    //                - 45442868181/3398467696 h k[8]
    //                + 3065993473/597172653   h k[9])
    arpra_set_d(&numerator, 1201146811.);
    arpra_set_d(&denominator, 1299019798.);
    arpra_div(&(scratch->c[10]), &numerator, &denominator);
    arpra_set_d(&numerator, -1028468189.);
    arpra_set_d(&denominator, 846180014.);
    arpra_div(&(scratch->a[10][0]), &numerator, &denominator);
    arpra_set_zero(&(scratch->a[10][1]));
    arpra_set_zero(&(scratch->a[10][2]));
    arpra_set_d(&numerator, 8478235783.);
    arpra_set_d(&denominator, 508512852.);
    arpra_div(&(scratch->a[10][3]), &numerator, &denominator);
    arpra_set_d(&numerator, 1311729495.);
    arpra_set_d(&denominator, 1432422823.);
    arpra_div(&(scratch->a[10][4]), &numerator, &denominator);
    arpra_set_d(&numerator, -10304129995.);
    arpra_set_d(&denominator, 1701304382.);
    arpra_div(&(scratch->a[10][5]), &numerator, &denominator);
    arpra_set_d(&numerator, -48777925059.);
    arpra_set_d(&denominator, 3047939560.);
    arpra_div(&(scratch->a[10][6]), &numerator, &denominator);
    arpra_set_d(&numerator, 15336726248.);
    arpra_set_d(&denominator, 1032824649.);
    arpra_div(&(scratch->a[10][7]), &numerator, &denominator);
    arpra_set_d(&numerator, -45442868181.);
    arpra_set_d(&denominator, 3398467696.);
    arpra_div(&(scratch->a[10][8]), &numerator, &denominator);
    arpra_set_d(&numerator, 3065993473.);
    arpra_set_d(&denominator, 597172653.);
    arpra_div(&(scratch->a[10][9]), &numerator, &denominator);

    // k[11] = f(t + h,
    //           x(t) + 185892177/718116043   h k[0]
    //                + 0                     h k[1]
    //                + 0                     h k[2]
    //                - 3185094517/667107341  h k[3]
    //                - 477755414/1098053517  h k[4]
    //                - 703635378/230739211   h k[5]
    //                + 5731566787/1027545527 h k[6]
    //                + 5232866602/850066563  h k[7]
    //                - 4093664535/808688257  h k[8]
    //                + 3962137247/1805957418 h k[9]
    //                + 65686358/487910083    h k[10])
    arpra_set_d(&(scratch->c[11]), 1.);
    arpra_set_d(&numerator, 185892177.);
    arpra_set_d(&denominator, 718116043.);
    arpra_div(&(scratch->a[11][0]), &numerator, &denominator);
    arpra_set_zero(&(scratch->a[11][1]));
    arpra_set_zero(&(scratch->a[11][2]));
    arpra_set_d(&numerator, -3185094517.);
    arpra_set_d(&denominator, 667107341.);
    arpra_div(&(scratch->a[11][3]), &numerator, &denominator);
    arpra_set_d(&numerator, -477755414.);
    arpra_set_d(&denominator, 1098053517.);
    arpra_div(&(scratch->a[11][4]), &numerator, &denominator);
    arpra_set_d(&numerator, -703635378.);
    arpra_set_d(&denominator, 230739211.);
    arpra_div(&(scratch->a[11][5]), &numerator, &denominator);
    arpra_set_d(&numerator, 5731566787.);
    arpra_set_d(&denominator, 1027545527.);
    arpra_div(&(scratch->a[11][6]), &numerator, &denominator);
    arpra_set_d(&numerator, 5232866602.);
    arpra_set_d(&denominator, 850066563.);
    arpra_div(&(scratch->a[11][7]), &numerator, &denominator);
    arpra_set_d(&numerator, -4093664535.);
    arpra_set_d(&denominator, 808688257.);
    arpra_div(&(scratch->a[11][8]), &numerator, &denominator);
    arpra_set_d(&numerator, 3962137247.);
    arpra_set_d(&denominator, 1805957418.);
    arpra_div(&(scratch->a[11][9]), &numerator, &denominator);
    arpra_set_d(&numerator, 65686358.);
    arpra_set_d(&denominator, 487910083.);
    arpra_div(&(scratch->a[11][10]), &numerator, &denominator);

    // k[12] = f(t + h,
    //           x(t) + 403863854/491063109    h k[0]
    //                + 0                      h k[1]
    //                + 0                      h k[2]
    //                - 5068492393/434740067   h k[3]
    //                - 411421997/543043805    h k[4]
    //                + 652783627/914296604    h k[5]
    //                + 11173962825/925320556  h k[6]
    //                - 13158990841/6184727034 h k[7]
    //                + 3936647629/1978049680  h k[8]
    //                - 160528059/685178525    h k[9]
    //                + 248638103/1413531060   h k[10]
    //                + 0                      h k[11])
    arpra_set_d(&(scratch->c[12]), 1.);
    arpra_set_d(&numerator, 403863854.);
    arpra_set_d(&denominator, 491063109.);
    arpra_div(&(scratch->a[12][0]), &numerator, &denominator);
    arpra_set_zero(&(scratch->a[12][1]));
    arpra_set_zero(&(scratch->a[12][2]));
    arpra_set_d(&numerator, -5068492393.);
    arpra_set_d(&denominator, 434740067.);
    arpra_div(&(scratch->a[12][3]), &numerator, &denominator);
    arpra_set_d(&numerator, -411421997.);
    arpra_set_d(&denominator, 543043805.);
    arpra_div(&(scratch->a[12][4]), &numerator, &denominator);
    arpra_set_d(&numerator, 652783627.);
    arpra_set_d(&denominator, 914296604.);
    arpra_div(&(scratch->a[12][5]), &numerator, &denominator);
    arpra_set_d(&numerator, 11173962825.);
    arpra_set_d(&denominator, 925320556.);
    arpra_div(&(scratch->a[12][6]), &numerator, &denominator);
    arpra_set_d(&numerator, -13158990841.);
    arpra_set_d(&denominator, 6184727034.);
    arpra_div(&(scratch->a[12][7]), &numerator, &denominator);
    arpra_set_d(&numerator, 3936647629.);
    arpra_set_d(&denominator, 1978049680.);
    arpra_div(&(scratch->a[12][8]), &numerator, &denominator);
    arpra_set_d(&numerator, -160528059.);
    arpra_set_d(&denominator, 685178525.);
    arpra_div(&(scratch->a[12][9]), &numerator, &denominator);
    arpra_set_d(&numerator, 248638103.);
    arpra_set_d(&denominator, 1413531060.);
    arpra_div(&(scratch->a[12][10]), &numerator, &denominator);
    arpra_set_zero(&(scratch->a[12][11]));

    // x_8(t + h) = x(t) + 14005451/335480064    h k[0]
    //                   + 0                     h k[1]
    //                   + 0                     h k[2]
    //                   + 0                     h k[3]
    //                   + 0                     h k[4]
    //                   - 59238493/1068277825   h k[5]
    //                   + 181606767/758867731   h k[6]
    //                   + 561292985/797845732   h k[7]
    //                   - 1041891430/1371343529 h k[8]
    //                   + 760417239/1151165299  h k[9]
    //                   + 118820643/751138087   h k[10]
    //                   - 528747749/2220607170  h k[11]
    //                   + 1/4                   h k[12]
    arpra_set_d(&numerator, 14005451.);
    arpra_set_d(&denominator, 335480064.);
    arpra_div(&(scratch->b_8[0]), &numerator, &denominator);
    arpra_set_zero(&(scratch->b_8[1]));
    arpra_set_zero(&(scratch->b_8[2]));
    arpra_set_zero(&(scratch->b_8[3]));
    arpra_set_zero(&(scratch->b_8[4]));
    arpra_set_d(&numerator, -59238493.);
    arpra_set_d(&denominator, 1068277825.);
    arpra_div(&(scratch->b_8[5]), &numerator, &denominator);
    arpra_set_d(&numerator, 181606767.);
    arpra_set_d(&denominator, 758867731.);
    arpra_div(&(scratch->b_8[6]), &numerator, &denominator);
    arpra_set_d(&numerator, 561292985.);
    arpra_set_d(&denominator, 797845732.);
    arpra_div(&(scratch->b_8[7]), &numerator, &denominator);
    arpra_set_d(&numerator, -1041891430.);
    arpra_set_d(&denominator, 1371343529.);
    arpra_div(&(scratch->b_8[8]), &numerator, &denominator);
    arpra_set_d(&numerator, 760417239.);
    arpra_set_d(&denominator, 1151165299.);
    arpra_div(&(scratch->b_8[9]), &numerator, &denominator);
    arpra_set_d(&numerator, 118820643.);
    arpra_set_d(&denominator, 751138087.);
    arpra_div(&(scratch->b_8[10]), &numerator, &denominator);
    arpra_set_d(&numerator, -528747749.);
    arpra_set_d(&denominator, 2220607170.);
    arpra_div(&(scratch->b_8[11]), &numerator, &denominator);
    arpra_set_d(&numerator, 1.);
    arpra_set_d(&denominator, 4.);
    arpra_div(&(scratch->b_8[12]), &numerator, &denominator);

    // x_7(t + h) = x(t) + 13451932/455176623    h k[0]
    //                   + 0                     h k[1]
    //                   + 0                     h k[2]
    //                   + 0                     h k[3]
    //                   + 0                     h k[4]
    //                   - 808719846/976000145   h k[5]
    //                   + 1757004468/5645159321 h k[6]
    //                   + 656045339/265891186   h k[7]
    //                   - 3867574721/1518517206 h k[8]
    //                   + 465885868/322736535   h k[9]
    //                   + 53011238/667516719    h k[10]
    //                   + 2/45                  h k[11]
    //                   + 0                     h k[12]
    arpra_set_d(&numerator, 13451932.);
    arpra_set_d(&denominator, 455176623.);
    arpra_div(&(scratch->b_7[0]), &numerator, &denominator);
    arpra_set_zero(&(scratch->b_7[1]));
    arpra_set_zero(&(scratch->b_7[2]));
    arpra_set_zero(&(scratch->b_7[3]));
    arpra_set_zero(&(scratch->b_7[4]));
    arpra_set_d(&numerator, -808719846.);
    arpra_set_d(&denominator, 976000145.);
    arpra_div(&(scratch->b_7[5]), &numerator, &denominator);
    arpra_set_d(&numerator, 1757004468.);
    arpra_set_d(&denominator, 5645159321.);
    arpra_div(&(scratch->b_7[6]), &numerator, &denominator);
    arpra_set_d(&numerator, 656045339.);
    arpra_set_d(&denominator, 265891186.);
    arpra_div(&(scratch->b_7[7]), &numerator, &denominator);
    arpra_set_d(&numerator, -3867574721.);
    arpra_set_d(&denominator, 1518517206.);
    arpra_div(&(scratch->b_7[8]), &numerator, &denominator);
    arpra_set_d(&numerator, 465885868.);
    arpra_set_d(&denominator, 322736535.);
    arpra_div(&(scratch->b_7[9]), &numerator, &denominator);
    arpra_set_d(&numerator, 53011238.);
    arpra_set_d(&denominator, 667516719.);
    arpra_div(&(scratch->b_7[10]), &numerator, &denominator);
    arpra_set_d(&numerator, 2.);
    arpra_set_d(&denominator, 45.);
    arpra_div(&(scratch->b_7[11]), &numerator, &denominator);
    arpra_set_zero(&(scratch->b_7[12]));

    // Clear temp vars.
    arpra_clear(&numerator);
    arpra_clear(&denominator);
}

static void dopri87_init (arpra_ode_stepper *stepper, arpra_ode_system *system)
{
    arpra_uint x_i, k_i, k_j;
    arpra_precision prec_x, prec_internal;
    dopri87_scratch *scratch;

    // Allocate scratch memory.
    scratch = malloc(sizeof(dopri87_scratch));
    for (k_i = 0; k_i < dopri87_stages; k_i++) {
        scratch->k[k_i] = malloc(system->dims * sizeof(arpra_range));
    }
    scratch->x_new_8 = malloc(system->dims * sizeof(arpra_range));
    scratch->x_new_7 = malloc(system->dims * sizeof(arpra_range));

    // Initialise scratch memory.
    prec_internal = arpra_get_internal_precision();
    for (x_i = 0; x_i < system->dims; x_i++) {
        prec_x = arpra_get_precision(&(system->x[x_i]));
        for (k_i = 0; k_i < dopri87_stages; k_i++) {
            arpra_init2(&(scratch->k[k_i][x_i]), prec_x);
        }
        arpra_init2(&(scratch->x_new_8[x_i]), prec_x);
        arpra_init2(&(scratch->x_new_7[x_i]), prec_x);
    }
    for (k_i = 0; k_i < dopri87_stages; k_i++) {
        scratch->a[k_i] = &(scratch->a_[(k_i * (k_i - 1)) / 2]);
        scratch->ah[k_i] = &(scratch->ah_[(k_i * (k_i - 1)) / 2]);
        for (k_j = 0; k_j < k_i; k_j++) {
            arpra_init2(&(scratch->a[k_i][k_j]), prec_internal);
            arpra_init2(&(scratch->ah[k_i][k_j]), prec_internal);
        }
        arpra_init2(&(scratch->b_8[k_i]), prec_internal);
        arpra_init2(&(scratch->bh_8[k_i]), prec_internal);
        arpra_init2(&(scratch->b_7[k_i]), prec_internal);
        arpra_init2(&(scratch->bh_7[k_i]), prec_internal);
        arpra_init2(&(scratch->c[k_i]), prec_internal);
        arpra_init2(&(scratch->ch[k_i]), prec_internal);
        arpra_init2(&(scratch->temp_t[k_i]), prec_internal);
    }
    arpra_init2(&(scratch->temp_x), prec_internal);

    // Precompute constants.
    dopri87_compute_constants(stepper, prec_internal);

    // Set stepper parameters.
    stepper->method = arpra_ode_dopri87;
    stepper->system = system;
    stepper->error = NULL;
    stepper->scratch = scratch;
}

static void dopri87_clear (arpra_ode_stepper *stepper)
{
    arpra_uint x_i, k_i, k_j;
    arpra_ode_system *system;
    dopri87_scratch *scratch;

    system = stepper->system;
    scratch = (dopri87_scratch *) stepper->scratch;

    // Clear scratch memory.
    for (x_i = 0; x_i < system->dims; x_i++) {
        for (k_i = 0; k_i < dopri87_stages; k_i++) {
            arpra_clear(&(scratch->k[k_i][x_i]));
        }
        arpra_clear(&(scratch->x_new_8[x_i]));
        arpra_clear(&(scratch->x_new_7[x_i]));
    }
    for (k_i = 0; k_i < dopri87_stages; k_i++) {
        for (k_j = 0; k_j < k_i; k_j++) {
            arpra_clear(&(scratch->a[k_i][k_j]));
            arpra_clear(&(scratch->ah[k_i][k_j]));
        }
        arpra_clear(&(scratch->b_8[k_i]));
        arpra_clear(&(scratch->bh_8[k_i]));
        arpra_clear(&(scratch->b_7[k_i]));
        arpra_clear(&(scratch->bh_7[k_i]));
        arpra_clear(&(scratch->c[k_i]));
        arpra_clear(&(scratch->ch[k_i]));
        arpra_clear(&(scratch->temp_t[k_i]));
    }
    arpra_clear(&(scratch->temp_x));

    // Free scratch memory.
    for (k_i = 0; k_i < dopri87_stages; k_i++) {
        free(scratch->k[k_i]);
    }
    free(scratch->x_new_8);
    free(scratch->x_new_7);
    free(scratch);
}

static void dopri87_step (arpra_ode_stepper *stepper, const arpra_range *h)
{
    arpra_uint x_i, k_i, k_j;
    arpra_precision prec_t, prec_x;
    arpra_range *x_new;
    arpra_ode_system *system;
    dopri87_scratch *scratch;

    system = stepper->system;
    scratch = (dopri87_scratch *) stepper->scratch;

    // Synchronise scratch precision and prepare step parameters.
    prec_t = arpra_get_precision(system->t);
    for (x_i = 0; x_i < system->dims; x_i++) {
        prec_x = arpra_get_precision(&(system->x[x_i]));
        for (k_i = 0; k_i < dopri87_stages; k_i++) {
            arpra_set_precision(&(scratch->k[k_i][x_i]), prec_x);
        }
        arpra_set_precision(&(scratch->x_new_8[x_i]), prec_x);
        arpra_set_precision(&(scratch->x_new_7[x_i]), prec_x);
    }
    for (k_i = 0; k_i < dopri87_stages; k_i++) {
        for (k_j = 0; k_j < k_i; k_j++) {
            arpra_set_precision(&(scratch->ah[k_i][k_j]), prec_t);
            arpra_mul(&(scratch->ah[k_i][k_j]), &(scratch->a[k_i][k_j]), h);
        }
        arpra_set_precision(&(scratch->bh_8[k_i]), prec_t);
        arpra_mul(&(scratch->bh_8[k_i]), &(scratch->b_8[k_i]), h);
        arpra_set_precision(&(scratch->bh_7[k_i]), prec_t);
        arpra_mul(&(scratch->bh_7[k_i]), &(scratch->b_7[k_i]), h);
        arpra_set_precision(&(scratch->ch[k_i]), prec_t);
        arpra_mul(&(scratch->ch[k_i]), &(scratch->c[k_i]), h);
        arpra_set_precision(&(scratch->temp_t[k_i]), prec_t);
        arpra_add(&(scratch->temp_t[k_i]), system->t, &(scratch->ch[k_i]));
    }

    // Begin step.
    for (x_i = 0; x_i < system->dims; x_i++) {
        prec_x = arpra_get_precision(&(system->x[x_i]));
        arpra_set_precision(&(scratch->temp_x), prec_x);

        // Compute k stages.
        for (k_i = 0; k_i < dopri87_stages; k_i++) {
            for (k_j = 0; k_j < k_i; k_j++) {
                x_new = (k_j == 0) ? system->x : scratch->x_new_8;
                arpra_mul(&(scratch->temp_x), &(scratch->ah[k_i][k_j]), &(scratch->k[k_j][x_i]));
                arpra_add(&(scratch->x_new_8[x_i]), &(x_new[x_i]), &(scratch->temp_x));
            }
            system->f(scratch->k[k_i],
                      &(scratch->temp_t[k_i]), scratch->x_new_8,
                      x_i, system->params);
        }

        // Compute fourth-order approximation.
        for (k_j = 0; k_j < dopri87_stages; k_j++) {
            x_new = (k_j == 0) ? system->x : scratch->x_new_7;
            arpra_mul(&(scratch->temp_x), &(scratch->bh_7[k_j]), &(scratch->k[k_j][x_i]));
            arpra_add(&(scratch->x_new_7[x_i]), &(x_new[x_i]), &(scratch->temp_x));
        } 
    }

    // Advance system.
    arpra_add(system->t, system->t, h);
    x_new = scratch->x_new_8;
    scratch->x_new_8 = system->x;
    system->x = x_new;
}

static const arpra_ode_method dopri87 =
{
    .init = &dopri87_init,
    .clear = &dopri87_clear,
    .step = &dopri87_step,
    .stages = dopri87_stages,
};

const arpra_ode_method *arpra_ode_dopri87 = &dopri87;
