/*
 * henon_map_mpfi.c -- Test Henon map model.
 *
 * Copyright 2020 James Paul Turner.
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

#include <stdlib.h>
#include <stdio.h>
#include <mpfi.h>

int main (int argc, char *argv[])
{
    mpfi_t one, x_new, y_new;
    mpfi_t a, b, x, y;
    mpfr_t uncertainty;
    FILE *x_out, *y_out;
    mpfr_prec_t prec, prec_internal;
    unsigned long n, i;

    n = 500;
    prec = 53;
    prec_internal = 128;

    // Initialise MPFI ranges
    mpfi_init2(one, 2);
    mpfi_init2(x_new, prec);
    mpfi_init2(y_new, prec);
    mpfi_init2(a, prec);
    mpfi_init2(b, prec);
    mpfi_init2(x, prec);
    mpfi_init2(y, prec);

    // Initialise MPFR vars
    mpfr_init2(uncertainty, prec_internal);

    // Set MPFI ranges (almost chaotic)
    mpfi_set_d(one, 1.0);
    mpfi_set_str(a, "1.057", 10);
    //mpfi_set_str(a, argv[1], 10);
    mpfi_set_str(b, "0.3", 10);
    mpfi_set_d(x, 0.0);
    mpfi_set_d(y, 0.0);
    mpfr_set_str(uncertainty, "1e-5", 10, MPFR_RNDU);
    mpfr_sub(&(x->left), &(x->left), uncertainty, MPFR_RNDD);
    mpfr_add(&(x->right), &(x->right), uncertainty, MPFR_RNDU);
    mpfr_sub(&(y->left), &(y->left), uncertainty, MPFR_RNDD);
    mpfr_add(&(y->right), &(y->right), uncertainty, MPFR_RNDU);

    // Open output files
    x_out = fopen("henon_x.dat", "w");
    y_out = fopen("henon_y.dat", "w");

    // Iterate Henon map
    for (i = 0; i < n; i++) {
        if (i % 10 == 0) {
            printf("%u\n", i);
        }

        // Compute new x
        mpfi_mul(x_new, x, x);
        mpfi_mul(x_new, x_new, a);
        mpfi_sub(x_new, one, x_new);
        mpfi_add(x_new, x_new, y);

        // Compute new y
        mpfi_mul(y_new, b, x);

        // Update x and y
        mpfi_set(x, x_new);
        mpfi_set(y, y_new);

        // Write output
        mpfr_out_str(x_out, 10, 40, &(x->left), MPFR_RNDN);
        fputs(" ", x_out);
        mpfr_out_str(x_out, 10, 40, &(x->right), MPFR_RNDN);
        fputs("\n", x_out);
        mpfr_out_str(y_out, 10, 40, &(y->left), MPFR_RNDN);
        fputs(" ", y_out);
        mpfr_out_str(y_out, 10, 40, &(y->right), MPFR_RNDN);
        fputs("\n", y_out);
    }

    // Clear MPFI ranges
    mpfi_clear(one);
    mpfi_clear(x_new);
    mpfi_clear(y_new);
    mpfi_clear(a);
    mpfi_clear(b);
    mpfi_clear(x);
    mpfi_clear(y);

    // Clear MPFR vars
    mpfr_clear(uncertainty);

    // Close output files
    fclose(x_out);
    fclose(y_out);

    // Cleanup
    mpfr_free_cache();
}
