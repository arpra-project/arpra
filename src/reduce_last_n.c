/*
 * reduce_last_n.c -- Reduce the last n deviation terms.
 *
 * Copyright 2017-2018 James Paul Turner.
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

void arpra_reduce_last_n (arpra_range *z, arpra_uint n)
{
    arpra_uint zTerm, zNext;
    arpra_mpfr **summands;
    arpra_mpfr temp;
    arpra_prec prec_internal;

    // Handle trivial cases.
    if (n > z->nTerms) n = z->nTerms;
    if (n < 2) return;

    // Handle domain violations.
    if (arpra_nan_p(z)) return;
    if (arpra_inf_p(z)) return;

    // Initialise vars.
    prec_internal = arpra_get_internal_precision();
    mpfr_init2(&temp, prec_internal);
    mpfr_set_prec(&(z->radius), prec_internal);
    mpfr_set_ui(&(z->radius), 0, MPFR_RNDU);
    zTerm = z->nTerms - n;
    summands = malloc(n * sizeof(arpra_mpfr *));

    // Merge the last n deviation terms.
    for (zNext = zTerm; zNext < z->nTerms; zNext++) {
        mpfr_abs(&(z->deviations[zNext]), &(z->deviations[zNext]), MPFR_RNDN);
        summands[zNext - zTerm] = &(z->deviations[zNext]);
    }
    mpfr_sum(&(z->deviations[zTerm]), summands, n, MPFR_RNDU);

    // Add the remaining deviation terms to radius.
    for (zNext = 0; zNext < zTerm; zNext++) {
        mpfr_abs(&temp, &(z->deviations[zNext]), MPFR_RNDU);
        mpfr_add(&(z->radius), &(z->radius), &temp, MPFR_RNDU);
    }

    // Store nonzero merged deviation term.
    if (!mpfr_zero_p(&(z->deviations[zTerm]))) {
        z->symbols[zTerm] = arpra_next_symbol();
        mpfr_add(&(z->radius), &(z->radius), &(z->deviations[zTerm]), MPFR_RNDU);
        zTerm++;
    }

    // Clear unused deviation terms.
    for (zNext = zTerm; zNext < z->nTerms; zNext++) {
        mpfr_clear(&(z->deviations[zNext]));
    }

    // Handle domain violations, and resize memory.
    z->nTerms = zTerm;
    if (mpfr_nan_p(&(z->radius))) {
        arpra_set_nan(z);
    }
    else if (mpfr_inf_p(&(z->radius))) {
        arpra_set_inf(z);
    }
    else {
        if (z->nTerms == 0) {
            free(z->symbols);
            free(z->deviations);
        }
        else {
            z->symbols = realloc(z->symbols, z->nTerms * sizeof(arpra_uint));
            z->deviations = realloc(z->deviations, z->nTerms * sizeof(arpra_mpfr));
        }
    }

    // Clear vars.
    mpfr_clear(&temp);
    free(summands);
}






/*
void new_arpra_reduce_last_n (arpra_range *z, const arpra_range *x, const arpra_uint n)
{
    arpra_uint xTerm, zTerm, unmerged;
    arpra_mpfr temp, error;
    arpra_prec prec, prec_internal;
    arpra_range zNew;

    // Handle trivial cases.
    if (n < 2) return;

    // Domain violations:
    // NaN  =  Nan
    // Inf  =  Inf

    // Handle domain violations.
    if (arpra_nan_p(x)) {
        arpra_set_nan(z);
        return;
    }
    if (arpra_inf_p(x)) {
        arpra_set_inf(z);
        return;
    }

    // Initialise vars.
    prec = arpra_get_precision(z);
    prec_internal = arpra_get_internal_precision();
    mpfr_init2(&temp, prec_internal);
    mpfr_init2(&error, prec_internal);
    mpfr_init2(&(zNew.centre), prec);
    mpfr_init2(&(zNew.radius), prec_internal);
    mpfr_set_si(&error, 0, MPFR_RNDU);
    unmerged = (n < x->nTerms) ? (x->nTerms - n) : 0;

    // z_0 = x_0
    if (mpfr_set(&(zNew.centre), &(x->centre), MPFR_RNDN)) {
        arpra_helper_error_half_ulp(&temp, &(zNew.centre));
        mpfr_add(&error, &error, &temp, MPFR_RNDU);
    }

    // Allocate memory for all possible deviation terms.
    zNew.nTerms = unmerged + 1;
    zNew.symbols = malloc(zNew.nTerms * sizeof(arpra_uint));
    zNew.deviations = malloc(zNew.nTerms * sizeof(arpra_mpfr));

    for (xTerm = 0, zTerm = 0; xTerm < unmerged; xTerm++) {
        zNew.symbols[zTerm] = x->symbols[xTerm];
        mpfr_init2(&(zNew.deviations[zTerm]), prec);

        // z_i = x_i

        // ============= WE DONT NEED TO SAVE SIGN UNTIL RADIUS CALCULATION =======

        if (mpfr_set(&(zNew.deviations[zTerm]), &(x->deviations[xTerm]), MPFR_RNDN)) {
            arpra_helper_error_half_ulp(&temp, &(zNew.deviations[zTerm]));
            mpfr_add(&error, &error, &temp, MPFR_RNDU);
        }

        // Store nonzero deviation terms.
        if (mpfr_zero_p(&(zNew.deviations[zTerm]))) {
            mpfr_clear(&(zNew.deviations[zTerm]));
        }
        else {


            // x cannot be used because x is const, and may still be needed
            // z cannot be used because z may have fewer terms than x

            // znew can be used because we can save abs x, then restore signs from x

            // but znew does not have room to store merge terms


            zTerm++;
        }

        // SET ZNEW AND DO ABS OF X
    }



    // Store nonzero reduced deviation term.
    if (!mpfr_zero_p(&(z->deviations[zTerm]))) {
        z->symbols[zTerm] = arpra_next_symbol();
        zTerm++;
    }

    // Clear unused deviation terms.
    //for (zNext = zTerm; zNext < z->nTerms; zNext++) {
    //    mpfr_clear(&(z->deviations[zNext]));
    //}

    // Handle domain violations, and resize memory.
    z->nTerms = zTerm;
    if (mpfr_nan_p(&(z->radius))) {
        arpra_set_nan(z);
    }
    else if (mpfr_inf_p(&(z->radius))) {
        arpra_set_inf(z);
    }
    else {
        if (z->nTerms == 0) {
            free(z->symbols);
            free(z->deviations);
        }
        else {
            z->symbols = realloc(z->symbols, z->nTerms * sizeof(arpra_uint));
            z->deviations = realloc(z->deviations, z->nTerms * sizeof(arpra_mpfr));
        }
    }

    // Compute new radius.
    arpra_helper_radius(z);
}
*/
