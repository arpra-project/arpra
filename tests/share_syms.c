/*
 * share_syms.c -- Share noise symbols between two affine forms.
 *
 * Copyright 2017-2020 James Paul Turner.
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

#include "arpra-test.h"

void test_share_all_syms (arpra_range *x1, arpra_range *x2)
{

    arpra_uint symbol, i;
    arpra_int x1_has_next, x2_has_next;

    i = 0;
    x1_has_next = x1->nTerms > 0;
    x2_has_next = x2->nTerms > 0;
    while (x1_has_next || x2_has_next) {
        symbol = arpra_helper_next_symbol();

        // Set x1 and x2 symbol if they exist.
        if (x1_has_next) {
            x1->symbols[i] = symbol;
        }
        if (x2_has_next) {
            x2->symbols[i] = symbol;
        }

        i++;
        x1_has_next = i < x1->nTerms;
        x2_has_next = i < x2->nTerms;
    }
}

void test_share_rand_syms (arpra_range *x1, arpra_range *x2)
{
    arpra_uint symbol, i;
    arpra_int x1_has_next, x2_has_next;

    i = 0;
    x1_has_next = x1->nTerms > 0;
    x2_has_next = x2->nTerms > 0;

    while (x1_has_next || x2_has_next) {
        symbol = arpra_helper_next_symbol();

        // Set x1 and x2 symbol if they exist.
        if (!x2_has_next) {
            x1->symbols[i] = symbol;
        }
        else if (!x1_has_next) {
            x2->symbols[i] = symbol;
        }

        // Else randomly share x1 and x2 symbols.
        else {
            if (gmp_urandomb_ui(test_randstate, 1)) {
                x1->symbols[i] = symbol;
                x2->symbols[i] = symbol;
            }
            else {
                x1->symbols[i] = symbol;
                x2->symbols[i] = arpra_helper_next_symbol();
            }
        }

        i++;
        x1_has_next = i < x1->nTerms;
        x2_has_next = i < x2->nTerms;
    }
}

void test_share_n_syms (arpra_range *x1, arpra_range *x2, arpra_uint n)
{
    arpra_uint symbol, i;
    arpra_int x1_has_next, x2_has_next;

    i = 0;
    x1_has_next = x1->nTerms > 0;
    x2_has_next = x2->nTerms > 0;

    while (x1_has_next || x2_has_next) {
        symbol = arpra_helper_next_symbol();

        // Set x1 and x2 symbol if they exist.
        if (!x2_has_next) {
            x1->symbols[i] = symbol;
        }
        else if (!x1_has_next) {
            x2->symbols[i] = symbol;
        }

        // Else share the first n symbols in x1 and x2.
        else {
            if (n > 0) {
                x1->symbols[i] = symbol;
                x2->symbols[i] = symbol;
                n--;
            }
            else {
                x1->symbols[i] = symbol;
                x2->symbols[i] = arpra_helper_next_symbol();
            }
        }

        i++;
        x1_has_next = i < x1->nTerms;
        x2_has_next = i < x2->nTerms;
    }
}
