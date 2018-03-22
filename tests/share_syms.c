/*
 * share_syms.c -- Share noise symbols between two affine forms.
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

#include "arpra-test.h"

void test_share_all_syms (arpra_ptr x, arpra_ptr y)
{

    arpra_uint symbol, term;
    arpra_int x_has_next, y_has_next;

    term = 0;
    x_has_next = x->nTerms > 0;
    y_has_next = y->nTerms > 0;
    while (x_has_next || y_has_next) {
        symbol = arpra_next_sym();

        // Set x and y symbol if they exist.
        if (x_has_next) {
            x->symbols[term] = symbol;
        }
        if (y_has_next) {
            y->symbols[term] = symbol;
        }

        term++;
        x_has_next = term < x->nTerms;
        y_has_next = term < y->nTerms;
    }
}

void test_share_rand_syms (arpra_ptr x, arpra_ptr y)
{
    arpra_uint symbol, term;
    arpra_int x_has_next, y_has_next;

    term = 0;
    x_has_next = x->nTerms > 0;
    y_has_next = y->nTerms > 0;
    while (x_has_next || y_has_next) {
        symbol = arpra_next_sym();

        // Set x symbol if y has no more terms.
        if (!y_has_next) {
            x->symbols[term] = symbol;
        }

        // Set y symbol if x has no more terms.
        else if (!x_has_next) {
            y->symbols[term] = symbol;
        }

        // x and y have symbols left.
        else {
            if (gmp_urandomb_ui(test_randstate, 1)) {
                // x and y share this symbol.
                x->symbols[term] = symbol;
                y->symbols[term] = symbol;
            }
            else {
                // x and y have different symbols.
                x->symbols[term] = symbol;
                y->symbols[term] = arpra_next_sym();
            }
        }

        term++;
        x_has_next = term < x->nTerms;
        y_has_next = term < y->nTerms;
    }
}
