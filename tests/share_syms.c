/*
 * share_syms.c -- Randomly share noise symbols between two affine forms.
 *
 * Copyright 2017 James Paul Turner.
 *
 * This file is part of the MPFA library.
 *
 * The MPFA library is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * The MPFA library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with the MPFA library. If not, see <http://www.gnu.org/licenses/>.
 */

#include "mpfa-test.h"

void test_share_syms (mpfa_ptr x, mpfa_ptr y, const mpfa_uint_t share_chance)
{
    mpfa_uint_t symbol, term;
    mpfa_int_t x_has_next, y_has_next;

    // Symbol share chance should be from 0 to 9.
    if (share_chance > 9) {
        fprintf(stderr, "Error: symbol share chance should be 0 to 9");
        exit(EXIT_FAILURE);
    }

    term = 0;
    x_has_next = x->nTerms > 0;
    y_has_next = y->nTerms > 0;
    while (x_has_next || y_has_next) {
        symbol = mpfa_next_sym();

        // Set x symbol if y has no more terms.
        if (!y_has_next) {
            x->symbols[term] = symbol;
            x_has_next = ++term < x->nTerms;
        }

        // Set y symbol if x has no more terms.
        else if (!x_has_next) {
            y->symbols[term] = symbol;
            y_has_next = ++term < y->nTerms;
        }

        else {
            if (test_rand_ui(10) < share_chance) {
                // x and y share this symbol.
                x->symbols[term] = symbol;
                y->symbols[term] = symbol;
            }
            else {
                // x and y have different symbols.
                x->symbols[term] = symbol;
                y->symbols[term] = mpfa_next_sym();
            }

            term++;
            x_has_next = term < x->nTerms;
            y_has_next = term < y->nTerms;
        }
    }
}
