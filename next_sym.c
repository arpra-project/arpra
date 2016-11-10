/*
 * next_symbol.c
 *
 *  Created on: 2 Nov 2016
 *      Author: jt273
 */

#include "mpfa.h"

unsigned long mpfa_next_sym() {
	static unsigned long mpfa_symbol;
	return mpfa_symbol++;
}
