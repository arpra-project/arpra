/*
 * div.c
 *
 *  Created on: 16 Nov 2016
 *      Author: jt273
 */

#include "mpfa.h"

/*
 * For now we just multiply the numerator with the reciprocal of the denominator.
 *
 * TODO: find a better way to do mpfa_div
 */

void mpfa_div (mpfa_ptr z, mpfa_srcptr x, mpfa_srcptr y) {
	mpfa_t zNew;

	mpfa_init(zNew);
	mpfa_inv(zNew, y);
	mpfa_mul(z, x, zNew);
	mpfa_clear(zNew);
}
