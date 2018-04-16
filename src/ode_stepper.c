/*
 * ode_stepper.c -- Initialise, clear and manipulate ODE steppers.
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

void arpra_ode_stepper_init (arpra_ode_stepper *stepper, arpra_ode_system *system,
                             const arpra_ode_method *method)
{
    method->init(stepper, system);
}

void arpra_ode_stepper_clear (arpra_ode_stepper *stepper)
{
    stepper->method->clear(stepper);
}

void arpra_ode_stepper_reset (arpra_ode_stepper *stepper)
{
    stepper->method->reset(stepper);
}

void arpra_ode_stepper_step (arpra_ode_stepper *stepper, const arpra_range *h)
{
    stepper->method->step(stepper, h);
}
