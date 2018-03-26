/*
 * arpra_ode.h -- Arpra public header for ordinary differential equations.
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

#ifndef ARPRA_ODE_H
#define ARPRA_ODE_H

#include <arpra.h>

// ODE system definition.
struct arpra_ode_system
{
    void (*function) (struct arpra_range *dxdt, struct arpra_range *t,
                      const struct arpra_range *y, void *parameters);
    void (*jacobian) (struct arpra_range *dfdy, struct arpra_range *dfdt, struct arpra_range *t,
                      const struct arpra_range *y, void *parameters);
    void *parameters;
};

// ODE stepper definition.
struct arpra_ode_stepper
{
    const struct arpra_ode_system *system;
    void *state;

    void (*init) (struct arpra_ode_stepper *stepper,
                  const struct arpra_ode_system *system);
    void (*init2) (struct arpra_ode_stepper *stepper,
                   const struct arpra_ode_system *system);
    void (*clear) (struct arpra_ode_stepper *stepper);
    void (*reset) (struct arpra_ode_stepper *stepper);
    void (*step) (struct arpra_ode_stepper *stepper,
                  struct arpra_range *dx, struct arpra_range *t, struct arpra_range *error,
                  struct arpra_range *x, const struct arpra_range *h);
};

// ODE stepper functions.
void arpra_ode_stepper_init (struct arpra_ode_stepper *stepper);
void arpra_ode_stepper_init2 (struct arpra_ode_stepper *stepper);
void arpra_ode_stepper_clear (struct arpra_ode_stepper *stepper);
void arpra_ode_stepper_reset (struct arpra_ode_stepper *stepper);
void arpra_ode_stepper_step (const struct arpra_ode_stepper *stepper,
                             struct arpra_range *dx, struct arpra_range *t, struct arpra_range *error,
                             struct arpra_range *x, const struct arpra_range *h);

// Built-in ODE steppers.
extern const struct arpra_ode_stepper *arpra_ode_euler;
extern const struct arpra_ode_stepper *arpra_ode_rk2;

#endif // ARPRA_ODE_H
