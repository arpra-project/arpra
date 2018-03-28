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

// Arpra ODE typedefs.
typedef struct arpra_ode_system_struct arpra_ode_system;
typedef struct arpra_ode_stepper_struct arpra_ode_stepper;
typedef struct arpra_ode_method_struct arpra_ode_method;

// System definition.
struct arpra_ode_system_struct
{
    void (*function) (arpra_range *dxdt,
                      const arpra_range *x, const arpra_range *t,
                      arpra_uint dimensions, void *parameters);
    const void *parameters;
    arpra_range *state;
    arpra_uint dimensions;
};

// Stepper definition.
struct arpra_ode_stepper_struct
{
    const arpra_ode_method *method;
    void *workspace;
    arpra_uint dimensions;
};

// Step method definition.
struct arpra_ode_method_struct
{
    void (*init) (arpra_ode_stepper *stepper,
                  const arpra_ode_system *system);
    void (*init2) (arpra_ode_stepper *stepper,
                   const arpra_ode_system *system);
    void (*clear) (arpra_ode_stepper *stepper);
    void (*reset) (arpra_ode_stepper *stepper);
    void (*step) (arpra_ode_stepper *stepper,
                  arpra_range *dx, arpra_range *t, arpra_range *error,
                  arpra_range *x, const arpra_range *h);
};

// Stepper functions.
void arpra_ode_stepper_init (arpra_ode_stepper *stepper);
void arpra_ode_stepper_init2 (arpra_ode_stepper *stepper);
void arpra_ode_stepper_clear (arpra_ode_stepper *stepper);
void arpra_ode_stepper_reset (arpra_ode_stepper *stepper);
void arpra_ode_stepper_step (const arpra_ode_stepper *stepper,
                             arpra_range *dx, arpra_range *t, arpra_range *error,
                             arpra_range *x, const arpra_range *h);

// Arpra built-in step methods.
extern const arpra_ode_method *arpra_ode_euler;
extern const arpra_ode_method *arpra_ode_rk2;

#endif // ARPRA_ODE_H
