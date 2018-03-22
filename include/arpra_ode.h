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
typedef struct __arpra_ode_system_struct arpra_ode_system_t;
typedef struct __arpra_ode_system_struct *arpra_ode_system_ptr;
typedef const struct __arpra_ode_system_struct *arpra_ode_system_srcptr;

struct __arpra_ode_system_struct
{
    void (*function) (arpra_ptr dxdt, arpra_ptr t, arpra_srcptr y, void *parameters);
    void (*jacobian) (arpra_ptr dfdy, arpra_ptr dfdt, arpra_ptr t, arpra_srcptr y, void *parameters);
    void *parameters;
};

// ODE stepper definition.
struct __arpra_ode_stepper_struct
{
    arpra_ode_system_srcptr system;
    void *state;

    void (*init) (struct __arpra_ode_stepper_struct *stepper,
                  const struct __arpra_ode_system_struct *system);
    void (*init2) (struct __arpra_ode_stepper_struct *stepper,
                   const struct __arpra_ode_system_struct *system);
    void (*clear) (struct __arpra_ode_stepper_struct *stepper);
    void (*reset) (struct __arpra_ode_stepper_struct *stepper);
    void (*step) (struct __arpra_ode_stepper_struct *stepper,
                  arpra_ptr dx, arpra_ptr t, arpra_ptr error,
                  arpra_srcptr x, arpra_srcptr h);
};

typedef struct __arpra_ode_stepper_struct arpra_ode_stepper_t;
typedef struct __arpra_ode_stepper_struct *arpra_ode_stepper_ptr;
typedef const struct __arpra_ode_stepper_struct *arpra_ode_stepper_srcptr;

// ODE stepper functions.
void arpra_ode_stepper_init (arpra_ode_stepper_ptr stepper);
void arpra_ode_stepper_init2 (arpra_ode_stepper_ptr stepper);
void arpra_ode_stepper_clear (arpra_ode_stepper_ptr stepper);
void arpra_ode_stepper_reset (arpra_ode_stepper_ptr stepper);
void arpra_ode_stepper_step (arpra_ode_stepper_srcptr stepper,
                             arpra_ptr dx, arpra_ptr t, arpra_ptr error,
                             arpra_srcptr x, arpra_srcptr h);

// Built-in ODE steppers.
extern arpra_ode_stepper_srcptr arpra_ode_euler;
extern arpra_ode_stepper_srcptr arpra_ode_rk2;

#endif // ARPRA_ODE_H
