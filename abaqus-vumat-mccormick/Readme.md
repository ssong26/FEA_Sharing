# Abaqus VUMAT — McCormick-Type Viscoplastic Material Model

This folder contains an Abaqus/Explicit `VUMAT` implementation for a McCormick-type viscoplastic constitutive model.

The code was originally developed for EN2340 finite element analysis coursework. It demonstrates how to implement a user-defined material model in Abaqus/Explicit using Fortran.

---

# Files

```text
abaqus_vumat_McCormick.for
```

A Fortran `VUMAT` subroutine for Abaqus/Explicit.

```text
Readme.md
```

Brief documentation for this example.

---

# Model Overview

The subroutine implements an elastic-viscoplastic material update with internal state variables.

The material response includes:

- isotropic elasticity
- deviatoric stress update
- equivalent von Mises stress calculation
- plastic strain increment update
- Newton-Raphson iteration for the plastic correction
- history-dependent internal variables

The implementation follows the Abaqus/Explicit `VUMAT` interface, where stress and state variables are updated at each material integration point.

---

# State Variables

The implementation stores several history variables:

```text
statev(1) = accumulated plastic strain
statev(2) = internal aging / time-like variable
statev(3) = previous plastic strain increment
```

These variables are updated during the stress integration procedure.

---

# Material Parameters

The material parameters are passed through `props`:

```fortran
E      = props(1)   ! Young's modulus
xnu    = props(2)   ! Poisson's ratio
Y      = props(3)   ! reference yield stress
e0     = props(4)   ! reference strain
m      = props(5)   ! hardening exponent
edot0  = props(6)   ! reference strain rate
S      = props(7)   ! stress scale
H      = props(8)   ! hardening / aging scale
td     = props(9)   ! characteristic time
Omega  = props(10)  ! aging parameter
alfa   = props(11)  ! aging exponent
```

---

# Constitutive Update

The code first computes the elastic trial stress. The equivalent von Mises stress is then evaluated as:

$$
\sigma_e =
\sqrt{\frac{3}{2} \, \mathbf{s} : \mathbf{s}}
$$

where $\mathbf{s}$ is the deviatoric stress tensor.

If the trial response remains elastic, the stress is updated directly. Otherwise, the plastic strain increment is solved iteratively.

---

# Newton-Raphson Plastic Correction

The plastic increment is obtained using a Newton-Raphson iteration.

The nonlinear residual contains contributions from:

- trial equivalent stress
- elastic unloading due to plastic flow
- strain hardening
- rate dependence
- aging / McCormick-type internal variable evolution

The iteration updates the plastic strain increment until the residual reaches the prescribed tolerance.

---

# Abaqus Usage

This file is intended to be compiled and linked with Abaqus/Explicit as a user material subroutine.

A typical command is:

```bash
abaqus job=your_job_name user=abaqus_vumat_McCormick.for
```

The input file should define:

- `*User Material`
- the 11 material constants listed above
- the required number of state variables through `*Depvar`

Example structure:

```abaqus
*User Material, constants=11
E, nu, Y, e0, m, edot0, S, H, td, Omega, alfa

*Depvar
3
```

---

# Educational Purpose

This example is intended for studying:

- Abaqus user material implementation
- explicit dynamics material updates
- viscoplastic constitutive modeling
- Newton-Raphson stress integration
- history-dependent material behavior

The implementation prioritizes transparency of the constitutive update rather than production-level robustness.

---

# Notes

- The code is written for Abaqus/Explicit `VUMAT`.
- The implementation assumes the Abaqus corotational stress convention.
- The subroutine is designed for educational and research purposes.
- Additional testing is recommended before use in production simulations.

---

# Related Topics

- Abaqus VUMAT
- user-defined material subroutines
- viscoplasticity
- McCormick effect
- dynamic strain aging
- nonlinear constitutive modeling
- explicit finite element analysis

---

# Author

Siyuan Song  
Brown University  
Computational Mechanics / Scientific Computing
