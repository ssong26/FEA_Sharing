# Abaqus UEL — 2D Beam Element Implementation

This folder contains a custom Abaqus `UEL` (User Element) implementation for a 2D beam element.

The project demonstrates how to develop user-defined finite elements in Abaqus using Fortran, including:

- element stiffness matrix construction
- numerical integration
- beam kinematics
- assembly of element residuals and tangent stiffness
- interaction with the Abaqus solver through the `UEL` interface

---

# Files

```text
UEL_2DBeam.for
```

Main Abaqus user element (`UEL`) implementation.

```text
*.inp
```

Example Abaqus input files for testing the custom beam element.

```text
Readme.md
```

Documentation for the project.

---

# Overview

The implementation develops a 2D beam finite element within the Abaqus user-element framework.

The element supports:

- axial deformation
- bending deformation
- nodal displacement DOFs
- element stiffness assembly
- numerical integration at Gauss points

The purpose of the project is to explore the internal structure of finite element formulations beyond built-in Abaqus elements.

---

# Beam Theory

The formulation is based on classical beam mechanics.

The elemental response includes:

- axial strain energy
- bending strain energy

The general beam stiffness contribution can be written as:

$$
K_e =
\int_{\Omega}
B^T D B \, d\Omega
$$

where:

- $B$ is the strain-displacement matrix
- $D$ is the constitutive matrix
- $\Omega$ is the element domain

---

# Degrees of Freedom

The beam element contains nodal degrees of freedom associated with:

- horizontal displacement
- vertical displacement
- rotational displacement

Typical nodal DOFs:

$$
[u_x,\; u_y,\; \theta]
$$

for each node.

---

# Numerical Method

The implementation includes:

- shape function evaluation
- Jacobian mapping
- Gaussian integration
- element residual assembly
- tangent stiffness matrix construction

The `UEL` subroutine communicates with Abaqus through:

- `RHS`
- `AMATRX`
- `SVARS`
- `ENERGY`

and other standard Abaqus user-element arrays.

---

# Abaqus User Element Workflow

The project demonstrates the standard workflow for developing user-defined finite elements in Abaqus:

1. Define custom element topology
2. Implement stiffness and residual calculation
3. Compile user subroutine
4. Link with Abaqus solver
5. Run custom-element simulations

---

# Running the Simulation

Typical Abaqus execution command:

```bash
abaqus job=beam_test user=UEL_2DBeam.for
```

The Abaqus input file should include:

```abaqus
*User Element
```

and corresponding custom element definitions.

---

# Educational Purpose

This project was developed for educational purposes in:

- finite element analysis
- computational mechanics
- custom element development
- Abaqus user subroutines
- numerical implementation of beam theory

The implementation prioritizes transparency and readability over industrial-level optimization.

---

# Topics Covered

- Abaqus UEL
- beam finite elements
- custom FEM implementation
- finite element stiffness assembly
- Gaussian integration
- computational solid mechanics

---

# Possible Extensions

Potential future improvements include:

- Timoshenko beam formulation
- geometric nonlinearity
- large rotation beam theory
- dynamic beam analysis
- mass matrix implementation
- nonlinear material behavior
- 3D beam elements
- higher-order interpolation

---

# Author

Siyuan Song  
Brown University  
Computational Mechanics / Scientific Computing
