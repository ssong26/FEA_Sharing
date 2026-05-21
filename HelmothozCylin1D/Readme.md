# 1D Cylindrical Helmholtz FEM

A compact Python implementation of a one-dimensional finite element solver for a cylindrical Helmholtz-type boundary value problem.

This example solves a radially symmetric PDE in cylindrical coordinates and studies the convergence of the numerical solution in both $L^2$ and $H^1$ norms.

---

## Problem

The governing equation is:

$$
\frac{1}{r}\frac{d}{dr}\left(r\frac{du}{dr}\right) - u = 0
$$

on the domain:

$$
0 \le r \le 1
$$

with boundary conditions:

```text
u(1) = 1
du/dr = 0 at r = 0
```

The exact solution is:

$$
u(r) = \frac{I_0(r)}{I_0(1)}
$$

where $I_0$ is the modified Bessel function of the first kind.

---

## Method

The code uses:

- 1D linear finite elements
- uniform radial mesh
- direct stiffness matrix assembly
- Dirichlet boundary condition enforcement
- comparison with the analytical Bessel-function solution
- convergence analysis in $L^2$ and $H^1$ norms

The weak-form functional is based on:

$$
\frac{1}{2}\int_0^1 r\left[\left(\frac{du}{dr}\right)^2 + u^2\right]dr
$$

---

## Files

```text
pde_hm3_fem.py
```

Main Python implementation.

```text
Homework 03 PDE.pdf
```

Original problem statement / derivation note.

```text
Readme.md
```

Documentation for this example.

---

## Requirements

```bash
pip install numpy scipy matplotlib
```

---

## Usage

Run:

```bash
python pde_hm3_fem.py
```

The script computes the FEM solution for different mesh sizes and saves the convergence plot as:

```text
pde.pdf
```

---

## Output

The generated plot shows:

```text
ln(number of elements) vs ln(error)
```

for:

- $L^2$ error
- $H^1$ error

This provides a simple convergence check for the finite element implementation.

---

## Topics

- Finite Element Method
- Cylindrical coordinates
- Helmholtz equation
- Bessel functions
- Error analysis
- Numerical PDEs

---

## Author

Siyuan Song
