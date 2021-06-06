# Symbolic Constructor for Lagrangian Polynomial and Finite Difference
a python script to symbolically construct Lagrange polynomials and finite difference equations

# Usage
donwload lagrangianpoly.py and import it like below:

```Python
from lagrangianpoly import LagrangianBasis, LagrangianPoly, Derivative
```

# Examples

## Lagrangian Basis

```Python
    x = sp.symbols('x')
    eq = LagrangianBasis(x, degreeOfPolynomial=5, pointAt=0)
```

```math
\frac{(x - x_1)(x - x_2)(x - x_3)(x - x_4)(x - x_5)}{(x_0 - x_1)(x_0 - x_2)(x_0 - x_3)(x_0 - x_4)(x_0 - x_5)}
```

## Lagrangian Polynomial

Given a set of data \((x_0, f_0), (x_1, f_1), (x_2, f_2)\) where \((x_0, x_1, x_2) = (-\varDelta x, 0, \varDelta x)\), Lagrangian Polynomial at \(x = 0\), e.g. \(x = x_1\), is constructed by the following script:

```Python
    x = sp.symbols('x')
    dx = sp.symbols('dx')

    # 3-point stencil lagrangian polynomial
    x_set = [-dx, 0, dx]
    f_set = sp.symbols('f0:{:d}'.format(len(x_set)))
    eq = LagrangianPoly(x, x_set, f_set)
```

We get `f0*x*(-dx + x)/(2*dx**2) - f1*(-dx + x)*(dx + x)/dx**2 + f2*x*(dx + x)/(2*dx**2)` corresponding to

```math
f_0\frac{x(x -\varDelta x)}{2\varDelta x^2} - f_1\frac{(x - \varDelta x)(x+\varDelta x)}{\varDelta x^2} + f_2\frac{x(\varDelta x + x)}{2\varDelta x^2}
```

## Finite Difference Equation
### 3-point stencil 1st order central difference on the regular grid

Given a set of data \((x_0, f_0), (x_1, f_1), (x_2, f_2)\) where \((x_0, x_1, x_2) = (-\varDelta x, 0, \varDelta x)\), 3-point stencil 1st order central difference on the regular grid at \(x = 0\), e.g. \(x = x_1\), is constructed by the following script:

```Python
    x_set = [-dx, 0, dx]
    f_set = sp.symbols('f0:{:d}'.format(len(x_set)))
    eq = Derivative(LagrangianPoly(x, x_set, f_set), x, orderOfDifference=1)
```

We get `(-f0 + f2)/(2*dx)` corresponding to well known expression

```math
\frac{f_2-f_0}{2\varDelta x}
```

### 3-point stencil 1st order central difference on the staggered grid

Given a set of data \((x_0, f_0), (x_1, f_1), (x_2, f_2)\) where \((x_0, x_1, x_2) = (\frac{-\varDelta x}{2}, 0, \frac{\varDelta x}{2})\), 3-point stencil 1st order central difference on the staggered grid at \(x = 0\), e.g. \(x = x_1\), is constructed by the following script:

```Python
    x_set = [-dx/2, 0, dx/2]
    f_set = sp.symbols('f0:{:d}'.format(len(x_set)))
    eq = Derivative(LagrangianPoly(x, x_set, f_set), x, orderOfDifference=1)
```

We get `(-f0 + f2)/dx`.
