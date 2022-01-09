# Symbolic Constructor for Lagrangian Polynomial and Finite Difference
a python script to symbolically construct Lagrange polynomials and finite difference equations

# Usage
donwload lagrangianpoly.py, SymbolicFiniteDifference.py, and TaylorExpansion.py and import those like below:

```Python
from lagrangianpoly import LagrangianBasis, LagrangianPoly, Derivative
import SymbolicFiniteDifference as fd
import TaylorExpansion as te
```

# Todo
- [ ] add comments and docstrings
- [ ] add comprehensive test
- [ ] packaging according to Python manner

# Examples
## lagrangianpoly
### Lagrangian Basis

```Python
    x = sp.symbols('x')
    eq = LagrangianBasis(x, degreeOfPolynomial=5, pointAt=0)
```

<!-- ```math
\frac{(x - x_1)(x - x_2)(x - x_3)(x - x_4)(x - x_5)}{(x_0 - x_1)(x_0 - x_2)(x_0 - x_3)(x_0 - x_4)(x_0 - x_5)}
``` -->
<img src="https://render.githubusercontent.com/render/math?math=\frac{(x - x_1)(x - x_2)(x - x_3)(x - x_4)(x - x_5)}{(x_0 - x_1)(x_0 - x_2)(x_0 - x_3)(x_0 - x_4)(x_0 - x_5)}">

### Lagrangian Polynomial

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

<!-- ```math
f_0\frac{x(x -\varDelta x)}{2\varDelta x^2} - f_1\frac{(x - \varDelta x)(x+\varDelta x)}{\varDelta x^2} + f_2\frac{x(\varDelta x + x)}{2\varDelta x^2}
``` -->
<img src="https://render.githubusercontent.com/render/math?math=f_0\frac{x(x -\varDelta x)}{2\varDelta x^2} - f_1\frac{(x - \varDelta x)(x+\varDelta x)}{\varDelta x^2} + f_2\frac{x(\varDelta x + x)}{2\varDelta x^2}">

### Finite Difference Equation
#### 3-point stencil central difference for 1st derivative on the regular grid

Given a set of data \((x_0, f_0), (x_1, f_1), (x_2, f_2)\) where \((x_0, x_1, x_2) = (-\varDelta x, 0, \varDelta x)\), 3-point stencil 1st order central difference on the regular grid at \(x = 0\), e.g. \(x = x_1\), is constructed by the following script:

```Python
    x_set = [-dx, 0, dx]
    f_set = sp.symbols('f0:{:d}'.format(len(x_set)))
    eq = Derivative(LagrangianPoly(x, x_set, f_set), x, orderOfDifference=1)
```

We get `(-f0 + f2)/(2*dx)` corresponding to well known expression

<!-- ```math
\frac{f_2-f_0}{2\varDelta x}
``` -->
<img src="https://render.githubusercontent.com/render/math?math=\frac{f_2-f_0}{2\varDelta x}">

#### 3-point stencil central difference for 1st derivative on the staggered grid

Given a set of data \((x_0, f_0), (x_1, f_1), (x_2, f_2)\) where \((x_0, x_1, x_2) = (\frac{-\varDelta x}{2}, 0, \frac{\varDelta x}{2})\), 3-point stencil 1st order central difference on the staggered grid at \(x = 0\), e.g. \(x = x_1\), is constructed by the following script:

```Python
    x_set = [-dx/2, 0, dx/2]
    f_set = sp.symbols('f0:{:d}'.format(len(x_set)))
    eq = Derivative(LagrangianPoly(x, x_set, f_set), x, orderOfDifference=1)
```

We get `(-f0 + f2)/dx`.

## Symbolic Finite Difference
### getFiniteDifferenceEquation

#### 5-point central difference for 1st derivative with a default interval symbol
```Python
    stencil = [-2, -1, 0, 1, 2]
    eq = fd.getFiniteDifferenceEquation(stencil, orderOfDifference=1)
```

The result is `(f_0 - 8*f_1 + 8*f_3 - f_4)/(12*h)`.

```Python
    stencil = [-2, -1, 0, 1, 2]
    eq = fd.getFiniteDifferenceEquation(stencil, 1, sameSubscriptsAsStencil=True)
```

`(-8*f_{-1} + f_{-2} + 8*f_{1} - f_{2})/(12*h)`
Sorting based on the order of subscripts does not seem to work.

#### 5-point central difference for 1st derivative with an user-defined interval symbol
```Python
    stencil = [-2, -1, 0, 1, 2]
    eq = fd.getFiniteDifferenceEquation(stencil, orderOfDifference=1, intervalSymbolStr='dx')
```

`(f_0 - 8*f_1 + 8*f_3 - f_4)/(12*dx)`

#### 3-point one-sided difference for 1st derivative with a default interval symbol
```Python
    stencil = [0, 1, 2]
    eq = fd.getFiniteDifferenceEquation(stencil, orderOfDifference=1)
```

`(-3*f_0 + 4*f_1 - f_2)/(2*h)`

#### 4-point one-sided difference for 2nd derivative with a default interval symbol
```Python
    stencil = [0, 1, 2, 3]
    eq = fd.getFiniteDifferenceEquation(stencil, orderOfDifference=2)
```

`(2*f_0 - 5*f_1 + 4*f_2 - f_3)/h**2`

#### 5-point central difference for 1st derivative on the staggered grid
```Python
    stencil = [-1.5, -0.5, 0, 0.5, 1.5]
    eq = fd.getFiniteDifferenceEquation(stencil, 1)
```

The result is `(f_0 - 27*f_1 + 27*f_3 - f_4)/(24*h)`.

```Python
    stencil = [-1.5, -0.5, 0, 0.5, 1.5]
    eq = fd.getFiniteDifferenceEquation(stencil, 1, sameSubscriptsAsStencil=True)
```

`(-27*f_{-0.5} + f_{-1.5} + 27*f_{0.5} - f_{1.5})/(24*h)`
Sorting based on the order of subscripts does not seem to work.

### getFiniteDifferenceCoefficients
#### 5-point central difference for 1st derivative
```Python
    stencil = [-2, -1, 0, 1, 2]
    coef = fd.getFiniteDifferenceCoefficients(stencil, orderOfDifference=1)
```

The result is `[1/12, -2/3, 0, 2/3, -1/12]`.

```Python
    # numerator and denominator
    stencil = [-2, -1, 0, 1, 2]
    numr, denom = fd.getFiniteDifferenceCoefficients(stencil, orderOfDifference=1, as_numr_denom=True)
```

`[1, -8, 0, 8, -1] 12`

#### 7-point central difference for 1st derivative
```Python
    stencil = [-3, -2, -1, 0, 1, 2, 3]
    coef = fd.getFiniteDifferenceCoefficients(stencil, orderOfDifference=1)
```

`[-1/60, 3/20, -3/4, 0, 3/4, -3/20, 1/60]`

```Python
    # numerator and denominator
    stencil = [-3, -2, -1, 0, 1, 2, 3]
    numr, denom = fd.getFiniteDifferenceCoefficients(stencil, orderOfDifference=1, as_numr_denom=True)
```

`[-1, 9, -45, 0, 45, -9, 1] 60`

#### 3-point one-sided difference for 1st derivative
```Python
    stencil = [0, 1, 2]
    coef = fd.getFiniteDifferenceCoefficients(stencil, orderOfDifference=1)
```

`[-3/2, 2, -1/2]`

```Python
    # numerator and denominator
    stencil = [0, 1, 2]
    numr, denom = fd.getFiniteDifferenceCoefficients(stencil, orderOfDifference=1, as_numr_denom=True)
```

`[-3, 4, -1] 2`

#### 4-point one-sided difference for 2nd derivative
```Python
    stencil = [0, 1, 2, 3]
    coef = fd.getFiniteDifferenceCoefficients(stencil, orderOfDifference=2)
```

`[2, -5, 4, -1]`

## TaylorExpansion
### TaylorExpansion

#### f(x+h) around x up to term including 6th order difference

```Python
    h = sp.symbols('h')
    f1 = te.TaylorExpansion(h, n=6)
```
`f + f^(1)*h + f^(2)*h**2/2 + f^(3)*h**3/6 + f^(4)*h**4/24 + f^(5)*h**5/120 + f^(6)*h**6/720 `

#### f(x-h) around x up to term including 7th order difference

```Python
    h = sp.symbols('h')
    f_1 = te.TaylorExpansion(-h, n=7)
```

`f - f^(1)*h + f^(2)*h**2/2 - f^(3)*h**3/6 + f^(4)*h**4/24 - f^(5)*h**5/120 + f^(6)*h**6/720 - f^(7)*h**7/5040`

#### f(x+2h) around x up to term including 4th order difference

```Python
    h = sp.symbols('h')
    f2 = te.TaylorExpansion(2*h, n=4)
```
`f + 2*f^(1)*h + 2*f^(2)*h**2 + 4*f^(3)*h**3/3 + 2*f^(4)*h**4/3`

### getTrunctaionError
#### 3-point central finite difference for 1st derivative
```Python
    stencil = [-1, 0, 1]
    err = te.getTruncationError(stencil, 1)
```

`-f^(3)*h**2/6`

#### 3-point central finite difference for 2nd derivative
```Python
    stencil = [-1, 0, 1]
    err = te.getTruncationError(stencil, 2)
```

`-f^(4)*h**2/12`

#### 5-point central finite difference for 1st derivative
```Python
    stencil = [-2, -1, 0, 1, 2]
    err = te.getTruncationError(stencil, 1)
```

`f^(5)*h**4/30`

#### 3-point central finite difference on the staggered grid for 1st derivative
```Python
    stencil = [-0.5, 0, 0.5]
    err = te.getTruncationError(stencil, 1)
```

`-f^(3)*h**2/24`

#### 5-point central finite difference on the staggered grid for 1st derivative
```Python
    stencil = [-1.5, -0.5, 0, 0.5, 1.5]
    err = te.getTruncationError(stencil, 1)
```

`3*f^(5)*h**4/640`
