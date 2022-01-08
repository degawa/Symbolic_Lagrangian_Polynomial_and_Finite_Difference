# Symbolic Constructor for Lagrangian Polynomial and Finite Difference
a python script to symbolically construct Lagrange polynomials and finite difference equations

# Usage
donwload lagrangianpoly.py and SymbolicFiniteDifference.py and import it like below:

```Python
from lagrangianpoly import LagrangianBasis, LagrangianPoly, Derivative
import SymbolicFiniteDifference as fd
```

# Todo
- [ ] add comments
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
#### 3-point stencil 1st order central difference on the regular grid

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

#### 3-point stencil 1st order central difference on the staggered grid

Given a set of data \((x_0, f_0), (x_1, f_1), (x_2, f_2)\) where \((x_0, x_1, x_2) = (\frac{-\varDelta x}{2}, 0, \frac{\varDelta x}{2})\), 3-point stencil 1st order central difference on the staggered grid at \(x = 0\), e.g. \(x = x_1\), is constructed by the following script:

```Python
    x_set = [-dx/2, 0, dx/2]
    f_set = sp.symbols('f0:{:d}'.format(len(x_set)))
    eq = Derivative(LagrangianPoly(x, x_set, f_set), x, orderOfDifference=1)
```

We get `(-f0 + f2)/dx`.

## Symbolic Finite Difference
### getFiniteDifferenceEquation

#### 5-point 1st order central difference with a default interval symbol
```Python
    stencil = [-2, -1, 0, 1, 2]
    eq = fd.getFiniteDifferenceEquation(stencil, orderOfDifference=1)
```

The result is `(f_0 - 8*f_1 + 8*f_3 - f_4)/(12*h)`.

#### 5-point 1st order central difference with an user-defined interval symbol
```Python
    stencil = [-2, -1, 0, 1, 2]
    eq = fd.getFiniteDifferenceEquation(stencil, orderOfDifference=1, intervalSymbolStr='dx')
```

`(f_0 - 8*f_1 + 8*f_3 - f_4)/(12*dx)`

#### 3-point 1st order one-sided difference with a default interval symbol
```Python
    stencil = [0, 1, 2]
    eq = fd.getFiniteDifferenceEquation(stencil, orderOfDifference=1)
```

`(-3*f_0 + 4*f_1 - f_2)/(2*h)`

#### 4-point 2nd order one-sided difference with a default interval symbol
```Python
    stencil = [0, 1, 2, 3]
    eq = fd.getFiniteDifferenceEquation(stencil, orderOfDifference=2)
```

`(2*f_0 - 5*f_1 + 4*f_2 - f_3)/h**2`

### getFiniteDifferenceCoefficients
#### 5-point 1st order central difference
```Python
    stencil = [-2, -1, 0, 1, 2]
    coef = fd.getFiniteDifferenceCoefficients(stencil, orderOfDifference=1)
```

The result is `[1/12, -2/3, 0, 2/3, -1/12]`.

#### 7-point 1st order central difference
```Python
    stencil = [-3, -2, -1, 0, 1, 2, 3]
    coef = fd.getFiniteDifferenceCoefficients(stencil, orderOfDifference=1)
```

`[-1/60, 3/20, -3/4, 0, 3/4, -3/20, 1/60]`

#### 3-point 1st order one-sided difference
```Python
    stencil = [0, 1, 2]
    coef = fd.getFiniteDifferenceCoefficients(stencil, orderOfDifference=1)
```

`[-3/2, 2, -1/2]`

#### 4-point 2nd order one-sided difference
```Python
    stencil = [0, 1, 2, 3]
    coef = fd.getFiniteDifferenceCoefficients(stencil, orderOfDifference=2)
```

`[2, -5, 4, -1]`
