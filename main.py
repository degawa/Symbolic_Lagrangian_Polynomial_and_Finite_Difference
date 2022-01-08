from lagrangianpoly import LagrangianBasis, LagrangianPoly, Derivative
import SymbolicFiniteDifference as fd


def main():
    import sympy as sp

    # lagrangian basis
    x = sp.symbols('x')
    eq = LagrangianBasis(x, degreeOfPolynomial=5, pointAt=0)
    print(eq)
    # (x - x1)*(x - x2)*(x - x3)*(x - x4)*(x - x5)/((x0 - x1)*(x0 - x2)*(x0 - x3)*(x0 - x4)*(x0 - x5))

    x = sp.symbols('x')
    dx = sp.symbols('dx')

    # 3-point stencil lagrangian polynomial
    x_set = [-dx, 0, dx]
    f_set = sp.symbols('f0:{:d}'.format(len(x_set)))
    eq = LagrangianPoly(x, x_set, f_set)
    print(eq)
    # f0*x*(-dx + x)/(2*dx**2) - f1*(-dx + x)*(dx + x)/dx**2 + f2*x*(dx + x)/(2*dx**2)

    # 3-point stencil 1st order central difference on regular grid
    x_set = [-dx, 0, dx]
    f_set = sp.symbols('f0:{:d}'.format(len(x_set)))
    eq = Derivative(LagrangianPoly(x, x_set, f_set), x, orderOfDifference=1)
    print(eq)

    # 3-point stencil 2nd order central difference on regular grid
    x_set = [-dx, 0, dx]
    f_set = sp.symbols('f0:{:d}'.format(len(x_set)))
    eq = Derivative(LagrangianPoly(x, x_set, f_set), x, orderOfDifference=2)
    print(eq)

    # 5-point stencil 1st order central difference on regular grid
    x_set = [-2*dx, -dx, 0, dx, 2*dx]
    f_set = sp.symbols('f0:{:d}'.format(len(x_set)))
    eq = Derivative(LagrangianPoly(x, x_set, f_set), x, orderOfDifference=1)
    print(eq)

    # 3-point stencil 1st order right-sided difference on regular grid
    x_set = [0, dx, 2*dx]
    f_set = sp.symbols('f0:{:d}'.format(len(x_set)))
    eq = Derivative(LagrangianPoly(x, x_set, f_set), x, orderOfDifference=1)
    print(eq)

    # 4-point stencil 2nd order right-sided difference on regular grid
    x_set = [0, dx, 2*dx, 3*dx]
    f_set = sp.symbols('f0:{:d}'.format(len(x_set)))
    eq = Derivative(LagrangianPoly(x, x_set, f_set), x, orderOfDifference=2)
    print(eq)

    # 3-point stencil 1st order central difference on staggered grid
    x_set = [-dx/2, 0, dx/2]
    f_set = sp.symbols('f0:{:d}'.format(len(x_set)))
    eq = Derivative(LagrangianPoly(x, x_set, f_set), x, orderOfDifference=1)
    print(eq)

    # 5-point 1st order central difference with a default interval symbol
    stencil = [-2, -1, 0, 1, 2]
    eq = fd.getFiniteDifferenceEquation(
        stencil, orderOfDifference=1)
    print(eq)

    # 5-point 1st order central difference with an user-defined interval symbol
    stencil = [-2, -1, 0, 1, 2]
    eq = fd.getFiniteDifferenceEquation(
        stencil, orderOfDifference=1, intervalSymbol='dx')
    print(eq)

    # 3-point one-sided 3-point 1st order central difference with a default interval symbol
    stencil = [0, 1, 2]
    eq = fd.getFiniteDifferenceEquation(
        stencil, orderOfDifference=1)
    print(eq)

    # 4-point one-sided 4-point 2nd order central difference with a default interval symbol
    stencil = [0, 1, 2, 3]
    eq = fd.getFiniteDifferenceEquation(
        stencil, orderOfDifference=2)
    print(eq)

    # coefficients for 5-point 1st order central difference
    stencil = [-2, -1, 0, 1, 2]
    coef = fd.getFiniteDifferenceCoefficients(
        stencil, orderOfDifference=1)
    print(coef)

    # coefficients for 7-point 1st order central difference
    stencil = [-3, -2, -1, 0, 1, 2, 3]
    coef = fd.getFiniteDifferenceCoefficients(
        stencil, orderOfDifference=1)
    print(coef)

    # coefficients for one-sided 3-point 1st order central difference
    stencil = [0, 1, 2]
    coef = fd.getFiniteDifferenceCoefficients(
        stencil, orderOfDifference=1)
    print(coef)

    # coefficients for one-sided 4-point 2nd order central difference
    stencil = [0, 1, 2, 3]
    coef = fd.getFiniteDifferenceCoefficients(
        stencil, orderOfDifference=2)
    print(coef)


if __name__ == "__main__":
    main()
