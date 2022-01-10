from lagrangianpoly import LagrangianBasis, LagrangianPoly, Derivative
import SymbolicFiniteDifference as fd
import SymbolicInterpolation as intp
import TaylorExpansion as te


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

    # 3-point stencil central difference for 1st derivative on regular grid
    x_set = [-dx, 0, dx]
    f_set = sp.symbols('f0:{:d}'.format(len(x_set)))
    eq = Derivative(LagrangianPoly(x, x_set, f_set), x, orderOfDifference=1)
    print(eq)

    # 3-point stencil central difference for 2nd derivative on regular grid
    x_set = [-dx, 0, dx]
    f_set = sp.symbols('f0:{:d}'.format(len(x_set)))
    eq = Derivative(LagrangianPoly(x, x_set, f_set), x, orderOfDifference=2)
    print(eq)

    # 5-point stencil central difference for 1st derivative on regular grid
    x_set = [-2*dx, -dx, 0, dx, 2*dx]
    f_set = sp.symbols('f0:{:d}'.format(len(x_set)))
    eq = Derivative(LagrangianPoly(x, x_set, f_set), x, orderOfDifference=1)
    print(eq)

    # 3-point stencil right-sided difference for 1st derivative on regular grid
    x_set = [0, dx, 2*dx]
    f_set = sp.symbols('f0:{:d}'.format(len(x_set)))
    eq = Derivative(LagrangianPoly(x, x_set, f_set), x, orderOfDifference=1)
    print(eq)

    # 4-point stencil right-sided difference for 2nd derivative on regular grid
    x_set = [0, dx, 2*dx, 3*dx]
    f_set = sp.symbols('f0:{:d}'.format(len(x_set)))
    eq = Derivative(LagrangianPoly(x, x_set, f_set), x, orderOfDifference=2)
    print(eq)

    # 3-point stencil central difference for 1st derivative on staggered grid
    x_set = [-dx/2, 0, dx/2]
    f_set = sp.symbols('f0:{:d}'.format(len(x_set)))
    eq = Derivative(LagrangianPoly(x, x_set, f_set), x, orderOfDifference=1)
    print(eq)

    # 5-point central difference for 1st derivative with a default interval symbol
    stencil = [-2, -1, 0, 1, 2]
    eq = fd.getFiniteDifferenceEquation(
        stencil, orderOfDifference=1)
    print(eq)

    # 5-point central difference for 1sit derivative with an user-defined interval symbol
    stencil = [-2, -1, 0, 1, 2]
    eq = fd.getFiniteDifferenceEquation(
        stencil, orderOfDifference=1, intervalSymbolStr='dx')
    print(eq)

    # 3-point one-sided difference for 1st derivative with a default interval symbol
    stencil = [0, 1, 2]
    eq = fd.getFiniteDifferenceEquation(
        stencil, orderOfDifference=1)
    print(eq)

    # 4-point one-sided difference for 2nd derivative with a default interval symbol
    stencil = [0, 1, 2, 3]
    eq = fd.getFiniteDifferenceEquation(
        stencil, orderOfDifference=2)
    print(eq)

    # 5-point central difference formula for 1st derivative with changing subscripts
    stencil = [-2, -1, 0, 1, 2]
    eq = fd.getFiniteDifferenceEquation(
        stencil, orderOfDifference=1, sameSubscriptsAsStencil=True)
    print(eq)

    # 5-point central difference formula for 1st derivative with changing subscripts
    stencil = [-1.5, -0.5, 0, 0.5, 1.5]
    eq = fd.getFiniteDifferenceEquation(
        stencil, orderOfDifference=1, sameSubscriptsAsStencil=True)
    print(eq)

    # coefficients for 5-point central difference for 1st derivative
    stencil = [-2, -1, 0, 1, 2]
    coef = fd.getFiniteDifferenceCoefficients(
        stencil, orderOfDifference=1)
    print(coef)

    # coefficients for 7-point central difference for 1st derivative
    stencil = [-3, -2, -1, 0, 1, 2, 3]
    coef = fd.getFiniteDifferenceCoefficients(
        stencil, orderOfDifference=1)
    print(coef)

    # coefficients for 3-point one-sided difference for 1st derivative
    stencil = [0, 1, 2]
    coef = fd.getFiniteDifferenceCoefficients(
        stencil, orderOfDifference=1)
    print(coef)

    # coefficients for 4-point one-sided difference for 2nd derivative
    stencil = [0, 1, 2, 3]
    coef = fd.getFiniteDifferenceCoefficients(
        stencil, orderOfDifference=2)
    print(coef)

    # numerator and denominator of coefficients for 5-point central difference for 1st derivative
    stencil = [-2, -1, 0, 1, 2]
    numr, denom = fd.getFiniteDifferenceCoefficients(
        stencil, orderOfDifference=1, as_numr_denom=True)
    print(numr, denom)

    # numerator and denominator of coefficients for 7-point central difference for 1st derivative
    stencil = [-3, -2, -1, 0, 1, 2, 3]
    numr, denom = fd.getFiniteDifferenceCoefficients(
        stencil, orderOfDifference=1, as_numr_denom=True)
    print(numr, denom)

    # numerator and denominator of coefficients for 3-point one-sided difference for 1st derivative
    stencil = [0, 1, 2]
    numr, denom = fd.getFiniteDifferenceCoefficients(
        stencil, orderOfDifference=1, as_numr_denom=True)
    print(numr, denom)

    # numerator and denominator of coefficients for 3-point one-sided difference for 1st derivative on staggered grid
    stencil = [-0.5, 0, 0.5]
    numr, denom = fd.getFiniteDifferenceCoefficients(
        stencil, orderOfDifference=1, as_numr_denom=True)
    print(numr, denom)

    # Tayloer expansion of f(x+h) around x up to term including 6th order difference
    h = sp.symbols('h')
    f1 = te.TaylorExpansion(h, n=6)
    print(f1)

    # Tayloer expansion of f(x-h) around x up to term including 7th order difference
    h = sp.symbols('h')
    f_1 = te.TaylorExpansion(-h, n=7)
    print(f_1)

    # Tayloer expansion of f(x+2h) around x up to term including 4th order difference
    h = sp.symbols('h')
    f2 = te.TaylorExpansion(2*h, n=4)
    print(f2)

    # truncation error of 3-point central finite difference for 1st derivative
    stencil = [-1, 0, 1]
    err = te.getTruncationError(stencil, 1)
    print(err)

    # truncation error of 3-point central finite difference for 2nd derivative
    stencil = [-1, 0, 1]
    err = te.getTruncationError(stencil, 2)
    print(err)

    # truncation error of 5-point central finite difference for 1st derivative
    stencil = [-2, -1, 0, 1, 2]
    err = te.getTruncationError(stencil, 1)
    print(err)

    # truncation error of 3-point central finite difference for 1st dervative on the staggered grid
    stencil = [-0.5, 0, 0.5]
    err = te.getTruncationError(stencil, 1)
    print(err)

    # truncation error of 3-point central finite difference for 1st derivative on the staggered grid
    stencil = [-1.5, -0.5, 0, 0.5, 1.5]
    err = te.getTruncationError(stencil, 1)
    print(err)

    # 2-point central interpolation
    stencil = [-1, 1]
    eq = intp.getInterpolationEquation(stencil)
    print(eq)

    # 4-point central interpolation
    stencil = [-2, -1, 1, 2]
    eq = intp.getInterpolationEquation(stencil)
    print(eq)

    # coefficients for 2-point central interpolation
    stencil = [-1, 1]
    coef = intp.getInterpolationCoefficients(stencil)
    print(coef)

    # coefficients 4-point central interpolation
    stencil = [-2, -1, 1, 2]
    coef = intp.getInterpolationCoefficients(stencil)
    print(coef)

    # numerator and denominator of coefficients for 2-point central interpolation
    stencil = [-1, 1]
    numr, denom = intp.getInterpolationCoefficients(
        stencil, as_numr_denom=True)
    print(numr, denom)

    # numerator and denominator of coefficients for 4-point central interpolation
    stencil = [-2, -1, 1, 2]
    numr, denom = intp.getInterpolationCoefficients(
        stencil, as_numr_denom=True)
    print(numr, denom)

    # truncation error of 2-point central interpolation
    stencil = [-1, 1]
    err = te.getTruncationErrorOfInterpolationEquation(stencil)
    print(err)

    # truncation error of 4-point central interpolation
    stencil = [-2, -1, 1, 2]
    err = te.getTruncationErrorOfInterpolationEquation(stencil)
    print(err)

    # linear extrapolation
    stencil = [1, 2]
    eq = intp.getInterpolationEquation(stencil)
    err = te.getTruncationErrorOfInterpolationEquation(stencil)
    print(eq)
    print(err)

    # quadratic extrapolation
    stencil = [1, 2, 3]
    eq = intp.getInterpolationEquation(stencil)
    err = te.getTruncationErrorOfInterpolationEquation(stencil)
    print(eq)
    print(err)

    return


if __name__ == "__main__":
    main()
