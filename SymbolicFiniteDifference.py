import sympy as sp
import lagrangianpoly as lp

_DefaultIndependentVariableSymbolStr = 'x'
_DefaultIntervalSymbolStr = 'h'
_DefaultFunctionSymbolStr = 'f'


def createXSetFromStencil(stencil, intervalSymbolStr=_DefaultIntervalSymbolStr):
    return [stencil[i]*sp.symbols(intervalSymbolStr) for i in range(len(stencil))]


def createSetOfFunctionSymbolsAtXSet(xSet, functionSymbolStr=_DefaultFunctionSymbolStr):
    fSet = sp.symbols(functionSymbolStr+'_0:{:d}'.format(len(xSet)))
    return fSet


def getFiniteDifferenceEquation(stencil, orderOfDifference=1,
                                intervalSymbolStr=_DefaultIntervalSymbolStr):
    xSet = createXSetFromStencil(stencil, intervalSymbolStr)
    fSet = createSetOfFunctionSymbolsAtXSet(xSet, _DefaultFunctionSymbolStr)

    x = sp.symbols(_DefaultIndependentVariableSymbolStr)
    return lp.Derivative(lp.LagrangianPoly(x, xSet, fSet), x, orderOfDifference)


def getFiniteDifferenceCoefficients(stencil, orderOfDifference=1, as_numr_denom=False):

    xSet = createXSetFromStencil(stencil, _DefaultIntervalSymbolStr)
    fSet = createSetOfFunctionSymbolsAtXSet(xSet, _DefaultFunctionSymbolStr)

    x = sp.symbols(_DefaultIndependentVariableSymbolStr)
    num, den = lp.LagrangianPoly(x, xSet, fSet).as_numer_denom()
    num_coef = num.as_poly(fSet).coeffs()
    den_coef = den.as_poly(fSet).coeffs()

    coef = [sp.diff(num/den_coef[0], x, orderOfDifference).subs(x, 0)
            for num in num_coef]

    coef_num = [c if c.is_number else sp.poly(c).coeffs()[0] for c in coef]

    import numpy as np
    import fractions as fr
    coef_rational = [sp.Rational(fr.Fraction(
        str(c)).limit_denominator(100000)) for c in coef_num]

    denom = [c.q for c in coef_rational]
    denom_lcm = np.lcm.reduce(np.array(denom))
    numr = [c*denom_lcm for c in coef_rational]

    if as_numr_denom:
        return numr, denom_lcm
    else:
        return [n/denom_lcm for n in numr]
