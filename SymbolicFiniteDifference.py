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


def getFiniteDifferenceCoefficients(stencil, orderOfDifference=1):
    xSet = createXSetFromStencil(stencil, _DefaultIntervalSymbolStr)
    fSet = createSetOfFunctionSymbolsAtXSet(xSet, _DefaultFunctionSymbolStr)

    x = sp.symbols(_DefaultIndependentVariableSymbolStr)
    num, den = lp.LagrangianPoly(x, xSet, fSet).as_numer_denom()
    num_coef = num.as_poly(fSet).coeffs()
    den_coef = den.as_poly(fSet).coeffs()

    coef = [sp.diff(num/den_coef[0], x, orderOfDifference).subs(x, 0)
            for num in num_coef]

    return [c if c.is_number else sp.poly(c).coeffs()[0] for c in coef]
