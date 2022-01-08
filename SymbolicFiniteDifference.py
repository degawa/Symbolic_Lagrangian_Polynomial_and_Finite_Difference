import sympy as sp
import lagrangianpoly as lp

_DefaultIndependentVariableSymbol = 'x'
_DefaultIntervalSymbol = 'h'
_DefaultFunctionSymbol = 'f'


def createXSetFromStencil(stencil, intervalSymbol=_DefaultIntervalSymbol):
    return [stencil[i]*sp.symbols(intervalSymbol) for i in range(len(stencil))]


def createSetOfFunctionSymbolsAtXSet(xSet, functionSymbol=_DefaultFunctionSymbol):
    fSet = sp.symbols(functionSymbol+'_0:{:d}'.format(len(xSet)))
    return fSet


def getFiniteDifferenceEquation(stencil, orderOfDifference=1,
                                intervalSymbol=_DefaultIntervalSymbol,
                                x=sp.symbols(_DefaultIndependentVariableSymbol)):
    xSet = createXSetFromStencil(stencil, intervalSymbol)
    fSet = createSetOfFunctionSymbolsAtXSet(xSet, _DefaultFunctionSymbol)

    return lp.Derivative(lp.LagrangianPoly(x, xSet, fSet), x, orderOfDifference)


def getFiniteDifferenceCoefficients(stencil, orderOfDifference=1,
                                    intervalSymbol=_DefaultIntervalSymbol,
                                    x=sp.symbols(_DefaultIndependentVariableSymbol)):
    xSet = createXSetFromStencil(stencil, intervalSymbol)

    fSet = createSetOfFunctionSymbolsAtXSet(xSet, _DefaultFunctionSymbol)
    num, den = lp.LagrangianPoly(x, xSet, fSet).as_numer_denom()
    num_coef = num.as_poly(fSet).coeffs()
    den_coef = den.as_poly(fSet).coeffs()

    coef = [sp.diff(num/den_coef[0], x, orderOfDifference).subs(x, 0)
            for num in num_coef]

    return [c if c.is_number else sp.poly(c).coeffs()[0] for c in coef]
