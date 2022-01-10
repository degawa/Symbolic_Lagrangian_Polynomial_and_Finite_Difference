import sympy as sp
import utils as util


def getFiniteDifferenceEquation(stencil, orderOfDifference=1,
                                intervalSymbolStr=util._DefaultIntervalSymbolStr,
                                sameSubscriptsAsStencil=False):
    xSet = util.createXSetFromStencil(stencil, intervalSymbolStr)
    fSet = util.createSetOfFunctionSymbolsAtXSet(xSet, util._DefaultFunctionSymbolStr,
                                                 sameSubscriptsAsStencil)

    coef = getFiniteDifferenceCoefficients(stencil, orderOfDifference)

    return sp.simplify(sum([coef[i]*fSet[i] for i in range(len(fSet))])
                       / sp.symbols(intervalSymbolStr)**orderOfDifference)


def getFiniteDifferenceCoefficients(stencil, orderOfDifference=1, as_numr_denom=False):
    import lagrangianpoly as lp

    xSet = util.createXSetFromStencil(stencil, util._DefaultIntervalSymbolStr)
    fSet = util.createSetOfFunctionSymbolsAtXSet(
        xSet, util._DefaultFunctionSymbolStr)

    x = sp.symbols(util._DefaultIndependentVariableSymbolStr)
    num, den = lp.LagrangianPoly(x, xSet, fSet).as_numer_denom()
    num_coef = num.as_poly(fSet).coeffs()
    den_coef = den.as_poly(fSet).coeffs()

    coef = [lp.Derivative(num/den_coef[0], x, orderOfDifference)
            for num in num_coef]

    return util.simplifyCoefficients(coef, as_numr_denom)
