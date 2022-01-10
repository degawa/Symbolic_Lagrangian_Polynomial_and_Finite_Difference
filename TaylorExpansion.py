import sympy as sp
import utils as util
import SymbolicFiniteDifference as fd


def TaylorExpansion(h, n):
    df_set = sp.symbols(util._DefaultFunctionSymbolStr +
                        '^((1:{:d}))'.format(n+1))

    coef = [h**i*sp.Rational(1, sp.factorial(i)) for i in range(1, n+1)]
    f = sp.symbols(util._DefaultFunctionSymbolStr)
    te = f
    for i in range(len(df_set)):
        te += df_set[i]*coef[i]

    return te


def _getDerivativeSymbol(functionSymbolStr, n):
    return sp.symbols(functionSymbolStr+"^(%d)" % n)


def getTruncationError(stencil, orderOfDifference,
                       intervalSymbolStr=util._DefaultIntervalSymbolStr):
    xSet = util.createXSetFromStencil(
        stencil, intervalSymbolStr=intervalSymbolStr)

    coef = fd.getFiniteDifferenceCoefficients(xSet, orderOfDifference)

    num_expterm = len(xSet)+orderOfDifference
    f_te = [TaylorExpansion(x, num_expterm) for x in xSet]

    eq = sum([coef[i]*f_te[i] for i in range(len(xSet))])

    intervalSymbol = sp.symbols(intervalSymbolStr)
    return sp.simplify(_getDerivativeSymbol(util._DefaultFunctionSymbolStr, orderOfDifference)
                       - sp.nsimplify(eq/intervalSymbol**orderOfDifference,
                                      rational=True, tolerance=1e-10))
