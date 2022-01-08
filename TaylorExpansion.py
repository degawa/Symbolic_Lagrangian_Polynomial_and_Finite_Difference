import sympy as sp
import SymbolicFiniteDifference as fd


def TaylorExpansion(h, n):
    df_set = sp.symbols(fd._DefaultFunctionSymbolStr+'^((1:{:d}))'.format(n+1))

    coef = [h**i*sp.Rational(1, sp.factorial(i)) for i in range(1, n+1)]
    f = sp.symbols(fd._DefaultFunctionSymbolStr)
    te = f
    for i in range(len(df_set)):
        te += df_set[i]*coef[i]

    return te


def _getDerivativeSymbol(functionSymbolStr, n):
    return sp.symbols(functionSymbolStr+"^(%d)" % n)


def getTruncationError(stencil, orderOfDifference,
                       intervalSymbolStr=fd._DefaultIntervalSymbolStr):
    xSet = fd.createXSetFromStencil(
        stencil, intervalSymbolStr=intervalSymbolStr)

    coef = fd.getFiniteDifferenceCoefficients(xSet, orderOfDifference)

    num_expterm = len(xSet)+orderOfDifference
    f_te = [TaylorExpansion(x, num_expterm) for x in xSet]

    eq = sum([coef[i]*f_te[i] for i in range(len(xSet))])

    intervalSymbol = sp.symbols(intervalSymbolStr)
    return _getDerivativeSymbol(fd._DefaultFunctionSymbolStr, orderOfDifference)\
        - sp.simplify(eq/intervalSymbol**orderOfDifference)
