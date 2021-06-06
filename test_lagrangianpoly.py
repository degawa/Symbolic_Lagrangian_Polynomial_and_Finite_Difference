import unittest
from lagrangianpoly import LagrangianBasis, LagrangianPoly, Derivative


class TestLagrangianPoly(unittest.TestCase):

    def test_LagrangianBasis_5stencil(self):
        import sympy as sp

        x = sp.symbols('x')
        acctual = LagrangianBasis(x, degreeOfPolynomial=5, pointAt=0)

        x0 = sp.symbols('x0')
        x1 = sp.symbols('x1')
        x2 = sp.symbols('x2')
        x3 = sp.symbols('x3')
        x4 = sp.symbols('x4')
        x5 = sp.symbols('x5')
        expected = (x - x1)*(x - x2)*(x - x3)*(x - x4)*(x - x5) \
            / ((x0 - x1)*(x0 - x2)*(x0 - x3)*(x0 - x4)*(x0 - x5))

        self.assertEqual(expected, acctual)

    def test_LagrangianPoly_3stencil(self):
        import sympy as sp

        x = sp.symbols('x')
        dx = sp.symbols('dx')

        x_set = [-dx, 0, dx]
        f_set = sp.symbols('f0:{:d}'.format(len(x_set)))
        acctual = LagrangianPoly(x, x_set, f_set)

        f0 = sp.symbols('f0')
        f1 = sp.symbols('f1')
        f2 = sp.symbols('f2')

        expected = f0*x*(-dx + x)/(2*dx**2) \
            - f1*(-dx + x) * (dx + x)/dx**2 \
            + f2*x*(dx + x)/(2*dx**2)

        self.assertEqual(expected, acctual)

    def test_Derivative_3stencil_regular_central_order1(self):
        import sympy as sp

        x = sp.symbols('x')
        dx = sp.symbols('dx')

        f0 = sp.symbols('f0')
        f1 = sp.symbols('f1')
        f2 = sp.symbols('f2')
        Lagrangian_poly = f0*x*(-dx + x)/(2*dx**2) \
            - f1*(-dx + x) * (dx + x)/dx**2 \
            + f2*x*(dx + x)/(2*dx**2)
        acctual = Derivative(Lagrangian_poly, x, orderOfDifference=1)

        f0 = sp.symbols('f0')
        f1 = sp.symbols('f1')
        f2 = sp.symbols('f2')

        expected = (-f0 + f2)/(2*dx)

        self.assertEqual(expected, acctual)

    def test_Derivative_3stencil_regular_rightSidede_order1(self):
        import sympy as sp

        x = sp.symbols('x')
        dx = sp.symbols('dx')

        f0 = sp.symbols('f0')
        f1 = sp.symbols('f1')
        f2 = sp.symbols('f2')
        Lagrangian_poly = f0*(-2*dx + x)*(-dx + x)/(2*dx**2) \
            - f1*x*(-2*dx + x)/dx**2 \
            + f2*x*(-dx + x)/(2*dx**2)
        acctual = Derivative(Lagrangian_poly, x, orderOfDifference=1)

        f0 = sp.symbols('f0')
        f1 = sp.symbols('f1')
        f2 = sp.symbols('f2')

        expected = (-3*f0 + 4*f1 - f2)/(2*dx)

        self.assertEqual(expected, acctual)

    def test_Derivative_3stencil_staggered_central_order1(self):
        import sympy as sp

        x = sp.symbols('x')
        dx = sp.symbols('dx')

        f0 = sp.symbols('f0')
        f1 = sp.symbols('f1')
        f2 = sp.symbols('f2')
        Lagrangian_poly = 2*f0*x*(-dx/2 + x)/dx**2 \
            - 4*f1*(-dx/2 + x)*(dx/2 + x)/dx**2 \
            + 2*f2*x*(dx/2 + x)/dx**2
        acctual = Derivative(Lagrangian_poly, x, orderOfDifference=1)

        f0 = sp.symbols('f0')
        f2 = sp.symbols('f2')

        expected = (-f0 + f2)/dx

        self.assertEqual(expected, acctual)


if __name__ == "__main__":
    unittest.main()
