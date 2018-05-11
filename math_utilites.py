from functools import lru_cache
from math import trunc
from mpmath import linspace, matrix, mpf, mp, sin, cos, sinh, cosh, lu_solve as solve, pi, det
from sympy import Expr, lambdify
# from LinAlg import matrix, GaussJordan as solve, det
# from numpy import linspace, matrix, sin, cos, sinh, cosh, pi
#from numpy.linalg import det, solve
# from scipy.linalg import solve


_prec = 64

# _vals = {}
#
# def _wr_(func):
#     func_name = [k for k, v in globals().items() if v is func][0]
#     _vals[func_name] = []
#     def f(*args, **kwargs):
#         res = func(*args, **kwargs)
#         _vals[func_name].append((args, kwargs, res))
#         return res
#     return f


# # lru_cache faster than memorize
# from mpmath import memoize, maxcalls
for func_name in 'sin', 'cos', 'sinh', 'cosh':
    globals()[func_name] = lru_cache(maxsize=None)(globals()[func_name])
    # globals()[func_name] = memoize(globals()[func_name])
    # globals()[func_name] = _wr_(globals()[func_name])


def get_prec():
    return _prec


def set_prec(n):
    global _prec
    _prec = n

# TODO: defining precision, and if it's low use numpy for computation, else mpmath (usilng sympy lambdify)
