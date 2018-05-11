############################################ IT WORKS!!! ###############################################################
from functools import lru_cache, partial

import time

from constants import Borders
from orthotropy.coeff_functions import get_coeff_function
import mpmath


class Timer:
    def start(self):
        self.start_time = time.process_time()
        if 'end_time' in dir(self):
            del self.end_time

    def stop(self):
        self.end_time = time.process_time()

    def time(self):
        if 'start_time' not in dir(self):
            return 0
        if 'end_time' in dir(self):
            return self.end_time - self.start_time
        return time.process_time() - self.start_time

    def __str__(self):
        total_seconds = self.time()
        days = total_seconds // (24*60*60)
        hours = (total_seconds % (24*60*60)) // (60*60)
        minutes = (total_seconds % (60*60)) // 60
        seconds = (total_seconds % 60)
        return (
            (f"{days:.0f} d " if days else "")
            + (f"{hours:.0f} h " if hours else "")
            + (f"{minutes:.0f} min " if minutes else "")
            + f"{seconds:.4f} sec"
        )


# TODO: (fix) if M % 2 then graph of sigma_x have bound near right corner
def find_dots(h, a, M):
    ravn_x = mp.linspace(0, h, M+2)[1:-1]
    ravn_y = mp.linspace(0, a, M+2)[1:-1]
    return {
        Borders.UPPER: {
            'sigma_x': [(0, y) for y in ravn_y],
            'tau': [(0, y) for y in ravn_y],
        },
        Borders.BOTTOM: {
            'sigma_x': [(h, y) for y in ravn_y],
            'tau': [(h, y) for y in ravn_y],
        },
        Borders.LEFT: {
            'u': [(x, 0) for x in mp.linspace(0, h, M+3)[1:-1]],
            'v': [(x, 0) for x in ravn_x],
        },
        Borders.RIGHT: {
            'u': [(x, a) for x in mp.linspace(0, h, M+3)[1:-1]],
            'v': [(x, a) for x in mp.linspace(0, h, M+4)[1:-1]],
        }
    }



def get_corner_dots(h, a, border=None):
    corner_dots = {
            (Borders.UPPER, Borders.LEFT): (0, 0),
            (Borders.UPPER, Borders.RIGHT): (0, a),
            (Borders.BOTTOM, Borders.LEFT): (h, 0),
            (Borders.BOTTOM, Borders.RIGHT): (h, a)
        }

    if border is None:
        return corner_dots

    return {dot for corner, dot in corner_dots.items() if border in corner}


# TODO: optimize
def build_system(conditions, get_coeff_func, h, a, M, mp):
    A, b = list(), list()

    dots = find_dots(h, a, M)

    coeff_func = {name: get_coeff_func(name) for name in ('u', 'v', 'sigma_x', 'sigma_y', 'tau')}
    for border in Borders:
        for name, inital_func in conditions[border].items():
            for dot in dots[border][name]:
                A.append(coeff_func[name](*dot))
                b.append(inital_func(dot[0] if border.is_vertical() else dot[1]))

                print(f'{border.name:6} {name:8}:', end='')
                print(*[f'{float(i):=14e}' for i in A[-1]], '|', b[-1])

    #return matrix(A), matrix([[i, ] for i in b])
    return A, b


def get_solution(M, h, a, E_x, E_y, G_xy, nu_xy, conditions, mp=None, *, report=False):
    mp = mpmath.mp if mp is None else mp
    if report:
        timer = Timer()

        print("Building system...")
        timer.start()

    alpha = [i * mp.pi / h for i in range(1, M + 1)]
    beta = [i * mp.pi / a for i in range(1, M + 1)]

    get_coeff_func = partial(get_coeff_function, E_x, E_y, G_xy, nu_xy, alpha, beta, mp)
    get_coeff_func = lru_cache(maxsize=None)(get_coeff_func)

    A, b = [mp.matrix(i) for i in build_system(conditions, get_coeff_func, h, a, M, mp)]

    if report:
        print(f"(done in {timer})")

        print("Solving system...")
        timer.start()

    coefficients = mp.lu_solve(A, b)

    if report:
        print(f"(done in {timer})")

    # TODO: solution -> dict
    # TODO: scalar product
    return [
        lambda x, y, coeff_func=coeff_func: sum(i*j for i, j in zip(coeff_func(x, y), coefficients))
        for coeff_func in [
            get_coeff_func(func_name)
            for func_name in ['u', 'v', 'sigma_x', 'sigma_y', 'tau']
        ]
    ]


########################################################################################################################
import matplotlib.pyplot as plt
from constants import Borders
from mpmath import MPContext
from functools import lru_cache

M = 14
mp = MPContext()
mp.prec = 1500
h = 5
a = 10

# Углеводородная пластина
E_x = 6.9 * 10**9
E_y = 220 * 10**9
G_xy = 5 * 10**9
nu_xy = 0.008

# # Березовая фанера
# E_x = 1.2 * 10**5
# E_y = 0.6 * 10**5
# G_xy = 0.07 * 10**5
# nu_xy = 0.036

for name in ['E_x', 'E_y', 'G_xy', 'nu_xy', 'h', 'a']:
    globals()[name] = mp.mpf(globals()[name])

conditions = {
    Borders.UPPER: {
        'sigma_x': lambda y: 1,
        'tau': lambda *_: 0
    },
    Borders.BOTTOM: {
        'sigma_x': lambda *_: 0,
        'tau': lambda *_: 0
    },
    Borders.LEFT: {
        'u': lambda *_: 0,
        'v': lambda *_: 0
    },
    Borders.RIGHT: {
        'u': lambda *_: 0,
        'v': lambda *_: 0
    }
}

##########################################
mp.sin = lru_cache(maxsize=None)(mp.sin)
mp.cos = lru_cache(maxsize=None)(mp.cos)
mp.sinh = lru_cache(maxsize=None)(mp.sinh)
mp.cosh = lru_cache(maxsize=None)(mp.cosh)

solution = get_solution(M, h, a, E_x, E_y, G_xy, nu_xy, conditions, mp, report=True)

# import tkinter
#
# from gui.graphs import create_graph_frame
#
# in_dots = [y for x, y in find_dots(h, a, M)[Borders.UPPER]['sigma_x']]
# dots = [
#     in_dots*2,
#     (
#         [conditions[Borders.BOTTOM]['sigma_x'](y) for y in in_dots]
#         + [conditions[Borders.BOTTOM]['sigma_x'](y) for y in in_dots]
#     )
# ]
#
# create_graph_frame(
#     tkinter.Tk(),
#     func_name='sigma_x',
#     var_name='y',
#     func=solution[2],
#     h=h,
#     a=a
# ).pack(expand=True, fill='both')
# tkinter.mainloop()


sigma_x = solution[2]

timer = Timer()
print("Plotting graph...")
timer.start()

y_dots = mp.linspace(0, a, 200)

list_of_vals = []
for x in [i*h/4 for i in range(5)]:
    list_of_vals.append(
        [sigma_x(x, y) for y in y_dots]
    )

for vals in list_of_vals:
    plt.plot(y_dots, vals)

dots = find_dots(h, a, M)[Borders.UPPER]['sigma_x']
plt.plot(
    [y for x, y in dots],
    [conditions[Borders.BOTTOM]['sigma_x'](y) for x, y in dots],
    '*'
)
plt.plot(
    [y for x, y in dots],
    [conditions[Borders.UPPER]['sigma_x'](y) for x, y in dots],
    '*'
)

print(f'(done in {timer})')
plt.show()