import matplotlib.pyplot as plt
from math_utilites import mp, mpf, linspace
from constants import Borders
from orthotropy import get_solution
from orthotropy.method import Timer, find_dots
from mpmath import MPContext
from functools import lru_cache

M = 15
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

y_dots = linspace(0, a, 200)

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