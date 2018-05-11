from constants import Borders
from isotropy import get_solution
from isotropy.method import find_dots
from mpmath import MPContext

M = 15
h = 5
a = 10

mp = MPContext()
mp.prec = 1500

# Алюминий
E = 70
nu = 0.34

for name in ['E', 'nu', 'h', 'a']:
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
solution = get_solution(M, h, a, E, nu, conditions, mp, report=True)

import tkinter

from gui.graphs import create_graph_frame

in_dots = [y for x, y in find_dots(h, a, M)[Borders.UPPER]['sigma_x']]
dots = [
    in_dots*2,
    (
        [conditions[Borders.BOTTOM]['sigma_x'](y) for y in in_dots]
        + [conditions[Borders.BOTTOM]['sigma_x'](y) for y in in_dots]
    )
]

create_graph_frame(
    tkinter.Tk(),
    func_name='sigma_x',
    var_name='y',
    func=solution[2],
    h=h,
    a=a
).pack(expand=True, fill='both')
tkinter.mainloop()
