from functools import lru_cache, partial

import time

from constants import Borders
from math_utilites import matrix, solve, pi
from .coeff_functions import get_coeff_function


def find_H(E_x, E_y, G_xy, nu_xy):
    nu_yx = nu_xy * E_y / E_x
    H = {(i, j): 0 for i in (1, 2, 6) for j in (1, 2, 6)}
    _d = 1 - nu_xy * nu_yx
    H[1, 1] = E_x / _d
    H[2, 2] = E_y / _d
    H[1, 2] = E_y * nu_xy / _d
    H[6, 6] = G_xy
    return H


def find_a_and_d(H):
    d1_0 = H[2, 2] * H[6, 6]
    d1_2 = H[1, 1] * H[2, 2] - 2 * H[1, 2] * H[6, 6] - H[1, 2]**2
    d1_4 = H[1, 1] * H[6, 6]

    d2_0 = d1_4
    d2_2 = d1_2
    d2_4 = d1_0

    diskr1 = d1_2 ** 2 - 4 * d1_0 * d1_4
    diskr2 = d2_2 ** 2 - 4 * d2_0 * d2_4
    if diskr1 <= 0 or diskr2 <= 0:
        raise Exception('invalid material characteristics')

    d1 = diskr1 ** (1/2)
    d2 = diskr2 ** (1/2)

    a11 = (2 * d1_4 / (d1_2 - d1)) ** (1 / 2)
    a12 = (2 * d1_4 / (d1_2 + d1)) ** (1 / 2)
    a21 = (2 * d2_4 / (d2_2 - d2)) ** (1 / 2)
    a22 = (2 * d2_4 / (d2_2 + d2)) ** (1 / 2)

    return [[a11, a12], [a21, a22]], [d1, d2]


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
    cos_x = [h * i / (M+1) for i in range(1, M+1)]
    sin_x = [h * (2*i - 1) / (2 * (M+1)) for i in range(1, M+1)]
    cos_y = [a * i / (M+1) for i in range(1, M+1)]
    sin_y = [a * (2*i - 1) / (2 * (M+1)) for i in range(1, M+1)]

    return {
        Borders.UPPER: {
            'u': [(0, y) for y in sin_y],
            'v': [(0, y) for y in cos_y],
            'sigma_x': [(0, y) for y in sin_y],
            'sigma_y': [(0, y) for y in sin_y],
            'tau': [(0, y) for y in cos_y],
        },
        Borders.BOTTOM: {
            'u': [(h, y) for y in sin_y],
            'v': [(h, y) for y in cos_y],
            'sigma_x': [(h, y) for y in sin_y],
            'sigma_y': [(h, y) for y in sin_y],
            'tau': [(h, y) for y in cos_y],
        },
        Borders.LEFT: {
            'u': [(x, 0) for x in cos_x],
            'v': [(x, 0) for x in sin_x],
            'sigma_x': [(x, 0) for x in sin_x],
            'sigma_y': [(x, 0) for x in sin_x],
            'tau': [(x, 0) for x in cos_x],
        },
        Borders.RIGHT: {
            'u': [(x, a) for x in cos_x],
            'v': [(x, a) for x in sin_x],
            'sigma_x': [(x, a) for x in sin_x],
            'sigma_y': [(x, a) for x in sin_x],
            'tau': [(x, a) for x in cos_x],
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
def build_system(conditions, get_coeff_func, h, a, M):
    A, b = list(), list()

    for func_name, two_eq_borders in [('u', (Borders.UPPER, Borders.BOTTOM)), ('v', (Borders.LEFT, Borders.RIGHT))]:
        coeff_func = get_coeff_func(func_name)
        this_funcs = {border: funcs[func_name] for border, funcs in conditions.items() if func_name in funcs}
        for border in two_eq_borders:
            if border in this_funcs:
                for dot in get_corner_dots(h, a, border):
                    A.append(coeff_func(*dot))
                    b.append(this_funcs[border](dot[0] if border.is_vertical() else dot[1]))
                break
        else:
            for border in this_funcs:
                dot = ({(0, 0), (h, a)} & get_corner_dots(h, a, border)).pop()
                A.append(coeff_func(*dot))
                b.append(this_funcs[border](dot[0] if border.is_vertical() else dot[1]))

    if len(A) == 3:
        border, tau_func = ({border: funcs['tau'] for border, funcs in conditions.items() if 'tau' in funcs}).popitem()
        dot = get_corner_dots(h, a, border).pop()
        A.append(get_coeff_func('tau')(*dot))
        b.append(tau_func(dot[0] if border.is_vertical() else dot[1]))

    # for row in A:
    #     print(' '*15, *[f'{float(i):=14e}' for i in row])

    dots = find_dots(h, a, M)
    coeff_func = {name: get_coeff_func(name) for name in ('u', 'v', 'sigma_x', 'sigma_y', 'tau')}
    for border in Borders:
        for name, inital_func in conditions[border].items():
            for dot in dots[border][name]:
                A.append(coeff_func[name](*dot))
                b.append(inital_func(dot[0] if border.is_vertical() else dot[1]))

                print(f'{border.name:6} {name:8}:', end='')
                print(*[f'{float(i):=14e}' for i in A[-1]], '|', b[-1])

    return matrix(A), matrix([[i, ] for i in b])


def get_solution(M, h, a, E_x, E_y, G_xy, nu_xy, conditions, *, report=False):
    if report:
        timer = Timer()

        print("Building system...")
        timer.start()

    H = find_H(E_x, E_y, G_xy, nu_xy)
    a_, d = find_a_and_d(H)

    alpha = [i * pi / h for i in range(1, M + 1)]
    beta = [i * pi / a for i in range(1, M + 1)]

    get_coeff_func = partial(get_coeff_function, H, a_, d, alpha, beta)
    get_coeff_func = lru_cache(maxsize=None)(get_coeff_func)

    A, b = build_system(conditions, get_coeff_func, h, a, M)

    if report:
        print(f"(done in {timer})")

        print("Solving system...")
        timer.start()

    coefficients = solve(A, b)

    if report:
        print(f"(done in {timer})")

    a_, d = find_a_and_d(H)

    return [
        lambda x, y, coeff_func=coeff_func: sum(i*j for i, j in zip(coeff_func(x, y), coefficients))
        for coeff_func in [
            get_coeff_func(func_name)
            for func_name in ['u', 'v', 'sigma_x', 'sigma_y', 'tau']
        ]
    ]
