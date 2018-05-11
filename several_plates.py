from orthotropy.coeff_functions import get_coeff_function as orthotropy_get_coeff_function
from isotropy.coeff_functions import get_coeff_function as isotropy_get_coeff_function
from orthotropy.method import Timer, find_dots, get_corner_dots
from functools import lru_cache, partial
from mpmath import MPContext
from constants import Borders, MaterialTypes
import matplotlib.pyplot as plt

#               a1          a2
#         <-----------><----------->
#          ___________  ___________ -----> y
#    ∧    |           ||           |
#    |    |   plate1  ||  plate2   |
#  h |    |           ||           |
#    ∨    |           ||           |
#        | ‾‾‾‾‾‾‾‾‾‾‾  ‾‾‾‾‾‾‾‾‾‾‾
#        |
#        |
#        ∨
#       x


M = 15
mp = MPContext()
mp.prec = 1500
h = 5
a1 = 5
a2 = 5
material_type = MaterialTypes.ORTHOTROPY


a = a1 + a2

plate1 = {'h': h, 'a': a1}
plate2 = {'h': h, 'a': a2}

if material_type == MaterialTypes.ORTHOTROPY:
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

    plate1['E_x'] = plate2['E_x'] = E_x
    plate1['E_y'] = plate2['E_y'] = E_y
    plate1['G_xy'] = plate2['G_xy'] = G_xy
    plate1['nu_xy'] = plate2['nu_xy'] = nu_xy

elif material_type == MaterialTypes.ISOTROPY:
    # Алюминий
    E = 70
    nu = 0.34

    plate1['E'] = plate2['E'] = E
    plate1['nu'] = plate2['nu'] = nu

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

h = mp.mpf(h)
a = mp.mpf(a)
a1 = mp.mpf(a1)
a2 = mp.mpf(a2)
for key in plate1:
    plate1[key] = mp.mpf(plate1[key])
    plate2[key] = mp.mpf(plate2[key])

########################################################################################################################
############################################# SYSTEM BUILDING: #########################################################
########################################################################################################################

report = True

mp.sin = lru_cache(maxsize=None)(mp.sin)
mp.cos = lru_cache(maxsize=None)(mp.cos)
mp.sinh = lru_cache(maxsize=None)(mp.sinh)
mp.cosh = lru_cache(maxsize=None)(mp.cosh)


if report:
    timer = Timer()

    print("Building system...")
    timer.start()

for plate in plate1, plate2:
    plate['alpha'] = [i * mp.pi / plate['h'] for i in range(1, M + 1)]
    plate['beta'] = [i * mp.pi / plate['a'] for i in range(1, M + 1)]

    if material_type is MaterialTypes.ORTHOTROPY:
        plate['get_coeff_func'] = partial(
            orthotropy_get_coeff_function,
            plate['E_x'],
            plate['E_y'],
            plate['G_xy'],
            plate['nu_xy'],
            plate['alpha'],
            plate['beta'],
            mp,
        )

    elif material_type is MaterialTypes.ISOTROPY:
        plate['get_coeff_func'] = partial(
            isotropy_get_coeff_function,
            plate['E'],
            plate['nu'],
            plate['alpha'],
            plate['beta'],
            mp,
        )

    plate['get_coeff_func'] = lru_cache(maxsize=None)(plate['get_coeff_func'])


A, b = list(), list()

################################# equations from conditions for plate1: ################################################

print('for plate1:')

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
            'tau': [(h, y) for y in mp.linspace(0, a, M+3)[1:-1]],
        },
        Borders.LEFT: {
            'u': [(x, 0) for x in mp.linspace(0, h, M+3)[1:-1]],
            'v': [(x, 0) for x in mp.linspace(0, h, M+4)[1:-1]],
        }
    }

a = plate1['a']

get_coeff_func = plate1['get_coeff_func']

'''
_conditions = {border: funcs for border, funcs in conditions.items() if border != Borders.RIGHT}

for func_name, two_eq_borders in [('u', (Borders.UPPER, Borders.BOTTOM)), ('v', (Borders.LEFT, Borders.RIGHT))]:
    coeff_func = get_coeff_func(func_name)
    this_funcs = {border: funcs[func_name] for border, funcs in _conditions.items() if func_name in funcs}
    for border in two_eq_borders:
        if border in this_funcs:
            for dot in get_corner_dots(h, a, border):
                A.append(coeff_func(*dot) + [0, ] * (8*M + 4))
                b.append(this_funcs[border](dot[0] if border.is_vertical() else dot[1]))
            break
    else:
        for border in this_funcs:
            dot = ({(0, 0), (h, a)} & get_corner_dots(h, a, border)).pop()
            A.append(coeff_func(*dot) + [0, ] * (8*M + 4))
            b.append(this_funcs[border](dot[0] if border.is_vertical() else dot[1]))

if len(A) == 3:
    border, tau_func = ({border: funcs['tau'] for border, funcs in _conditions.items() if 'tau' in funcs}).popitem()
    dot = get_corner_dots(h, a, border).pop()
    A.append(get_coeff_func('tau')(*dot) + [0, ] * (8*M + 4))
    b.append(tau_func(dot[0] if border == Borders.LEFT else dot[1]))
'''

'''
A.append(get_coeff_func('u')(0, 0) + [0, ] * (8*M + 4)); b.append(conditions[Borders.LEFT]['u'](0, 0))
A.append(get_coeff_func('v')(0, 0) + [0, ] * (8*M + 4)); b.append(conditions[Borders.LEFT]['v'](0, 0))
A.append(get_coeff_func('v')(h, 0) + [0, ] * (8*M + 4)); b.append(conditions[Borders.LEFT]['v'](h, 0))
A.append(get_coeff_func('tau')(0, 0) + [0, ] * (8*M + 4)); b.append(conditions[Borders.UPPER]['tau'](0, 0))
'''

# for row in A:
#     print(' '*15, *[f'{float(i):=14e}' for i in row])

dots = find_dots(h, a, M)
coeff_func = {name: get_coeff_func(name) for name in ('u', 'v', 'sigma_x', 'sigma_y', 'tau')}
for border in filter(lambda b: b != Borders.RIGHT, Borders):
    for name, inital_func in conditions[border].items():
        for dot in dots[border][name]:
            A.append(coeff_func[name](*dot) + [0, ] * (8*M + 4))
            b.append(inital_func(dot[0] if border.is_vertical() else dot[1]))

            # print(f'{border.name:6} {name:8}:', end='')
            # print(*[f'{float(i):=14e}' for i in A[-1]], '|', b[-1])

################################# equations from conditions for plate2: ################################################

print('for plate2:')

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
            'tau': [(h, y) for y in mp.linspace(0, a, M+3)[1:-1]],
        },
        Borders.RIGHT: {
            'u': [(x, a) for x in mp.linspace(0, h, M+3)[1:-1]],
            'v': [(x, a) for x in mp.linspace(0, h, M+4)[1:-1]],
        }
    }

a = plate2['a']
get_coeff_func = plate2['get_coeff_func']

'''
_conditions = {border: funcs for border, funcs in conditions.items() if border != Borders.LEFT}

for func_name, two_eq_borders in [('u', (Borders.UPPER, Borders.BOTTOM)), ('v', (Borders.RIGHT, Borders.LEFT))]:
    coeff_func = get_coeff_func(func_name)
    this_funcs = {border: funcs[func_name] for border, funcs in _conditions.items() if func_name in funcs}
    for border in two_eq_borders:
        if border in this_funcs:
            for dot in get_corner_dots(h, a, border):
                A.append([0, ] * (8*M + 4) + coeff_func(*dot))
                b.append(this_funcs[border](dot[0] if border.is_vertical() else dot[1] + plate2['a']))
            break
    else:
        for border in this_funcs:
            dot = ({(0, 0), (h, a)} & get_corner_dots(h, a, border)).pop()
            A.append([0, ] * (8*M + 4) + coeff_func(*dot))
            b.append(this_funcs[border](dot[0] if border.is_vertical() else dot[1] + plate2['a']))

if len(A) == 3:
    border, tau_func = ({border: funcs['tau'] for border, funcs in _conditions.items() if 'tau' in funcs}).popitem()
    dot = get_corner_dots(h, a, border).pop()
    A.append([0, ] * (8*M + 4) + get_coeff_func('tau')(*dot))
    b.append(tau_func(dot[0] if border == Borders.RIGHT else dot[1]))
'''

'''
A.append([0, ] * (8*M + 4) + get_coeff_func('u')(0, a)); b.append(conditions[Borders.RIGHT]['u'](0, a1 + a2))
A.append([0, ] * (8*M + 4) + get_coeff_func('v')(0, a)); b.append(conditions[Borders.RIGHT]['v'](0, a1 + a2))
A.append([0, ] * (8*M + 4) + get_coeff_func('v')(h, a)); b.append(conditions[Borders.RIGHT]['v'](h, a1 + a2))
A.append([0, ] * (8*M + 4) + get_coeff_func('tau')(0, a)); b.append(conditions[Borders.UPPER]['tau'](0, a1 + a2))
'''

# for row in A:
#     print(' '*15, *[f'{float(i):=14e}' for i in row])

dots = find_dots(h, a, M)
coeff_func = {name: get_coeff_func(name) for name in ('u', 'v', 'sigma_x', 'sigma_y', 'tau')}
for border in filter(lambda b: b != Borders.LEFT, Borders):
    for name, inital_func in conditions[border].items():
        for dot in dots[border][name]:
            A.append([0, ] * (8*M + 4) + coeff_func[name](*dot))
            b.append(inital_func(dot[0] if border.is_vertical() else dot[1] + plate2['a']))

            # print(f'{border.name:6} {name:8}:', end='')
            # print(*[f'{float(i):=14e}' for i in A[-1]], '|', b[-1])

####################################### equations from touching border: ################################################

print('for both plates:')

a = a1 + a2

considered_func_names = ['u', 'v', 'sigma_y', 'tau']
x_dots = {name: mp.linspace(0, h, M+2)[1:-1] for name in considered_func_names}
x_dots

x_dots = mp.linspace(0, h, M+2)[1:-1]

for name in considered_func_names:
    for x in x_dots:
        coeffs1 = plate1['get_coeff_func'](name)(x, plate1['a'])
        coeffs2 = [-c for c in plate2['get_coeff_func'](name)(x, 0)]

        A.append(coeffs1 + coeffs2)
        b.append(0)
        # print(f'{"":6} {name:8}:', end='')
        # print(*[f'{float(i):=14e}' for i in A[-1]], '|', b[-1])

if report:
    print(f"(done in {timer})")

########################################################################################################################
############################################## CREATING SOLUTION: ######################################################
########################################################################################################################

if report:
    print("Solving system...")
    timer.start()

coefficients = mp.lu_solve(A, b)

if report:
    print(f"(done in {timer})")


def solution_func(x, y, func_name):
    plate = plate1

    if y > a1:
        y -= a1
        plate = plate2

    coeffs_for_coeffs = plate['get_coeff_func'](func_name)(x, y)
    coeffs = coefficients[(8*M+4):] if plate == plate2 else coefficients[:(8*M+4)]

    return sum(cc*c for cc, c in zip(coeffs_for_coeffs, coeffs))

solution = [
    lambda x, y, func_name=func_name: solution_func(x, y, func_name)
    for func_name in ['u', 'v', 'sigma_x', 'sigma_y', 'tau']
]

########################################################################################################################
################################################# PLOTTING GRAPH: ######################################################
########################################################################################################################
a = a1 + a2

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

print(f'(done in {timer})')
plt.show()
