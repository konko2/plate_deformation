from orthotropy.coeff_functions import get_coeff_function as orthotropy_get_coeff_function
from isotropy.coeff_functions import get_coeff_function as isotropy_get_coeff_function
from orthotropy.method import Timer, find_dots, get_corner_dots
from functools import lru_cache, partial
from mpmath import MPContext
from constants import Borders, MaterialTypes
import matplotlib.pyplot as plt

#                a1           a2           a3
#          <-----------><-----------><----------->
#           ___________  ___________  ___________ -----> y
#     ∧    |           ||           ||           |
#     |    |  plate11  ||  plate12  ||  plate13  |
#  h1 |    |           ||           ||           |
#     ∨    |___________||___________||___________|
#     ∧    |‾‾‾‾‾‾‾‾‾‾‾||‾‾‾‾‾‾‾‾‾‾‾||‾‾‾‾‾‾‾‾‾‾‾|
#     |    |  plate21  ||  plate22  ||  plate23  |
#  h2 |    |           ||           ||           |
#     ∨    |___________||___________||___________|
#     ∧    |‾‾‾‾‾‾‾‾‾‾‾||‾‾‾‾‾‾‾‾‾‾‾||‾‾‾‾‾‾‾‾‾‾‾|
#     |    |  plate31  ||  plate32  ||  plate33  |
#  h3 |    |           ||           ||           |
#     ∨    |           ||           ||           |
#         | ‾‾‾‾‾‾‾‾‾‾‾  ‾‾‾‾‾‾‾‾‾‾‾  ‾‾‾‾‾‾‾‾‾‾‾
#         |
#         ∨
#        x
#

M = 3
mp = MPContext()
mp.prec = 1500
h1 = h2 = h3 = mp.mpf(5)/3
a1 = a2 = a3 = mp.mpf(10)/3

material_type = MaterialTypes.ISOTROPY


a = a1 + a2 + a3
h = h1 + h2 + h3

plate11 = {'h': h1, 'a': a1}
plate12 = {'h': h1, 'a': a2}
plate13 = {'h': h1, 'a': a3}
plate21 = {'h': h2, 'a': a1}
plate22 = {'h': h2, 'a': a2}
plate23 = {'h': h2, 'a': a3}
plate31 = {'h': h3, 'a': a1}
plate32 = {'h': h3, 'a': a2}
plate33 = {'h': h3, 'a': a3}

plates = [
    plate11, plate12, plate13,
    plate21, plate22, plate23,
    plate31, plate32, plate33
]

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

    for plate in plates:
        plate['E_x'] = E_x
        plate['E_y'] = E_y
        plate['G_xy'] = G_xy
        plate['nu_xy'] = nu_xy

elif material_type == MaterialTypes.ISOTROPY:
    # Алюминий
    E = 70
    nu = 0.34

    for plate in plates:
        plate['E'] = E
        plate['nu'] = nu

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

for plate in plates:
    for key in plate:
        plate[key] = mp.mpf(plate[key])

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

for plate in plates:
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

################################# equations from conditions for plate11: ###############################################

print('for plate11:')

dots = {
        Borders.UPPER: {
            'sigma_x': [(0, y) for y in mp.linspace(0, a1, M+2)[1:-1]],
            'tau': [(0, y) for y in mp.linspace(0, a1, M+2)[1:-1]],
        },
        Borders.LEFT: {
            'u': [(x, 0) for x in mp.linspace(0, h1, M+3)[1:-1]],
            'v': [(x, 0) for x in mp.linspace(0, h1, M+2)[1:-1]],
        }
    }

get_coeff_func = plate11['get_coeff_func']

coeff_func = {name: get_coeff_func(name) for name in ('u', 'v', 'sigma_x', 'sigma_y', 'tau')}
for border in dots.keys():
    for name, inital_func in conditions[border].items():
        for x, y in dots[border][name]:
            A.append(coeff_func[name](x, y) + [0, ] * (8*M + 4) * 8)
            b.append(inital_func(x if border.is_vertical() else y))

            # print(f'{border.name:6} {name:8}:', end='')
            # print(*[f'{float(i):=14e}' for i in A[-1]], '|', b[-1])

################################# equations from conditions for plate12: ###############################################

print('for plate12:')

dots = {
        Borders.UPPER: {
            'sigma_x': [(0, y) for y in mp.linspace(0, a2, M+2)[1:-1]],
            'tau': [(0, y) for y in mp.linspace(0, a2, M+2)[1:-1]],
        }
    }

get_coeff_func = plate12['get_coeff_func']

coeff_func = {name: get_coeff_func(name) for name in ('u', 'v', 'sigma_x', 'sigma_y', 'tau')}
for border in dots.keys():
    for name, inital_func in conditions[border].items():
        for x, y in dots[border][name]:
            A.append([0, ] * (8*M + 4) + coeff_func[name](x, y) + [0, ] * (8*M + 4) * 7)
            b.append(inital_func(x if border.is_vertical() else y + a1))

            # print(f'{border.name:6} {name:8}:', end='')
            # print(*[f'{float(i):=14e}' for i in A[-1]], '|', b[-1])

################################# equations from conditions for plate13: ###############################################

print('for plate13:')

dots = {
        Borders.UPPER: {
            'sigma_x': [(0, y) for y in mp.linspace(0, a3, M+2)[1:-1]],
            'tau': [(0, y) for y in mp.linspace(0, a3, M+2)[1:-1]],
        },
        Borders.RIGHT: {
            'u': [(x, 0) for x in mp.linspace(0, h1, M+3)[1:-1]],
            'v': [(x, 0) for x in mp.linspace(0, h1, M+4)[1:-1]],
        }
    }

get_coeff_func = plate13['get_coeff_func']

coeff_func = {name: get_coeff_func(name) for name in ('u', 'v', 'sigma_x', 'sigma_y', 'tau')}
for border in dots.keys():
    for name, inital_func in conditions[border].items():
        for x, y in dots[border][name]:
            A.append([0, ] * (8*M + 4) * 2 + coeff_func[name](x, y) + [0, ] * (8*M + 4) * 6)
            b.append(inital_func(x if border.is_vertical() else y + a1 + a2))

            # print(f'{border.name:6} {name:8}:', end='')
            # print(*[f'{float(i):=14e}' for i in A[-1]], '|', b[-1])

################################# equations from conditions for plate21: ###############################################

print('for plate21:')

dots = {
        Borders.LEFT: {
            'u': [(x, 0) for x in mp.linspace(0, h2, M+3)[1:-1]],
            'v': [(x, 0) for x in mp.linspace(0, h2, M+2)[1:-1]],
        }
    }

get_coeff_func = plate21['get_coeff_func']

coeff_func = {name: get_coeff_func(name) for name in ('u', 'v', 'sigma_x', 'sigma_y', 'tau')}
for border in dots.keys():
    for name, inital_func in conditions[border].items():
        for x, y in dots[border][name]:
            A.append([0, ] * (8*M + 4) * 3 + coeff_func[name](x, y) + [0, ] * (8*M + 4) * 5)
            b.append(inital_func(x + h1 if border.is_vertical() else y))

            print(f'{border.name:6} {name:8}:', end='')
            print(*[f'{float(i):=14e}' for i in A[-1]], '|', b[-1])

################################# equations from conditions for plate23: ###############################################

print('for plate23:')

dots = {
        Borders.RIGHT: {
            'u': [(x, 0) for x in mp.linspace(0, h2, M+3)[1:-1]],
            'v': [(x, 0) for x in mp.linspace(0, h2, M+4)[1:-1]],
        }
    }

get_coeff_func = plate23['get_coeff_func']

coeff_func = {name: get_coeff_func(name) for name in ('u', 'v', 'sigma_x', 'sigma_y', 'tau')}
for border in dots.keys():
    for name, inital_func in conditions[border].items():
        for x, y in dots[border][name]:
            A.append([0, ] * (8*M + 4) * 5 + coeff_func[name](x, y) + [0, ] * (8*M + 4) * 3)
            b.append(inital_func(x + h1 if border.is_vertical() else y + a1 + a2))

            print(f'{border.name:6} {name:8}:', end='')
            print(*[f'{float(i):=14e}' for i in A[-1]], '|', b[-1])

################################# equations from conditions for plate31: ###############################################

print('for plate31:')

dots = {
        Borders.BOTTOM: {
            'sigma_x': [(0, y) for y in mp.linspace(0, a1, M+2)[1:-1]],
            'tau': [(0, y) for y in mp.linspace(0, a1, M+2)[1:-1]],
        },
        Borders.LEFT: {
            'u': [(x, 0) for x in mp.linspace(0, h3, M+3)[1:-1]],
            'v': [(x, 0) for x in mp.linspace(0, h3, M+2)[1:-1]],
        }
    }

get_coeff_func = plate31['get_coeff_func']

coeff_func = {name: get_coeff_func(name) for name in ('u', 'v', 'sigma_x', 'sigma_y', 'tau')}
for border in dots.keys():
    for name, inital_func in conditions[border].items():
        for x, y in dots[border][name]:
            A.append([0, ] * (8*M + 4) * 6 + coeff_func[name](x, y) + [0, ] * (8*M + 4) * 2)
            b.append(inital_func(x + h1 + h2 if border.is_vertical() else y))

            print(f'{border.name:6} {name:8}:', end='')
            print(*[f'{float(i):=14e}' for i in A[-1]], '|', b[-1])

################################# equations from conditions for plate32: ###############################################

print('for plate32:')

dots = {
        Borders.BOTTOM: {
            'sigma_x': [(0, y) for y in mp.linspace(0, a2, M+2)[1:-1]],
            'tau': [(0, y) for y in mp.linspace(0, a2, M+2)[1:-1]],
        }
    }

get_coeff_func = plate32['get_coeff_func']

coeff_func = {name: get_coeff_func(name) for name in ('u', 'v', 'sigma_x', 'sigma_y', 'tau')}
for border in dots.keys():
    for name, inital_func in conditions[border].items():
        for x, y in dots[border][name]:
            A.append([0, ] * (8*M + 4) * 7 + coeff_func[name](x, y) + [0, ] * (8*M + 4))
            b.append(inital_func(x + h1 + h2 if border.is_vertical() else y + a1))

            print(f'{border.name:6} {name:8}:', end='')
            print(*[f'{float(i):=14e}' for i in A[-1]], '|', b[-1])

################################# equations from conditions for plate33: ###############################################

print('for plate33:')

dots = {
    Borders.BOTTOM: {
        'sigma_x': [(0, y) for y in mp.linspace(0, a3, M + 2)[1:-1]],
        'tau': [(0, y) for y in mp.linspace(0, a3, M + 2)[1:-1]],
    },
    Borders.RIGHT: {
        'u': [(x, 0) for x in mp.linspace(0, h3, M + 3)[1:-1]],
        'v': [(x, 0) for x in mp.linspace(0, h3, M + 4)[1:-1]],
    }
}

get_coeff_func = plate33['get_coeff_func']

coeff_func = {name: get_coeff_func(name) for name in ('u', 'v', 'sigma_x', 'sigma_y', 'tau')}
for border in dots.keys():
    for name, inital_func in conditions[border].items():
        for x, y in dots[border][name]:
            A.append([0, ] * (8 * M + 4) * 8 + coeff_func[name](x, y))
            b.append(inital_func(x + h1 + h2 if border.is_vertical() else y + a1 + a2))

            print(f'{border.name:6} {name:8}:', end='')
            print(*[f'{float(i):=14e}' for i in A[-1]], '|', b[-1])

################################## equations for plate11 & plate12: ####################################################

print('for plate11 & plate12:')

considered_func_names = ['u', 'v', 'sigma_y', 'tau']

x_dots = {
    'u': mp.linspace(0, h1, M+3)[1:-1],
    'v': mp.linspace(0, h1, M+4)[1:-1],
    'sigma_y': mp.linspace(0, h1, M+2)[1:-1],
    'tau': mp.linspace(0, h1, M+3)[1:-1]
}

for name in considered_func_names:
    for x in x_dots[name]:
        coeffs11 = plate11['get_coeff_func'](name)(x, a1)
        coeffs12 = [-c for c in plate12['get_coeff_func'](name)(x, 0)]

        A.append(coeffs11 + coeffs12 + [0, ] * (8*M + 4) * 7)
        b.append(0)
        # print(f'{"":6} {name:8}:', end='')
        # print(*[f'{float(i):=14e}' for i in A[-1]], '|', b[-1])

################################## equations for plate12 & plate13: ####################################################

print('for plate12 & plate13:')

considered_func_names = ['u', 'v', 'sigma_y', 'tau']

x_dots = {
    'u': mp.linspace(0, h1, M+3)[1:-1],
    'v': mp.linspace(0, h1, M+4)[1:-1],
    'sigma_y': mp.linspace(0, h1, M+2)[1:-1],
    'tau': mp.linspace(0, h1, M+3)[1:-1]
}

for name in considered_func_names:
    for x in x_dots[name]:
        coeffs12 = plate12['get_coeff_func'](name)(x, a2)
        coeffs13 = [-c for c in plate13['get_coeff_func'](name)(x, 0)]

        A.append([0, ] * (8*M + 4) + coeffs12 + coeffs13 + [0, ] * (8*M + 4) * 6)
        b.append(0)
        # print(f'{"":6} {name:8}:', end='')
        # print(*[f'{float(i):=14e}' for i in A[-1]], '|', b[-1])

################################## equations for plate21 & plate22: ####################################################

print('for plate21 & plate22:')

considered_func_names = ['u', 'v', 'sigma_y', 'tau']

x_dots = {
    'u': mp.linspace(0, h2, M+3)[1:-1],
    'v': mp.linspace(0, h2, M+4)[1:-1],
    'sigma_y': mp.linspace(0, h2, M+2)[1:-1],
    'tau': mp.linspace(0, h2, M+3)[1:-1]
}

for name in considered_func_names:
    for x in x_dots[name]:
        coeffs21 = plate21['get_coeff_func'](name)(x, a1)
        coeffs22 = [-c for c in plate22['get_coeff_func'](name)(x, 0)]

        A.append([0, ] * (8*M + 4) * 3 + coeffs21 + coeffs22 + [0, ] * (8*M + 4) * 4)
        b.append(0)
        # print(f'{"":6} {name:8}:', end='')
        # print(*[f'{float(i):=14e}' for i in A[-1]], '|', b[-1])

################################## equations for plate22 & plate23: ####################################################

print('for plate22 & plate23:')

considered_func_names = ['u', 'v', 'sigma_y', 'tau']

x_dots = {
    'u': mp.linspace(0, h2, M+3)[1:-1],
    'v': mp.linspace(0, h2, M+4)[1:-1],
    'sigma_y': mp.linspace(0, h2, M+2)[1:-1],
    'tau': mp.linspace(0, h2, M+3)[1:-1]
}

for name in considered_func_names:
    for x in x_dots[name]:
        coeffs22 = plate22['get_coeff_func'](name)(x, a2)
        coeffs23 = [-c for c in plate23['get_coeff_func'](name)(x, 0)]

        A.append([0, ] * (8*M + 4) * 4 + coeffs22 + coeffs23 + [0, ] * (8*M + 4) * 3)
        b.append(0)
        # print(f'{"":6} {name:8}:', end='')
        # print(*[f'{float(i):=14e}' for i in A[-1]], '|', b[-1])

################################## equations for plate31 & plate32: ####################################################

print('for plate31 & plate32:')

considered_func_names = ['u', 'v', 'sigma_y', 'tau']

x_dots = {
    'u': mp.linspace(0, h3, M + 3)[1:-1],
    'v': mp.linspace(0, h3, M + 4)[1:-1],
    'sigma_y': mp.linspace(0, h3, M + 2)[1:-1],
    'tau': mp.linspace(0, h3, M + 3)[1:-1]
}

for name in considered_func_names:
    for x in x_dots[name]:
        coeffs31 = plate31['get_coeff_func'](name)(x, a1)
        coeffs32 = [-c for c in plate32['get_coeff_func'](name)(x, 0)]

        A.append([0, ] * (8 * M + 4) * 6 + coeffs31 + coeffs32 + [0, ] * (8 * M + 4) * 1)
        b.append(0)
        # print(f'{"":6} {name:8}:', end='')
        # print(*[f'{float(i):=14e}' for i in A[-1]], '|', b[-1])

################################## equations for plate32 & plate33: ####################################################

print('for plate32 & plate33:')

considered_func_names = ['u', 'v', 'sigma_y', 'tau']

x_dots = {
    'u': mp.linspace(0, h3, M + 3)[1:-1],
    'v': mp.linspace(0, h3, M + 4)[1:-1],
    'sigma_y': mp.linspace(0, h3, M + 2)[1:-1],
    'tau': mp.linspace(0, h3, M + 3)[1:-1]
}

for name in considered_func_names:
    for x in x_dots[name]:
        coeffs32 = plate32['get_coeff_func'](name)(x, a2)
        coeffs33 = [-c for c in plate33['get_coeff_func'](name)(x, 0)]

        A.append([0, ] * (8 * M + 4) * 7 + coeffs32 + coeffs33)
        b.append(0)
        # print(f'{"":6} {name:8}:', end='')
        # print(*[f'{float(i):=14e}' for i in A[-1]], '|', b[-1])
################################## equations for plate11 & plate21: ####################################################

print('for plate11 & plate21:')

considered_func_names = ['u', 'v', 'sigma_x', 'tau']

y_dots = mp.linspace(0, a1, M+2)[1:-1]

for name in considered_func_names:
    for y in y_dots:
        coeffs11 = plate11['get_coeff_func'](name)(h1, y)
        coeffs21 = [-c for c in plate21['get_coeff_func'](name)(0, y)]

        A.append(coeffs11 + [0, ] * (8 * M + 4) * 2 + coeffs21 + [0, ] * (8 * M + 4) * 5)
        b.append(0)
        # print(f'{"":6} {name:8}:', end='')
        # print(*[f'{float(i):=14e}' for i in A[-1]], '|', b[-1])

################################## equations for plate12 & plate22: ####################################################

print('for plate12 & plate22:')

considered_func_names = ['u', 'v', 'sigma_x', 'tau']

y_dots = mp.linspace(0, a2, M + 2)[1:-1]

for name in considered_func_names:
    for y in y_dots:
        coeffs12 = plate12['get_coeff_func'](name)(h1, y)
        coeffs22 = [-c for c in plate22['get_coeff_func'](name)(0, y)]

        A.append([0, ] * (8 * M + 4) + coeffs12 + [0, ] * (8 * M + 4) * 2 + coeffs22 + [0, ] * (8 * M + 4) * 4)
        b.append(0)
        # print(f'{"":6} {name:8}:', end='')
        # print(*[f'{float(i):=14e}' for i in A[-1]], '|', b[-1])

################################## equations for plate13 & plate23: ####################################################

print('for plate13 & plate23:')

considered_func_names = ['u', 'v', 'sigma_x', 'tau']

y_dots = mp.linspace(0, a3, M + 2)[1:-1]

for name in considered_func_names:
    for y in y_dots:
        coeffs13 = plate13['get_coeff_func'](name)(h1, y)
        coeffs23 = [-c for c in plate23['get_coeff_func'](name)(0, y)]

        A.append([0, ] * (8 * M + 4) * 2 + coeffs13 + [0, ] * (8 * M + 4) * 2 + coeffs23 + [0, ] * (8 * M + 4) * 3)
        b.append(0)
        # print(f'{"":6} {name:8}:', end='')
        # print(*[f'{float(i):=14e}' for i in A[-1]], '|', b[-1])
################################## equations for plate21 & plate31: ####################################################

print('for plate21 & plate31:')

considered_func_names = ['u', 'v', 'sigma_x', 'tau']

y_dots = mp.linspace(0, a1, M + 2)[1:-1]

for name in considered_func_names:
    for y in y_dots:
        coeffs21 = plate21['get_coeff_func'](name)(h2, y)
        coeffs31 = [-c for c in plate31['get_coeff_func'](name)(0, y)]

        A.append([0, ] * (8 * M + 4) * 3 + coeffs21 + [0, ] * (8 * M + 4) * 2 + coeffs31 + [0, ] * (8 * M + 4) * 2)
        b.append(0)
        # print(f'{"":6} {name:8}:', end='')
        # print(*[f'{float(i):=14e}' for i in A[-1]], '|', b[-1])

################################## equations for plate22 & plate32: ####################################################

print('for plate22 & plate32:')

considered_func_names = ['u', 'v', 'sigma_x', 'tau']

y_dots = mp.linspace(0, a2, M + 2)[1:-1]

for name in considered_func_names:
    for y in y_dots:
        coeffs22 = plate22['get_coeff_func'](name)(h2, y)
        coeffs32 = [-c for c in plate32['get_coeff_func'](name)(0, y)]

        A.append([0, ] * (8 * M + 4) * 4 + coeffs22 + [0, ] * (8 * M + 4) * 2 + coeffs32 + [0, ] * (8 * M + 4))
        b.append(0)
        # print(f'{"":6} {name:8}:', end='')
        # print(*[f'{float(i):=14e}' for i in A[-1]], '|', b[-1])

################################## equations for plate23 & plate33: ####################################################

print('for plate23 & plate33:')

considered_func_names = ['u', 'v', 'sigma_x', 'tau']

y_dots = mp.linspace(0, a3, M + 2)[1:-1]

for name in considered_func_names:
    for y in y_dots:
        coeffs23 = plate23['get_coeff_func'](name)(h2, y)
        coeffs33 = [-c for c in plate33['get_coeff_func'](name)(0, y)]

        A.append([0, ] * (8 * M + 4) * 5 + coeffs23 + [0, ] * (8 * M + 4) * 2 + coeffs33)
        b.append(0)
        # print(f'{"":6} {name:8}:', end='')
        # print(*[f'{float(i):=14e}' for i in A[-1]], '|', b[-1])


if report:
    print(f"(done in {timer})")

for cached_func in mp.sin, mp.cos, mp.sinh, mp.cosh:
    cached_func.cache_clear()

########################################################################################################################
############################################## CREATING SOLUTION: ######################################################
########################################################################################################################
print(len(A), len(A[0]))

if report:
    print("Solving system...")
    timer.start()

coefficients = mp.lu_solve(A, b)

if report:
    print(f"(done in {timer})")


def solution_func(x, y, func_name):
    if y < a1:
        col = 0
    elif y < a1 + a2:
        y -= a1
        col = 1
    else:
        y -= a1 + a2
        col = 2

    if x < h1:
        row = 0
    elif x < h1 + h2:
        x -= h1
        row = 1
    else:
        x -= h1 + h2
        row = 2

    plate = {
        (0, 0): plate11,
        (0, 1): plate12,
        (0, 2): plate13,
        (1, 0): plate21,
        (1, 1): plate22,
        (1, 2): plate23,
        (2, 0): plate31,
        (2, 1): plate32,
        (2, 2): plate33
    }[row, col]

    _begining = (8 * M + 4) * col + (8 * M + 4) * 3 * row
    coeffs = coefficients[_begining:_begining + (8 * M + 4)]

    coeffs_for_coeffs = plate['get_coeff_func'](func_name)(x, y)

    return sum(cc*c for cc, c in zip(coeffs_for_coeffs, coeffs))

solution = [
    lambda x, y, func_name=func_name: solution_func(x, y, func_name)
    for func_name in ['u', 'v', 'sigma_x', 'sigma_y', 'tau']
]

########################################################################################################################
################################################# PLOTTING GRAPH: ######################################################
########################################################################################################################

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

# tau = solution[4]
#
# timer = Timer()
# print("Plotting graph...")
# timer.start()
#
# y_dots = mp.linspace(0, a, 200)
#
# list_of_vals = []
# for x in [i*h/4 for i in range(5)]:
#     list_of_vals.append(
#         [tau(x, y) for y in y_dots]
#     )
#
# for vals in list_of_vals:
#     plt.plot(y_dots, vals)

print(f'(done in {timer})')
plt.show()

