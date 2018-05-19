from orthotropy.coeff_functions import get_coeff_function as orthotropy_get_coeff_function
from isotropy.coeff_functions import get_coeff_function as isotropy_get_coeff_function
from orthotropy.method import Timer, find_dots, get_corner_dots
from functools import lru_cache, partial
from mpmath import MPContext
from constants import Borders, MaterialTypes
import matplotlib.pyplot as plt

#               a1           a2           a3
#         <-----------><-----------><----------->
#          ___________  ___________  ___________ -----> y
#    ∧    |           ||           ||           |
#    |    |  plate1   ||  plate2   ||  plate3   |
#  h |    |           ||           ||           |
#    ∨    |           ||           ||           |
#        | ‾‾‾‾‾‾‾‾‾‾‾  ‾‾‾‾‾‾‾‾‾‾‾  ‾‾‾‾‾‾‾‾‾‾‾
#        |
#        |
#        ∨
#       x


M = 6
mp = MPContext()
mp.prec = 1500
h = 5
a1 = mp.mpf(10)/3
a2 = mp.mpf(10)/3
a3 = mp.mpf(10)/3
material_type = MaterialTypes.ORTHOTROPY


a = a1 + a2 + a3

plate1 = {'h': h, 'a': a1}
plate2 = {'h': h, 'a': a2}
plate3 = {'h': h, 'a': a3}

plates = [plate1, plate2, plate3]

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

h = mp.mpf(h)
a = mp.mpf(a)
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

################################# equations from conditions for plate1: ################################################

print('for plate1:')

dots = {
    Borders.UPPER: {
        'sigma_x': [(0, y) for y in mp.linspace(0, a1, M+2)[1:-1]],
        'tau': [(0, y) for y in mp.linspace(0, a1, M+2)[1:-1]],
    },
    Borders.BOTTOM: {
        'sigma_x': [(h, y) for y in mp.linspace(0, a1, M+2)[1:-1]],
        'tau': [(h, y) for y in mp.linspace(0, a1, M+2)[1:-1]],
    },
    Borders.LEFT: {
        'u': [(x, 0) for x in mp.linspace(0, h, M+3)[1:-1]],
        'v': [(x, 0) for x in mp.linspace(0, h, M+2)[1:-1]],
    }
}


get_coeff_func = plate1['get_coeff_func']

coeff_func = {name: get_coeff_func(name) for name in ('u', 'v', 'sigma_x', 'sigma_y', 'tau')}
for border in dots.keys():
    for name, inital_func in conditions[border].items():
        for x, y in dots[border][name]:
            A.append(coeff_func[name](x, y) + [0, ] * (8*M + 4) * 2)
            b.append(inital_func(x if border.is_vertical() else y))

            # print(f'{border.name:6} {name:8}:', end='')
            # print(*[f'{float(i):=14e}' for i in A[-1]], '|', b[-1])

################################# equations from conditions for plate2: ################################################

print('for plate2:')

dots = {
    Borders.UPPER: {
        'sigma_x': [(0, y) for y in mp.linspace(0, a2, M+2)[1:-1]],
        'tau': [(0, y) for y in mp.linspace(0, a2, M+2)[1:-1]],
    },
    Borders.BOTTOM: {
        'sigma_x': [(h, y) for y in mp.linspace(0, a2, M+2)[1:-1]],
        'tau': [(h, y) for y in mp.linspace(0, a2, M+2)[1:-1]],
    }
}

get_coeff_func = plate2['get_coeff_func']

coeff_func = {name: get_coeff_func(name) for name in ('u', 'v', 'sigma_x', 'sigma_y', 'tau')}
for border in dots.keys():
    for name, inital_func in conditions[border].items():
        for x, y in dots[border][name]:
            A.append([0, ] * (8*M + 4) + coeff_func[name](x, y) + [0, ] * (8*M + 4))
            b.append(inital_func(x if border.is_vertical() else y + a1))

            # print(f'{border.name:6} {name:8}:', end='')
            # print(*[f'{float(i):=14e}' for i in A[-1]], '|', b[-1])

####################################### equations from 1 touching border: ##############################################

print('for 1 touching border:')

considered_func_names = ['u', 'v', 'sigma_y', 'tau']

x_dots = {
    'u': mp.linspace(0, h, M+3)[1:-1],
    'v': mp.linspace(0, h, M+4)[1:-1],
    'sigma_y': mp.linspace(0, h, M+2)[1:-1],
    'tau': mp.linspace(0, h, M+3)[1:-1]
}

for name in considered_func_names:
    for x in x_dots[name]:
        coeffs1 = plate1['get_coeff_func'](name)(x, a1)
        coeffs2 = [-c for c in plate2['get_coeff_func'](name)(x, 0)]

        A.append(coeffs1 + coeffs2 + [0, ] * (8*M + 4))
        b.append(0)
        # print(f'{"":6} {name:8}:', end='')
        # print(*[f'{float(i):=14e}' for i in A[-1]], '|', b[-1])

################################# equations from conditions for plate3: ################################################

print('for plate3:')

dots = {
    Borders.UPPER: {
        'sigma_x': [(0, y) for y in mp.linspace(0, a3, M+2)[1:-1]],
        'tau': [(0, y) for y in mp.linspace(0, a3, M+2)[1:-1]],
    },
    Borders.BOTTOM: {
        'sigma_x': [(h, y) for y in mp.linspace(0, a3, M+2)[1:-1]],
        'tau': [(h, y) for y in mp.linspace(0, a3, M+2)[1:-1]],
    },
    Borders.RIGHT: {
        'u': [(x, a3) for x in mp.linspace(0, h, M+3)[1:-1]],
        'v': [(x, a3) for x in mp.linspace(0, h, M+4)[1:-1]],
    }
}

get_coeff_func = plate3['get_coeff_func']

coeff_func = {name: get_coeff_func(name) for name in ('u', 'v', 'sigma_x', 'sigma_y', 'tau')}
for border in dots.keys():
    for name, inital_func in conditions[border].items():
        for x, y in dots[border][name]:
            A.append([0, ] * (8*M + 4) * 2 + coeff_func[name](x, y))
            b.append(inital_func(x if border.is_vertical() else y + a1 + a2))

            # print(f'{border.name:6} {name:8}:', end='')
            # print(*[f'{float(i):=14e}' for i in A[-1]], '|', b[-1])

####################################### equations from 2 touching border: ##############################################

print('for 2 touching border:')

considered_func_names = ['u', 'v', 'sigma_y', 'tau']

x_dots = {
    'u': mp.linspace(0, h, M+3)[1:-1],
    'v': mp.linspace(0, h, M+4)[1:-1],
    'sigma_y': mp.linspace(0, h, M+2)[1:-1],
    'tau': mp.linspace(0, h, M+3)[1:-1]
}

for name in considered_func_names:
    for x in x_dots[name]:
        coeffs2 = plate2['get_coeff_func'](name)(x, a2)
        coeffs3 = [-c for c in plate3['get_coeff_func'](name)(x, 0)]

        A.append([0, ] * (8*M + 4) + coeffs2 + coeffs3)
        b.append(0)
        # print(f'{"":6} {name:8}:', end='')
        # print(*[f'{float(i):=14e}' for i in A[-1]], '|', b[-1])

if report:
    print(f"(done in {timer})")

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
        plate = plate1
        coeffs = coefficients[:(8*M+4)]
    elif y < a1 + a2:
        y -= a1
        plate = plate2
        coeffs = coefficients[(8*M+4):2*(8*M+4)]
    else:
        y -= a1 + a2
        plate = plate3
        coeffs = coefficients[2*(8*M+4):]

    coeffs_for_coeffs = plate['get_coeff_func'](func_name)(x, y)

    return sum(cc*c for cc, c in zip(coeffs_for_coeffs, coeffs))

solution = [
    lambda x, y, func_name=func_name: solution_func(x, y, func_name)
    for func_name in ['u', 'v', 'sigma_x', 'sigma_y', 'tau']
]

import pickle
file = open('coeffs.pickle', 'wb')
pickle.dump([str(c) for c in coefficients], file)
file.close()


########################################################################################################################
################################################# PLOTTING GRAPH: ######################################################
########################################################################################################################
#
# sigma_x = solution[2]
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
#         [sigma_x(x, y) for y in y_dots]
#     )
#
# for vals in list_of_vals:
#     plt.plot(y_dots, vals)
#
# # tau = solution[4]
# #
# # timer = Timer()
# # print("Plotting graph...")
# # timer.start()
# #
# # y_dots = mp.linspace(0, a, 200)
# #
# # list_of_vals = []
# # for x in [i*h/4 for i in range(5)]:
# #     list_of_vals.append(
# #         [tau(x, y) for y in y_dots]
# #     )
# #
# # for vals in list_of_vals:
# #     plt.plot(y_dots, vals)
#
# print(f'(done in {timer})')
# plt.show()
#
