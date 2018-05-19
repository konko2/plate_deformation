from orthotropy.coeff_functions import get_coeff_function as orthotropy_get_coeff_function
from isotropy.coeff_functions import get_coeff_function as isotropy_get_coeff_function
from orthotropy.method import Timer, find_dots, get_corner_dots
from functools import lru_cache, partial
from mpmath import MPContext
from constants import Borders, MaterialTypes
import matplotlib.pyplot as plt


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


import pickle
f = open('coeffs.pickle', 'rb')
coefficients = [mp.mpf(str_c) for str_c in pickle.load(f)]
f.close()

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


