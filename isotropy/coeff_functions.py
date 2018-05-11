from functools import partial


def create_L_funcs(nu, E, alpha, beta, mp):
    L = dict()
    for name in ('u', 'v', 'sigma_x', 'sigma_y', 'tau'):
        L[name] = [[None, ] * 4, [None, ] * 4]

    L['u'][0][0] = lambda y: [
        (
            (
                -y * alpha_i * mp.sinh(alpha_i * y)
                + 2 * (nu - 1) * mp.cosh(alpha_i * y)
            ) / (2 * (nu - 1))
        ) for alpha_i in alpha
    ]
    L['u'][0][1] = lambda y: [
        (
            - (
                (2 * nu - 1) * mp.sinh(alpha_i * y)
                - alpha_i * y * mp.cosh(alpha_i * y)
            ) / (2 * (nu - 1))
        ) for alpha_i in alpha
    ]
    L['u'][0][2] = lambda y: [
        (
            - y * mp.sinh(alpha_i * y) / (2 * (nu - 1)**2 * E)
        ) for alpha_i in alpha
    ]
    L['u'][0][3] = lambda y: [
        (
            - (
                (4 * nu - 3) * mp.sinh(alpha_i * y)
                - alpha_i * y * mp.cosh(alpha_i * y)
            ) / (2 * (nu - 1)**2 * E * alpha_i)
        ) for alpha_i in alpha
    ]

    L['v'][0][0] = lambda y: [
        (
            - (
                (2 * nu - 1) * mp.sinh(alpha_i * y)
                + alpha_i * y * mp.cosh(alpha_i * y)
            ) / (2 * (nu - 1))
        ) for alpha_i in alpha
    ]
    L['v'][0][1] = lambda y: [
        (
            (
                alpha_i * y * mp.sinh(alpha_i * y)
                + 2 * (nu - 1) * mp.cosh(alpha_i * y)
            ) / (2 * (nu - 1))
        ) for alpha_i in alpha
    ]
    L['v'][0][2] = lambda y: [
        (
            - (
                (4 * nu - 3) * mp.sinh(alpha_i * y)
                + alpha_i * y * mp.cosh(alpha_i * y)
            ) / (2 * (nu - 1)**2 * E * alpha_i)
        ) for alpha_i in alpha
    ]
    L['v'][0][3] = lambda y: [-i for i in L['u'][0][2](y)]

    L['sigma_x'][0][0] = lambda y: [
        (
            -E * alpha_i / 2 * (
                alpha_i * y * mp.sinh(alpha_i * y) + 2 * mp.cosh(alpha_i * y)
            )
        ) for alpha_i in alpha
    ]
    L['sigma_x'][0][1] = lambda y: [
        (
            E * alpha_i / 2 * (
                mp.sinh(alpha_i * y) + alpha_i * y * mp.cosh(alpha_i * y)
            )
        ) for alpha_i in alpha
    ]
    L['sigma_x'][0][2] = lambda y: [
        (
            - (
                alpha_i * y * mp.sinh(alpha_i * y)
                + 2 * nu * mp.cosh(alpha_i * y)
            ) / (2 * (nu - 1))
        ) for alpha_i in alpha
    ]
    L['sigma_x'][0][3] = lambda y: [
        (
            - (
                (2 * nu - 3) * mp.sinh(alpha_i * y)
                - alpha_i * y * mp.cosh(alpha_i * y)
            ) / (2 * (nu - 1))
        ) for alpha_i in alpha
    ]

    L['sigma_y'][0][0] = lambda y: [
        (
            alpha_i**2 * y * E * mp.sinh(alpha_i * y) / 2
        ) for alpha_i in alpha
    ]
    L['sigma_y'][0][1] = lambda y: [
        (
            E * alpha_i / 2 * (
                mp.sinh(alpha_i * y) - alpha_i * y * mp.cosh(alpha_i * y)
            )
        ) for alpha_i in alpha
    ]
    L['sigma_y'][0][2] = L['v'][0][1]
    L['sigma_y'][0][3] = lambda y: [-i for i in L['u'][0][1](y)]

    L['tau'][0][0] = L['sigma_x'][0][1]
    L['tau'][0][1] = lambda y: [-i for i in L['sigma_y'][0][0](y)]
    L['tau'][0][2] = lambda y: [-i for i in L['v'][0][0](y)]
    L['tau'][0][3] = L['u'][0][0]

    L['u'][1][0] = lambda x: [
        (
            (
                beta_i * x * mp.sinh(beta_i * x)
                + 2 * (nu - 1) * mp.cosh(beta_i * x)
            ) / (2 * (nu - 1))
        ) for beta_i in beta
    ]
    L['u'][1][1] = lambda x: [
        (
            - (
                (2 * nu - 1) * mp.sinh(beta_i * x)
                + beta_i * x * mp.cosh(beta_i * x)
            ) / (2 * (nu - 1))
        ) for beta_i in beta
    ]
    L['u'][1][2] = lambda x: [
        (
            - (
                (4 * nu - 3) * mp.sinh(beta_i * x)
                + beta_i * x * mp.cosh(beta_i * x)
            ) / (2 * (nu - 1)**2 * E * beta_i)
        ) for beta_i in beta
    ]
    L['u'][1][3] = lambda x: [
        (
            x * mp.sinh(beta_i * x) / (2 * (nu - 1)**2 * E)
        ) for beta_i in beta
    ]

    L['v'][1][0] = lambda x: [
        (
            - (
                (2 * nu - 1) * mp.sinh(beta_i * x)
                - beta_i * x * mp.cosh(beta_i * x)
            ) / (2 * (nu - 1))
        ) for beta_i in beta
    ]
    L['v'][1][1] = lambda x: [
        (
            (
                - beta_i * x * mp.sinh(beta_i * x)
                + 2 * (nu - 1) * mp.cosh(beta_i * x)
            ) / (2 * (nu - 1))
        ) for beta_i in beta
    ]
    L['v'][1][2] = lambda x: [-i for i in L['u'][1][3](x)]
    L['v'][1][3] = lambda x: [
        (
            - (
                (4 * nu - 3) * mp.sinh(beta_i * x)
                - beta_i * x * mp.cosh(beta_i * x)
            ) / (2 * (nu-1)**2 * E * beta_i)
        ) for beta_i in beta
    ]

    L['sigma_x'][1][0] = lambda x: [
        (
            E * beta_i / 2 * (
                mp.sinh(beta_i * x) - beta_i * x * mp.cosh(beta_i * x)
            )
        ) for beta_i in beta
    ]
    L['sigma_x'][1][1] = lambda x: [
        (
            E * beta_i**2 * x * mp.sinh(beta_i * x) / 2
        ) for beta_i in beta
    ]
    L['sigma_x'][1][2] = L['u'][1][0]
    L['sigma_x'][1][3] = lambda x: [-i for i in L['v'][1][0](x)]

    L['sigma_y'][1][0] = lambda x: [
        (
            E * beta_i / 2 * (
                mp.sinh(beta_i * x) + beta_i * x * mp.cosh(beta_i * x)
            )
        ) for beta_i in beta
    ]
    L['sigma_y'][1][1] = lambda x: [
        (
            - E * beta_i / 2 * (
                beta_i * x * mp.sinh(beta_i * x) + 2 * mp.cosh(beta_i * x)
            )
        ) for beta_i in beta
    ]
    L['sigma_y'][1][2] = lambda x: [
        (
            - (
                beta_i * x * mp.sinh(beta_i * x) + 2 * nu * mp.cosh(beta_i * x)
            ) / (2 * (nu - 1))
        ) for beta_i in beta
    ]
    L['sigma_y'][1][3] = lambda x: [
        (
            - (
                (2 * nu - 3) * mp.sinh(beta_i * x) - beta_i * x * mp.cosh(beta_i * x)
            ) / (2 * (nu - 1))
        ) for beta_i in beta
    ]

    L['tau'][1][0] = lambda x: [-i for i in L['sigma_x'][1][1](x)]
    L['tau'][1][1] = L['sigma_y'][1][0]
    L['tau'][1][2] = lambda x: [-i for i in L['u'][1][1](x)]
    L['tau'][1][3] = L['v'][1][1]

    return L


def coeff_u_func(x, y, alpha, beta, _L, E, nu, mp):
    result = list()
    L = [[L_i_func(var) for L_i_func in _L[i]] for i, var in enumerate([y, x])]

    result.extend([1, -2*y/(E*(nu-1))])
    for i, alpha_i in enumerate(alpha):
        result.extend([L[0][k][i] * mp.cos(alpha_i * x) for k in range(4)])
    result.extend([0, 0])
    for i, beta_i in enumerate(beta):
        result.extend([L[1][k][i] * mp.sin(beta_i * y) for k in range(4)])

    return result


def coeff_v_func(x, y, alpha, beta, _L, E, nu, mp):
    result = list()
    L = [[L_i_func(var) for L_i_func in _L[i]] for i, var in enumerate([y, x])]

    result.extend([0, 0])
    for i, alpha_i in enumerate(alpha):
        result.extend([L[0][k][i] * mp.sin(alpha_i * x) for k in range(4)])
    result.extend([1, -2*x/(E*(nu-1))])
    for i, beta_i in enumerate(beta):
        result.extend([L[1][k][i] * mp.cos(beta_i * y) for k in range(4)])

    return result


def coeff_sigma_func(x, y, alpha, beta, _L, E, nu, mp):
    result = list()
    L = [[L_i_func(var) for L_i_func in _L[i]] for i, var in enumerate([y, x])]

    result.extend([0, 0])
    for i, alpha_i in enumerate(alpha):
        result.extend([L[0][k][i] * mp.sin(alpha_i * x) for k in range(4)])
    result.extend([0, 0])
    for i, beta_i in enumerate(beta):
        result.extend([L[1][k][i] * mp.sin(beta_i * y) for k in range(4)])

    return result

coeff_sigma_x_func = coeff_sigma_func
coeff_sigma_y_func = coeff_sigma_func


def coeff_tau_func(x, y, alpha, beta, _L, E, nu, mp):
    result = list()
    L = [[L_i_func(var) for L_i_func in _L[i]] for i, var in enumerate([y, x])]

    result.extend([0, 1])
    for i, alpha_i in enumerate(alpha):
        result.extend([L[0][k][i] * mp.cos(alpha_i * x) for k in range(4)])
    result.extend([0, 1])
    for i, beta_i in enumerate(beta):
        result.extend([L[1][k][i] * mp.cos(beta_i * y) for k in range(4)])

    return result


def get_coeff_function(nu, E, alpha, beta, mp, func_name):
    L = create_L_funcs(nu, E, alpha, beta, mp)
    return partial(globals()[f'coeff_{func_name}_func'], alpha=alpha, beta=beta, _L=L[func_name], E=E, nu=nu, mp=mp)


# TODO:
# 1) Сделать аккуратнее
# 2) Добавить возможность считать используя конкретный контекст
#   В чем проблема?
#   Хотелось бы сослать все триг функции в одно место, чтобы все что здесь написано ссылалось на то место
#   то есть достать их из контекста в get_coeff_function а дальше не думать об этом и просто пользоваться
#   Т. к. может быть создано множество окон с графиками, использовать глобальное пространство имен нельзя
#   (или придется каждый раз переопределять имена на уровне модуля, это не красиво)
#   Поскольку create_L_funcs создается в глобальном пространстве имен, то дерево выглядит так:
#   локальное пр. имен -> глобальное пр. имен
#   То есть обычным поиском в дереве имен не обойтись. Какой может быть выход?
#   а) Можно передать контекст в аргументы каждой фунции, и дальше использовать функции контекста
#       +:
#       -:
#   б) Можно передавать триг функции в аргументы каждой функции.
#       +:
#       -:
#   в) Сделать функцию, которая будет выдавать триг функции, вызывать ее внутри функций этого модуля
#       +:
#       -:
#   г) Использовать синтаксис классов (т е для каждой задачи будет свой экземпляр этих функци)
#       +:
#       -:
#   д) Использовать декораторы функций (как замена синтаксиса классов)
#       +:
#       -:
#   е) Использовать модуль для хранения триг функций
#       +:
#       -:
#
# 3) Добавить к функциям кэш, отчистку кэша
#   Очищать кэш нужно в методе, после того как построена матрица (именно метод понимает когда матрица построена, а не
#   представленные в этом модуле функции). Значит у метода должен быть доступ к триг функциям, которые здесь
#   используются.
# 4) Сделать так, чтобы get_coeff_function возвращал словарь
#
