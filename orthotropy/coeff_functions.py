from functools import partial


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


def create_L_funcs(H, a_, d, alpha, beta, mp):
    L = dict()
    for name in ('u', 'v', 'sigma_x', 'sigma_y', 'tau'):
        L[name] = [[None, ] * 4, [None, ] * 4]

    L['u'][0][0] = lambda y: [
        H[6, 6] / d[0] * (
            + (H[1, 2] + H[2, 2] * a_[0][0]**2) * mp.cosh(a_[0][0] * alpha_i * y)
            - (H[1, 2] + H[2, 2] * a_[0][1]**2) * mp.cosh(a_[0][1] * alpha_i * y)
        ) for alpha_i in alpha
    ]
    L['u'][0][1] = lambda y: [
        - H[6, 6] / d[0] * (
            + (H[1, 2] + H[2, 2] * a_[0][0]**2) * mp.sinh(a_[0][0] * alpha_i * y) / a_[0][0]
            - (H[1, 2] + H[2, 2] * a_[0][1]**2) * mp.sinh(a_[0][1] * alpha_i * y) / a_[0][1]
        ) for alpha_i in alpha
    ]
    L['u'][0][2] = lambda y: [
        - (H[1, 2] + H[6, 6])/ (d[0] * alpha_i) * (
            mp.cosh(a_[0][0] * alpha_i * y) - mp.cosh(a_[0][1] * alpha_i * y)
        ) for alpha_i in alpha
    ]
    L['u'][0][3] = lambda y: [
        (
            - (H[6, 6] - H[2, 2] * a_[0][0]**2) * mp.sinh(a_[0][0] * alpha_i * y) / a_[0][0]
            + (H[6, 6] - H[2, 2] * a_[0][1]**2) * mp.sinh(a_[0][1] * alpha_i * y) / a_[0][1]
        ) / (d[0] * alpha_i) for alpha_i in alpha
    ]

    L['v'][0][0] = lambda y: [
        H[6, 6] / d[0] * (
            + (H[1, 1] + H[1, 2] * a_[0][0]**2) * mp.sinh(a_[0][0] * alpha_i * y) / a_[0][0]
            - (H[1, 1] + H[1, 2] * a_[0][1]**2) * mp.sinh(a_[0][1] * alpha_i * y) / a_[0][1]
        ) for alpha_i in alpha
    ]
    L['v'][0][1] = lambda y: [
        (
            (
                H[1, 2] * H[6, 6] + H[1, 2]**2 -
                H[1, 1] * H[2, 2] + H[2, 2] * H[6, 6] * a_[0][0]**2
            ) * mp.cosh(a_[0][0] * alpha_i * y)
            - (
                H[1, 2] * H[6, 6] + H[1, 2]**2 -
                H[1, 1] * H[2, 2] + H[2, 2] * H[6, 6] * a_[0][1]**2
            ) * mp.cosh(a_[0][1] * alpha_i * y)
        ) / d[0] for alpha_i in alpha
    ]
    L['v'][0][2] = lambda y: [
        (
            - (H[1, 1] - H[6, 6] * a_[0][0]**2) * mp.sinh(a_[0][0] * alpha_i * y) / a_[0][0]
            + (H[1, 1] - H[6, 6] * a_[0][1]**2) * mp.sinh(a_[0][1] * alpha_i * y) / a_[0][1]
        ) / (d[0] * alpha_i) for alpha_i in alpha
    ]
    L['v'][0][3] = lambda y: [-i for i in L['u'][0][2](y)]

    L['sigma_x'][0][0] = lambda y: [
        H[6, 6] * (H[1, 2]**2 - H[1, 1] * H[2, 2]) * alpha_i / d[0] * (
            + mp.cosh(a_[0][0] * alpha_i * y) * a_[0][0]**2
            - mp.cosh(a_[0][1] * alpha_i * y) * a_[0][1]**2
        ) for alpha_i in alpha
    ]
    L['sigma_x'][0][1] = lambda y: [
        (
            (
                H[1, 1] * H[1, 2] * H[6, 6]
                + a_[0][0]**2 * (
                    -H[1, 1] * H[1, 2] * H[2, 2] + H[6, 6] * H[1, 2]**2 +
                    H[1, 2]**3 + H[1, 1] * H[2, 2] * H[6, 6]
                )
                + a_[0][0]**4 * H[1, 2] * H[2, 2] * H[6, 6]
            ) * mp.sinh(a_[0][0] * alpha_i * y) / a_[0][0]
             - (
                H[1, 1] * H[1, 2] * H[6, 6]
                + a_[0][1]**2 * (
                    -H[1, 1] * H[1, 2] * H[2, 2] + H[6, 6] * H[1, 2]**2 +
                    H[1, 2]**3 + H[1, 1] * H[2, 2] * H[6, 6]
                )
                + a_[0][1]**4 * H[1, 2] * H[2, 2] * H[6, 6]
            ) * mp.sinh(a_[0][1] * alpha_i * y) / a_[0][1]
        ) * alpha_i / d[0] for alpha_i in alpha
    ]
    L['sigma_x'][0][2] = lambda y: [
        H[6, 6] / d[0] * (
            + (H[1, 1] + H[1, 2] * a_[0][0]**2) * mp.cosh(a_[0][0] * alpha_i * y)
            - (H[1, 1] + H[1, 2] * a_[0][1]**2) * mp.cosh(a_[0][1] * alpha_i * y)
        ) for alpha_i in alpha
    ]
    L['sigma_x'][0][3] = lambda y: [
        (
            (
                H[1, 1] * H[6, 6] + a_[0][0]**2 *
                (-H[1, 1] * H[2, 2] + H[1, 2] * H[6, 6] + H[1, 2]**2)
            ) * mp.sinh(a_[0][0] * alpha_i * y) / a_[0][0]
            - (
                H[1, 1] * H[6, 6] + a_[0][1]**2 *
                (-H[1, 1] * H[2, 2] + H[1, 2] * H[6, 6] + H[1, 2]**2)
            ) *  mp.sinh(a_[0][1] * alpha_i * y) / a_[0][1]
        ) / d[0] for alpha_i in alpha
    ]

    L['sigma_y'][0][0] = lambda y: [
        H[6, 6] * (H[1, 1] * H[2, 2] - H[1, 2]**2) * alpha_i / d[0] * (
            mp.cosh(a_[0][0] * alpha_i * y) - mp.cosh(a_[0][1] * alpha_i * y)
        ) for alpha_i in alpha
    ]
    L['sigma_y'][0][1] = lambda y: [
        (
            (
                H[6, 6] * H[1, 2]**2
                + a_[0][0]**2 * (
                    2 * H[1, 2] * H[2, 2] * H[6, 6]
                    - H[1, 1] * H[2, 2]**2 + H[2, 2] * H[1, 2]**2
                )
                + a_[0][0]**4 * H[6, 6] * H[2, 2]**2
            ) * mp.sinh(a_[0][0] * alpha_i * y) / a_[0][0]
            - (
                H[6, 6] * H[1, 2]**2
                + a_[0][1]**2 * (
                    2 * H[1, 2] * H[2, 2] * H[6, 6]
                    - H[1, 1] * H[2, 2]**2 + H[2, 2] * H[1, 2]**2
                )
                + a_[0][1]**4 * H[6, 6] * H[2, 2]**2
            ) * mp.sinh(a_[0][1] * alpha_i * y) / a_[0][1]
        ) * alpha_i / d[0] for alpha_i in alpha
    ]
    L['sigma_y'][0][2] = L['v'][0][1]
    L['sigma_y'][0][3] = lambda y: [-i for i in L['u'][0][1](y)]

    L['tau'][0][0] = lambda y: [
        H[6, 6]**2 * alpha_i / d[0] * (
            (
                H[1, 1] + 2 * H[1, 2] * a_[0][0]**2 + H[2, 2] * a_[0][0]**4
            ) * mp.sinh(a_[0][0] * alpha_i * y) / a_[0][0]
            - (
                H[1, 1] + 2 * H[1, 2] * a_[0][1]**2 + H[2, 2] * a_[0][1]**4
            ) * mp.sinh(a_[0][1] * alpha_i * y) / a_[0][1]
        ) for alpha_i in alpha
    ]
    L['tau'][0][1] = lambda y: [-i for i in L['sigma_y'][0][0](y)]
    L['tau'][0][2] = lambda y: [-i for i in L['v'][0][0](y)]
    L['tau'][0][3] = L['u'][0][0]

    L['u'][1][0] = lambda x: [
        (
            - (
                H[1, 1] * H[2, 2] - H[1, 2]**2 -
                H[1, 2] * H[6, 6] - H[1, 1] * H[6, 6] * a_[1][0]**2
            ) * mp.cosh(a_[1][0] * beta_i * x)
            + (
                H[1, 1] * H[2, 2] - H[1, 2]**2 -
                H[1, 2] * H[6, 6] - H[1, 1] * H[6, 6] * a_[1][1]**2
            ) * mp.cosh(a_[1][1] * beta_i * x)
        ) / d[1] for beta_i in beta
    ]
    L['u'][1][1] = lambda x: [
        H[6, 6] / d[1] * (
            + (H[2, 2] + H[1, 2] * a_[1][0]**2) * mp.sinh(a_[1][0] * beta_i * x) / a_[1][0]
            - (H[2, 2] + H[1, 2] * a_[1][1]**2) * mp.sinh(a_[1][1] * beta_i * x) / a_[1][1]
        ) for beta_i in beta
    ]
    L['u'][1][2] = lambda x: [
        (
            - (H[2, 2] - H[6, 6] * a_[1][0]**2) * mp.sinh(a_[1][0] * beta_i * x) / a_[1][0]
            + (H[2, 2] - H[6, 6] * a_[1][1]**2) * mp.sinh(a_[1][1] * beta_i * x) / a_[1][1]
        ) / (d[1] * beta_i) for beta_i in beta
    ]
    L['u'][1][3] = lambda x: [
        (H[1, 2] + H[6, 6]) / (beta_i * d[1]) * (
            mp.cosh(a_[1][0] * beta_i * x) - mp.cosh(a_[1][1] * beta_i * x)
        ) for beta_i in beta
    ]

    L['v'][1][0] = lambda x: [
        H[6, 6] * (
            - (H[1, 2] + H[1, 1] * a_[1][0]**2) * mp.sinh(a_[1][0] * beta_i * x) / a_[1][0]
            + (H[1, 2] + H[1, 1] * a_[1][1]**2) * mp.sinh(a_[1][1] * beta_i * x) / a_[1][1]
        ) / d[1] for beta_i in beta
    ]
    L['v'][1][1] = lambda x: [
        H[6, 6] / d[1] * (
            + (H[1, 2] + H[1, 1] * a_[1][0]**2) * mp.cosh(a_[1][0] * beta_i * x)
            - (H[1, 2] + H[1, 1] * a_[1][1]**2) * mp.cosh(a_[1][1] * beta_i * x)
        ) for beta_i in beta
    ]
    L['v'][1][2] = lambda x: [-i for i in L['u'][1][3](x)]
    L['v'][1][3] = lambda x: [
        (
            - (H[6, 6] - H[1, 1] * a_[1][0]**2) * mp.sinh(a_[1][0] * beta_i * x) / a_[1][0]
            + (H[6, 6] - H[1, 1] * a_[1][1]**2) * mp.sinh(a_[1][1] * beta_i * x) / a_[1][1]
        ) / (d[1] * beta_i) for beta_i in beta
    ]

    L['sigma_x'][1][0] = lambda x: [
        (
            (
                H[6, 6] * H[1, 2]**2
                + a_[1][0]**2 * (
                    2 * H[1, 1] * H[1, 2] * H[6, 6] +
                    H[1, 1] * H[1, 2]**2 - H[2, 2] * H[1, 1]**2
                )
                + a_[1][0]**4 * H[6, 6] * H[1, 1]**2
            ) * mp.sinh(a_[1][0] * beta_i * x) / a_[1][0]
            - (
                H[6, 6] * H[1, 2]**2
                + a_[1][1]**2 * (
                    2 * H[1, 1] * H[1, 2] * H[6, 6] +
                    H[1, 1] * H[1, 2]**2 - H[2, 2] * H[1, 1]**2
                )
                + a_[1][1]**4 * H[6, 6] * H[1, 1]**2
            ) * mp.sinh(a_[1][1] * beta_i * x) / a_[1][1]
        ) * beta_i / d[1] for beta_i in beta
    ]
    L['sigma_x'][1][1] = lambda x: [
        H[6, 6] * (H[1, 1] * H[2, 2] - H[1, 2]**2) * beta_i / d[1] * (
            mp.cosh(a_[1][0] * beta_i * x) - mp.cosh(a_[1][1] * beta_i * x)
        ) for beta_i in beta
    ]
    L['sigma_x'][1][2] = L['u'][1][0]
    L['sigma_x'][1][3] = lambda x: [-i for i in L['v'][1][0](x)]

    L['sigma_y'][1][0] = lambda x: [
        (
            (
                H[1, 2] * H[2, 2] * H[6, 6]
                + a_[1][0]**2 * (
                    H[1, 1] * H[2, 2] * H[6, 6] + H[6, 6] * H[1, 2]**2 +
                    H[1, 2]**3 - H[1, 1] * H[1, 2] * H[2, 2]
                )
                + a_[1][0]**4 * H[1, 1] * H[1, 2] * H[6, 6]
            ) * mp.sinh(a_[1][0] * beta_i * x) / a_[1][0]
            - (
                H[1, 2] * H[2, 2] * H[6, 6]
                + a_[1][1]**2 * (
                    H[1, 1] * H[2, 2] * H[6, 6] + H[6, 6] * H[1, 2]**2 +
                    H[1, 2]**3 - H[1, 1] * H[1, 2] * H[2, 2]
                )
                + a_[1][1]**4 * H[1, 1] * H[1, 2] * H[6, 6]
            ) * mp.sinh(a_[1][1] * beta_i * x) / a_[1][1]
        ) * beta_i / d[1] for beta_i in beta
    ]
    L['sigma_y'][1][1] = lambda x: [
        H[6, 6] * (H[1, 2]**2 - H[1, 1] * H[2, 2]) * beta_i / d[1] * (
            + mp.cosh(a_[1][0] * beta_i * x) * a_[1][0]**2
            - mp.cosh(a_[1][1] * beta_i * x) * a_[1][1]**2
        ) for beta_i in beta
    ]
    L['sigma_y'][1][2] = lambda x: [
        H[6, 6] / d[1] * (
            + (H[2, 2] + H[1, 2] * a_[1][0]**2) * mp.cosh(a_[1][0] * beta_i * x)
            - (H[2, 2] + H[1, 2] * a_[1][1]**2) * mp.cosh(a_[1][1] * beta_i * x)
        ) for beta_i in beta
    ]
    L['sigma_y'][1][3] = lambda x: [
        (
            (
                H[2, 2] * H[6, 6] + a_[1][0]**2 *
                (-H[1, 1] * H[2, 2] + H[1, 2] * H[6, 6] + H[1, 2]**2)
            ) * mp.sinh(a_[1][0] * beta_i * x) / a_[1][0]
            - (
                H[2, 2] * H[6, 6] + a_[1][1]**2 *
                (-H[1, 1] * H[2, 2] + H[1, 2] * H[6, 6] + H[1, 2]**2)
            ) * mp.sinh(a_[1][1] * beta_i * x) / a_[1][1]
        ) / d[1] for beta_i in beta
    ]

    L['tau'][1][0] = lambda x: [-i for i in L['sigma_x'][1][1](x)]
    L['tau'][1][1] = lambda x: [
        H[6, 6]**2 * (
            (
                H[2, 2] + 2 * H[1, 2] * a_[1][0]**2 + H[1, 1] * a_[1][0]**4
            ) * mp.sinh(a_[1][0] * beta_i * x) / a_[1][0]
            - (
                H[2, 2] + 2 * H[1, 2] * a_[1][1]**2 + H[1, 1] * a_[1][1]**4
            ) * mp.sinh(a_[1][1] * beta_i * x) / a_[1][1]
        ) * beta_i / d[1] for beta_i in beta
    ]
    L['tau'][1][2] = lambda x: [-i for i in L['u'][1][1](x)]
    L['tau'][1][3] = L['v'][1][1]

    return L


def coeff_u_func(x, y, H, a_, d, _L, alpha, beta, mp):
    result = list()
    L = [[L_i_func(var) for L_i_func in _L[i]] for i, var in enumerate([y, x])]

    result.append(H[2, 2] * H[6, 6] * (a_[0][0]**2 - a_[0][1]**2) / d[0])
    result.append(H[2, 2] * (a_[0][0]**2 - a_[0][1]**2) * y / d[0])
    for i, alpha_i in enumerate(alpha):
        result.extend([L[0][k][i] * mp.cos(alpha_i * x) for k in range(4)])
    result.extend([0, 0])
    for i, beta_i in enumerate(beta):
        result.extend([L[1][k][i] * mp.sin(beta_i * y) for k in range(4)])

    return result


def coeff_v_func(x, y, H, a_, d, _L, alpha, beta, mp):
    result = list()
    L = [[L_i_func(var) for L_i_func in _L[i]] for i, var in enumerate([y, x])]

    result.extend([0, 0])
    for i, alpha_i in enumerate(alpha):
        result.extend([L[0][k][i] * mp.sin(alpha_i * x) for k in range(4)])
    result.append(H[1, 1] * H[6, 6] * (a_[1][0]**2 - a_[1][1]**2) / d[1])
    result.append(H[1, 1] * (a_[1][0]**2 - a_[1][1]**2) * x / d[1])
    for i, beta_i in enumerate(beta):
        result.extend([L[1][k][i] * mp.cos(beta_i * y) for k in range(4)])

    return result


def coeff_sigma_func(x, y, H, a_, d, _L, alpha, beta, mp):
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


def coeff_tau_func(x, y, H, a_, d, _L, alpha, beta, mp):
    result = list()
    L = [[L_i_func(var) for L_i_func in _L[i]] for i, var in enumerate([y, x])]
    result.append(0)
    result.append(H[2, 2] * H[6, 6] * (a_[0][0]**2 - a_[0][1]**2) / d[0])
    for i, alpha_i in enumerate(alpha):
        result.extend([L[0][k][i] * mp.cos(alpha_i * x) for k in range(4)])

    result.append(0)
    result.append(H[1, 1] * H[6, 6] * (a_[1][0]**2 - a_[1][1]**2) / d[1])
    for i, beta_i in enumerate(beta):
        result.extend([L[1][k][i] * mp.cos(beta_i * y) for k in range(4)])

    return result


def get_coeff_function(E_x, E_y, G_xy, nu_xy, alpha, beta, mp, func_name):
    H = find_H(E_x, E_y, G_xy, nu_xy)
    a_, d = find_a_and_d(H)
    L = create_L_funcs(H, a_, d, alpha, beta, mp)
    return partial(
        globals()[f'coeff_{func_name}_func'],
        alpha=alpha,
        beta=beta,
        _L=L[func_name],
        H=H,
        a_=a_,
        d=d,
        mp=mp
    )

"""
def _wrapper_u(H, a_, d, alpha, beta):
    M = len(alpha)

    _c = dict()
    _c['h_0'] = {
        0: H[2, 2] * H[6, 6] * (a_[0][0]**2 - a_[0][1]**2) / d[0],
        3: H[2, 2] * (a_[0][0]**2 - a_[0][1]**2) / d[0]
    }

    _c['h_n'] = [
        [
            H[6, 6] * (H[1, 2] + H[2, 2] * a_[0][0]**2) / d[0],
            - H[6, 6] * (H[1, 2] + H[2, 2] * a_[0][1]**2) / d[0]
        ],
        [
            - H[6, 6] * (H[1, 2] + H[2, 2] * a_[0][0]**2) / (d[0] * a_[0][0]),
            H[6, 6] * (H[1, 2] + H[2, 2] * a_[0][1]**2) / (d[0] * a_[0][1])
        ],
        [
            - (H[1, 2] + H[6, 6]) / d[0],
            (H[1, 2] + H[6, 6]) / d[0]
        ],
        [
            (- H[6, 6] + H[2, 2] * a_[0][0]**2) / (d[0] * a_[0][0]),
            (H[6, 6] - H[2, 2] * a_[0][1]**2) / (d[0] * a_[0][1])
        ]
    ]

    _c['g_n'] = [
        [
            (
                -H[1, 1] * H[2, 2] + H[1, 2]**2 +
                H[1, 2] * H[6, 6] + H[1, 1] * H[6, 6] * a_[1][0]**2
            ) / d[1],
            (
                H[1, 1] * H[2, 2] - H[1, 2]**2 -
                H[1, 2] * H[6, 6] - H[1, 1] * H[6, 6] * a_[1][1]**2
            ) / d[1]
        ],
        [
            H[6, 6] * (H[2, 2] + H[1, 2] * a_[1][0]**2) / (d[1] * a_[1][0]),
            -H[6, 6] * (H[2, 2] + H[1, 2] * a_[1][1]**2) / (d[1] * a_[1][1])
        ],
        [
            (-H[2, 2] + H[6, 6] * a_[1][0]**2) / (d[1] * a_[1][0]),
            (H[2, 2] - H[6, 6] * a_[1][1]**2) / (d[1] * a_[1][1])
        ],
        [
            (H[1, 2] + H[6, 6]) / d[1],
            -(H[1, 2] + H[6, 6]) / d[1]
        ]
    ]

    def coeff_u_eq(x, y):
        result = [_c['h_0'][0], _c['h_0'][3]*y]

        if x == 0:
            for i in range(M):
                result.append(
                    _c['h_n'][0][0] * cosh(a_[0][0] * alpha[i] * y) +
                    _c['h_n'][0][1] * cosh(a_[0][1] * alpha[i] * y)
                )
                result.append(
                    _c['h_n'][1][0] * sinh(a_[0][0] * alpha[i] * y) +
                    _c['h_n'][1][1] * sinh(a_[0][1] * alpha[i] * y)
                )
                result.append(
                    (
                        _c['h_n'][2][0] * cosh(a_[0][0] * alpha[i] * y) +
                        _c['h_n'][2][1] * cosh(a_[0][1] * alpha[i] * y)
                    ) / alpha[i]
                )
                result.append(
                    (
                        _c['h_n'][3][0] * sinh(a_[0][0] * alpha[i] * y) +
                        _c['h_n'][3][1] * sinh(a_[0][1] * alpha[i] * y)
                    ) / alpha[i]
                )

            result.extend([0, 0])

            for i in range(M):
                result.extend(
                    [sum(_c['g_n'][0]) * sin(beta[i]*y), 0, 0, 0]
                )
        elif y == 0:
            for i in range(M):
                result.extend(
                    [sum(_c['h_n'][0]) * cos(alpha[i] * x), 0, 0, 0]
                )

            result.extend([0, ] * (2 + 4*M))
        else:
            for i in range(M):
                result.append(
                    (
                        _c['h_n'][0][0] * cosh(a_[0][0] * alpha[i] * y) +
                        _c['h_n'][0][1] * cosh(a_[0][1] * alpha[i] * y)
                    ) * cos(alpha[i]*x)
                )
                result.append(
                    (
                        _c['h_n'][1][0] * sinh(a_[0][0] * alpha[i] * y) +
                        _c['h_n'][1][1] * sinh(a_[0][1] * alpha[i] * y)
                    ) * cos(alpha[i]*x)
                )
                result.append(
                    (
                        _c['h_n'][2][0] * cosh(a_[0][0] * alpha[i] * y) +
                        _c['h_n'][2][1] * cosh(a_[0][1] * alpha[i] * y)
                    ) * cos(alpha[i]*x) / alpha[i]
                )
                result.append(
                    (
                        _c['h_n'][3][0] * sinh(a_[0][0] * alpha[i] * y) +
                        _c['h_n'][3][1] * sinh(a_[0][1] * alpha[i] * y)
                    ) * cos(alpha[i]*x) / alpha[i]
                )

            result.extend([0, 0])

            for i in range(M):
                result.append(
                    (
                        _c['g_n'][0][0] * cosh(a_[1][0] * beta[i] * x) +
                        _c['g_n'][0][1] * cosh(a_[1][1] * beta[i] * x)
                    ) * sin(beta[i]*y)
                )
                result.append(
                    (
                        _c['g_n'][1][0] * sinh(a_[1][0] * beta[i] * x) +
                        _c['g_n'][1][1] * sinh(a_[1][1] * beta[i] * x)
                    ) * sin(beta[i]*y)
                )
                result.append(
                    (
                        _c['g_n'][2][0] * sinh(a_[1][0] * beta[i] * x) +
                        _c['g_n'][2][1] * sinh(a_[1][1] * beta[i] * x)
                    ) * sin(beta[i]*y) / beta[i]
                )
                result.append(
                    (
                        _c['g_n'][3][0] * cosh(a_[1][0] * beta[i] * x) +
                        _c['g_n'][3][1] * cosh(a_[1][1] * beta[i] * x)
                    ) * sin(beta[i]*y) / beta[i]
                )

        return result
    return coeff_u_eq


def _wrapper_v(H, a_, d, alpha, beta):
    M = len(alpha)

    _c = dict()
    _c['h_n'] = [
        [
            H[6, 6] * (H[1, 1] + H[1, 2] * a_[0][0]**2) / (d[0] * a_[0][0]),
            -H[6, 6] * (H[1, 1] + H[1, 2] * a_[0][1]**2) / (d[0] * a_[0][1]),
        ],
        [
            (
                H[1, 2] * H[6, 6] + H[1, 2]**2 -
                H[1, 1] * H[2, 2] + H[2, 2] * H[6, 6] * a_[0][0]**2
            ) / d[0],
            (
                -H[1, 2] * H[6, 6] - H[1, 2]**2 +
                H[1, 1] * H[2, 2] - H[2, 2] * H[6, 6] * a_[0][1]**2
            ) / d[0]

        ],
        [
            (-H[1, 1] + H[6, 6] * a_[0][0]**2) / (d[0] * a_[0][0]),
            (H[1, 1] - H[6, 6] * a_[0][1]**2) / (d[0] * a_[0][1])
        ],
        [
            (H[1, 2] + H[6, 6]) / d[0],
            -(H[1, 2] + H[6, 6]) / d[0]
        ]
    ]

    _c['g_0'] = {
        1: H[1, 1] * H[6, 6] * (a_[1][0]**2 - a_[1][1]**2) / d[1],
        3: H[1, 1] * (a_[1][0]**2 - a_[1][1]**2) / d[1]
    }

    _c['g_n'] = [
        [
            -H[6, 6] * (H[1, 2] + H[1, 1] * a_[1][0]**2) / (d[1] * a_[1][0]),
            H[6, 6] * (H[1, 2] + H[1, 1] * a_[1][1]**2) / (d[1] * a_[1][1])
        ],
        [
            H[6, 6] * (H[1, 2] + H[1, 1] * a_[1][0]**2) / d[1],
            -H[6, 6] * (H[1, 2] + H[1, 1] * a_[1][1]**2) / d[1]
        ],
        [
            -(H[1, 2] + H[6, 6]) / d[1],
            (H[1, 2] + H[6, 6]) / d[1]
        ],
        [
            (-H[6, 6] + H[1, 1] * a_[1][0]**2) / (d[1] * a_[1][0]),
            (H[6, 6] - H[1, 1] * a_[1][1]**2) / (d[1] * a_[1][1])
        ]
    ]

    def coeff_v_eq(x, y):
        result = [0, 0]

        if x == 0:
            result.extend([0, ] * 4*M + [_c['g_0'][1], 0])
            for i in range(M):
                result.extend(
                    [0, sum(_c['g_n'][1]) * cos(beta[i]*y), 0, 0]
                )
        elif y == 0:
            for i in range(M):
                result.extend(
                    [0, sum(_c['h_n'][1]) * sin(alpha[i]*x), 0, 0]
                )

            result.extend([_c['g_0'][1], _c['g_0'][3]*x])

            for i in range(M):
                result.append(
                    _c['g_n'][0][0] * sinh(a_[1][0] * beta[i] * x) +
                    _c['g_n'][0][1] * sinh(a_[1][1] * beta[i] * x)
                )
                result.append(
                    _c['g_n'][1][0] * cosh(a_[1][0] * beta[i] * x) +
                    _c['g_n'][1][1] * cosh(a_[1][1] * beta[i] * x)
                )
                result.append(
                    (
                        _c['g_n'][2][0] * cosh(a_[1][0] * beta[i] * x) +
                        _c['g_n'][2][1] * cosh(a_[1][1] * beta[i] * x)
                    ) / beta[i]
                )
                result.append(
                    (
                        _c['g_n'][3][0] * sinh(a_[1][0] * beta[i] * x) +
                        _c['g_n'][3][1] * sinh(a_[1][1] * beta[i] * x)
                    ) / beta[i]
                )
        else:
            for i in range(M):
                result.append(
                    (
                        _c['h_n'][0][0] * sinh(a_[0][0] * alpha[i] * y) +
                        _c['h_n'][0][1] * sinh(a_[0][1] * alpha[i] * y)
                    ) * sin(alpha[i]*x)
                )
                result.append(
                    (
                        _c['h_n'][1][0] * cosh(a_[0][0] * alpha[i] * y) +
                        _c['h_n'][1][1] * cosh(a_[0][1] * alpha[i] * y)
                    ) * sin(alpha[i]*x)
                )
                result.append(
                    (
                        _c['h_n'][2][0] * sinh(a_[0][0] * alpha[i] * y) +
                        _c['h_n'][2][1] * sinh(a_[0][1] * alpha[i] * y)
                    ) * sin(alpha[i]*x) / alpha[i]
                )
                result.append(
                    (
                        _c['h_n'][3][0] * cosh(a_[0][0] * alpha[i] * y) +
                        _c['h_n'][3][1] * cosh(a_[0][1] * alpha[i] * y)
                    ) * sin(alpha[i]*x) / alpha[i]
                )

            result.extend([_c['g_0'][1], _c['g_0'][3]*x])

            for i in range(M):
                result.append(
                    (
                        _c['g_n'][0][0] * sinh(a_[1][0] * beta[i] * x) +
                        _c['g_n'][0][1] * sinh(a_[1][1] * beta[i] * x)
                    ) * cos(beta[i]*y)
                )
                result.append(
                    (
                        _c['g_n'][1][0] * cosh(a_[1][0] * beta[i] * x) +
                        _c['g_n'][1][1] * cosh(a_[1][1] * beta[i] * x)
                    ) * cos(beta[i] * y)
                )
                result.append(
                    (
                        _c['g_n'][2][0] * cosh(a_[1][0] * beta[i] * x) +
                        _c['g_n'][2][1] * cosh(a_[1][1] * beta[i] * x)
                    ) * cos(beta[i]*y) / beta[i]
                )
                result.append(
                    (
                        _c['g_n'][3][0] * sinh(a_[1][0] * beta[i] * x) +
                        _c['g_n'][3][1] * sinh(a_[1][1] * beta[i] * x)
                    ) * cos(beta[i]*y) / beta[i]
                )

        return result
    return coeff_v_eq


def _wrapper_sigma_x(H, a_, d, alpha, beta):
    M = len(alpha)

    _c = dict()
    _c['h_n'] = [
        [
            H[6, 6] * (H[1, 2]**2 - H[1, 1] * H[2, 2]) * a_[0][0]**2 / d[0],
            -H[6, 6] * (H[1, 2]**2 - H[1, 1] * H[2, 2]) * a_[0][1]**2 / d[0]
        ],
        [
            (
                H[1, 1] * H[1, 2] * H[6, 6]
                + a_[0][0]**2 * (
                    -H[1, 1] * H[1, 2] * H[2, 2] + H[6, 6] * H[1, 2]**2 +
                    H[1, 2]**3 + H[1, 1] * H[2, 2] * H[6, 6]
                )
                + a_[0][0]**4 * H[1, 2] * H[2, 2] * H[6, 6]
            ) / (d[0] * a_[0][0]),
            (
                - H[1, 1] * H[1, 2] * H[6, 6]
                - a_[0][1]**2 * (
                    -H[1, 1] * H[1, 2] * H[2, 2] + H[6, 6] * H[1, 2]**2 +
                    H[1, 2]**3 + H[1, 1] * H[2, 2] * H[6, 6]
                )
                - a_[0][1]**4 * H[1, 2] * H[2, 2] * H[6, 6]
            ) / (d[0] * a_[0][1])
        ],
        [
            H[6, 6] * (H[1, 1] + H[1, 2] * a_[0][0]**2) / d[0],
            - H[6, 6] * (H[1, 1] + H[1, 2] * a_[0][1]**2) / d[0]
        ],
        [
            (
                H[1, 1] * H[6, 6] + a_[0][0]**2 *
                (-H[1, 1] * H[2, 2] + H[1, 2] * H[6, 6] + H[1, 2]**2)
            ) / (d[0] * a_[0][0]),
            (
                -H[1, 1] * H[6, 6] + a_[0][1]**2 *
                (H[1, 1] * H[2, 2] - H[1, 2] * H[6, 6] - H[1, 2]**2)
            ) / (d[0] * a_[0][1])
        ]
    ]

    _c['g_n'] = [
        [
            (
                H[6, 6] * H[1, 2]**2
                + a_[1][0]**2 * (
                    2 * H[1, 1] * H[1, 2] * H[6, 6] +
                    H[1, 1] * H[1, 2]**2 - H[2, 2] * H[1, 1]**2
                )
                + a_[1][0]**4 * H[6, 6] * H[1, 1]**2
            ) / (d[1] * a_[1][0]),
            (
                - H[6, 6] * H[1, 2]**2
                - a_[1][1]**2 * (
                    2 * H[1, 1] * H[1, 2] * H[6, 6] +
                    H[1, 1] * H[1, 2]**2 - H[2, 2] * H[1, 1]**2
                )
                - a_[1][1]**4 * H[6, 6] * H[1, 1]**2
            ) / (d[1] * a_[1][1])
        ],
        [
            H[6, 6] * (H[1, 1] * H[2, 2] - H[1, 2]**2) / d[1],
            - H[6, 6] * (H[1, 1] * H[2, 2] - H[1, 2]**2) / d[1],
        ],
        [
            (
                - H[1, 1] * H[2, 2] + H[1, 2]**2
                + H[1, 2] * H[6, 6] + H[1, 1] * H[6, 6] * a_[1][0]**2
            ) / d[1],
            (
                H[1, 1] * H[2, 2] - H[1, 2]**2
                - H[1, 2] * H[6, 6] - H[1, 1] * H[6, 6] * a_[1][1]**2
            ) / d[1]
        ],
        [
            H[6, 6] * (H[1, 2] + H[1, 1] * a_[1][0]**2) / (d[1] * a_[1][0]),
            -H[6, 6] * (H[1, 2] + H[1, 1] * a_[1][1]**2) / (d[1] * a_[1][1])
        ]
    ]

    def coeff_sigma_x_eq(x, y):
        result = [0, 0]

        '''
        if x == 0:
            result.extend([0,] * (4*M + 2))

            for i in range(M):
                result.extend(
                    [0, 0, sum(_c['g_n'][2]) * sin(beta[i]*y), 0]
                )
        elif y == 0:
            for i in range(M):
                result.extend(
                    [
                        sum(_c['h_n'][0]) * alpha[i] * sin(alpha[i]*x),
                        0,
                        sum(_c['h_n'][2]) * sin(alpha[i]*x),
                        0
                    ]
                )

            result.extend([0, ] * (4*M+2))
        else:
            for i in range(M):
                result.append(
                    (
                        _c['h_n'][0][0] * cosh(a_[0][0] * alpha[i] * y) +
                        _c['h_n'][0][1] * cosh(a_[0][1] * alpha[i] * y)
                    ) * sin(alpha[i]*x) * alpha[i]
                )
                result.append(
                    (
                        _c['h_n'][1][0] * sinh(a_[0][0] * alpha[i] * y) +
                        _c['h_n'][1][1] * sinh(a_[0][1] * alpha[i] * y)
                    ) * sin(alpha[i]*x) * alpha[i]
                )
                result.append(
                    (
                        _c['h_n'][2][0] * cosh(a_[0][0] * alpha[i] * y) +
                        _c['h_n'][2][1] * cosh(a_[0][1] * alpha[i] * y)
                    ) * sin(alpha[i]*x)
                )
                result.append(
                    (
                        _c['h_n'][3][0] * sinh(a_[0][0] * alpha[i] * y) +
                        _c['h_n'][3][1] * sinh(a_[0][1] * alpha[i] * y)
                    ) * sin(alpha[i]*x)
                )

            result.extend([0, 0])

            for i in range(M):
                result.append(
                    (
                        _c['g_n'][0][0] * sinh(a_[1][0] * beta[i] * x) +
                        _c['g_n'][0][1] * sinh(a_[1][1] * beta[i] * x)
                    ) * sin(beta[i]*y) * beta[i]
                )
                result.append(
                    (
                        _c['g_n'][1][0] * cosh(a_[1][0] * beta[i] * x) +
                        _c['g_n'][1][1] * cosh(a_[1][1] * beta[i] * x)
                    ) * sin(beta[i]*y) * beta[i]
                )
                result.append(
                    (
                        _c['g_n'][2][0] * cosh(a_[1][0] * beta[i] * x) +
                        _c['g_n'][2][1] * cosh(a_[1][1] * beta[i] * x)
                    ) * sin(beta[i]*y)
                )
                result.append(
                    (
                        _c['g_n'][3][0] * sinh(a_[1][0] * beta[i] * x) +
                        _c['g_n'][3][1] * sinh(a_[1][1] * beta[i] * x)
                    ) * sin(beta[i]*y)
                )
        '''
        for i in range(M):
            result.append(
                (
                    _c['h_n'][0][0] * cosh(a_[0][0] * alpha[i] * y) +
                    _c['h_n'][0][1] * cosh(a_[0][1] * alpha[i] * y)
                ) * sin(alpha[i] * x) * alpha[i]
            )
            result.append(
                (
                    _c['h_n'][1][0] * sinh(a_[0][0] * alpha[i] * y) +
                    _c['h_n'][1][1] * sinh(a_[0][1] * alpha[i] * y)
                ) * sin(alpha[i] * x) * alpha[i]
            )
            result.append(
                (
                    _c['h_n'][2][0] * cosh(a_[0][0] * alpha[i] * y) +
                    _c['h_n'][2][1] * cosh(a_[0][1] * alpha[i] * y)
                ) * sin(alpha[i] * x)
            )
            result.append(
                (
                    _c['h_n'][3][0] * sinh(a_[0][0] * alpha[i] * y) +
                    _c['h_n'][3][1] * sinh(a_[0][1] * alpha[i] * y)
                ) * sin(alpha[i] * x)
            )

        result.extend([0, 0])

        for i in range(M):
            result.append(
                (
                    _c['g_n'][0][0] * sinh(a_[1][0] * beta[i] * x) +
                    _c['g_n'][0][1] * sinh(a_[1][1] * beta[i] * x)
                ) * sin(beta[i] * y) * beta[i]
            )
            result.append(
                (
                    _c['g_n'][1][0] * cosh(a_[1][0] * beta[i] * x) +
                    _c['g_n'][1][1] * cosh(a_[1][1] * beta[i] * x)
                ) * sin(beta[i] * y) * beta[i]
            )
            result.append(
                (
                    _c['g_n'][2][0] * cosh(a_[1][0] * beta[i] * x) +
                    _c['g_n'][2][1] * cosh(a_[1][1] * beta[i] * x)
                ) * sin(beta[i] * y)
            )
            result.append(
                (
                    _c['g_n'][3][0] * sinh(a_[1][0] * beta[i] * x) +
                    _c['g_n'][3][1] * sinh(a_[1][1] * beta[i] * x)
                ) * sin(beta[i] * y)
            )

        return result
    return coeff_sigma_x_eq


def _wrapper_sigma_y(H, a_, d, alpha, beta):
    M = len(alpha)

    _c = dict()
    _c['h_n'] = [
        [
            H[6, 6] * (H[1, 1] * H[2, 2] - H[1, 2]**2) / d[0],
            -H[6, 6] * (H[1, 1] * H[2, 2] - H[1, 2]**2) / d[0]
        ],
        [
            (
                H[6, 6] * H[1, 2]**2
                + a_[0][0]**2 * (
                    2 * H[1, 2] * H[2, 2] * H[6, 6]
                    - H[1, 1] * H[2, 2]**2 + H[2, 2] * H[1, 2]**2
                )
                + a_[0][0]**4 * H[6, 6] * H[2, 2]**2
            ) / (d[0] * a_[0][0]),
            (
                - H[6, 6] * H[1, 2]**2
                - a_[0][1]**2 * (
                    2 * H[1, 2] * H[2, 2] * H[6, 6]
                    - H[1, 1] * H[2, 2]**2 + H[2, 2] * H[1, 2]**2
                )
                - a_[0][1]**4 * H[6, 6] * H[2, 2]**2
            ) / (d[0] * a_[0][1]),
        ],
        [
            (
                H[1, 2] * H[6, 6] + H[1, 2]**2
                - H[1, 1] * H[2, 2] + H[2, 2] * H[6, 6] * a_[0][0]**2
            ) / d[0],
            (
                - H[1, 2] * H[6, 6] - H[1, 2]**2
                + H[1, 1] * H[2, 2] - H[2, 2] * H[6, 6] * a_[0][1]**2
            ) / d[0]
        ],
        [
            H[6, 6] * (H[1, 2] + H[2, 2] * a_[0][0]**2) / (d[0] * a_[0][0]),
            - H[6, 6] * (H[1, 2] + H[2, 2] * a_[0][1]**2) / (d[0] * a_[0][1])
        ]
    ]

    _c['g_n'] = [
        [
            (
                H[1, 2] * H[2, 2] * H[6, 6]
                + a_[1][0]**2 * (
                    H[1, 1] * H[2, 2] * H[6, 6] + H[6, 6] * H[1, 2]**2 +
                    H[1, 2]**3 - H[1, 1] * H[1, 2] * H[2, 2]
                )
                + a_[1][0]**4 * H[1, 1] * H[1, 2] * H[6, 6]
            ) / (d[1] * a_[1][0]),
            (
                - H[1, 2] * H[2, 2] * H[6, 6]
                - a_[1][1]**2 * (
                    H[1, 1] * H[2, 2] * H[6, 6] + H[6, 6] * H[1, 2]**2 +
                    H[1, 2]**3 - H[1, 1] * H[1, 2] * H[2, 2]
                )
                - a_[1][1]**4 * H[1, 1] * H[1, 2] * H[6, 6]
            ) / (d[1] * a_[1][1])
        ],
        [
            H[6, 6] * (H[1, 2]**2 - H[1, 1] * H[2, 2]) * a_[1][0]**2 / d[1],
            - H[6, 6] * (H[1, 2]**2 - H[1, 1] * H[2, 2]) * a_[1][1]**2 / d[1]
        ],
        [
            H[6, 6] * (H[2, 2] + H[1, 2] * a_[1][0]**2) / d[1],
            - H[6, 6] * (H[2, 2] + H[1, 2] * a_[1][1]**2) / d[1]
        ],
        [
            (
                H[2, 2] * H[6, 6] + a_[1][0]**2 *
                (-H[1, 1] * H[2, 2] + H[1, 2] * H[6, 6] + H[1, 2]**2)
            ) / (d[1] * a_[1][0]),
            (
                - H[2, 2] * H[6, 6] + a_[1][1]**2 *
                (H[1, 1] * H[2, 2] - H[1, 2] * H[6, 6] - H[1, 2]**2)
            ) / (d[1] * a_[1][1])
        ]
    ]

    def coeff_sigma_y_eq(x, y):
        result = [0, 0]
        '''
        if x == 0:
            result.extend([0, ] * (4*M + 2))

            for i in range(M):
                result.extend(
                    [
                        0,
                        sum(_c['g_n'][1]) * beta[i] * sin(beta[i]*y),
                        sum(_c['g_n'][2]) * sin(beta[i]*y),
                        0
                    ]
                )
        elif y == 0:
            for i in range(M):
                result.extend(
                    [0, 0, sum(_c['h_n'][2]) * sin(alpha[i]*x), 0]
                )

            result.extend([0, ] * (4*M + 2))
        else:
            for i in range(M):
                result.append(
                    (
                        _c['h_n'][0][0] * cosh(a_[0][0] * alpha[i] * y) +
                        _c['h_n'][0][1] * cosh(a_[0][1] * alpha[i] * y)
                    ) * sin(alpha[i]*x) * alpha[i]
                )
                result.append(
                    (
                        _c['h_n'][1][0] * sinh(a_[0][0] * alpha[i] * y) +
                        _c['h_n'][1][1] * sinh(a_[0][1] * alpha[i] * y)
                    ) * sin(alpha[i]*x) * alpha[i]
                )
                result.append(
                    (
                        _c['h_n'][2][0] * cosh(a_[0][0] * alpha[i] * y) +
                        _c['h_n'][2][1] * cosh(a_[0][1] * alpha[i] * y)
                    ) * sin(alpha[i]*x)
                )
                result.append(
                    (
                        _c['h_n'][3][0] * sinh(a_[0][0] * alpha[i] * y) +
                        _c['h_n'][3][1] * sinh(a_[0][1] * alpha[i] * y)
                    ) * sin(alpha[i]*x)
                )

            result.extend([0, 0])

            for i in range(M):
                result.append(
                    (
                        _c['g_n'][0][0] * sinh(a_[1][0] * beta[i] * x) +
                        _c['g_n'][0][1] * sinh(a_[1][1] * beta[i] * x)
                    ) * sin(beta[i]*y) * beta[i]
                )
                result.append(
                    (
                        _c['g_n'][1][0] * cosh(a_[1][0] * beta[i] * x) +
                        _c['g_n'][1][1] * cosh(a_[1][1] * beta[i] * x)
                    ) * sin(beta[i]*y) * beta[i]
                )
                result.append(
                    (
                        _c['g_n'][2][0] * cosh(a_[1][0] * beta[i] * x) +
                        _c['g_n'][2][1] * cosh(a_[1][1] * beta[i] * x)
                    ) * sin(beta[i]*y)
                )
                result.append(
                    (
                        _c['g_n'][3][0] * sinh(a_[1][0] * beta[i] * x) +
                        _c['g_n'][3][1] * sinh(a_[1][1] * beta[i] * x)
                    ) * sin(beta[i]*y)
                )
        '''
        for i in range(M):
            result.append(
                (
                    _c['h_n'][0][0] * cosh(a_[0][0] * alpha[i] * y) +
                    _c['h_n'][0][1] * cosh(a_[0][1] * alpha[i] * y)
                ) * sin(alpha[i] * x) * alpha[i]
            )
            result.append(
                (
                    _c['h_n'][1][0] * sinh(a_[0][0] * alpha[i] * y) +
                    _c['h_n'][1][1] * sinh(a_[0][1] * alpha[i] * y)
                ) * sin(alpha[i] * x) * alpha[i]
            )
            result.append(
                (
                    _c['h_n'][2][0] * cosh(a_[0][0] * alpha[i] * y) +
                    _c['h_n'][2][1] * cosh(a_[0][1] * alpha[i] * y)
                ) * sin(alpha[i] * x)
            )
            result.append(
                (
                    _c['h_n'][3][0] * sinh(a_[0][0] * alpha[i] * y) +
                    _c['h_n'][3][1] * sinh(a_[0][1] * alpha[i] * y)
                ) * sin(alpha[i] * x)
            )

        result.extend([0, 0])

        for i in range(M):
            result.append(
                (
                    _c['g_n'][0][0] * sinh(a_[1][0] * beta[i] * x) +
                    _c['g_n'][0][1] * sinh(a_[1][1] * beta[i] * x)
                ) * sin(beta[i] * y) * beta[i]
            )
            result.append(
                (
                    _c['g_n'][1][0] * cosh(a_[1][0] * beta[i] * x) +
                    _c['g_n'][1][1] * cosh(a_[1][1] * beta[i] * x)
                ) * sin(beta[i] * y) * beta[i]
            )
            result.append(
                (
                    _c['g_n'][2][0] * cosh(a_[1][0] * beta[i] * x) +
                    _c['g_n'][2][1] * cosh(a_[1][1] * beta[i] * x)
                ) * sin(beta[i] * y)
            )
            result.append(
                (
                    _c['g_n'][3][0] * sinh(a_[1][0] * beta[i] * x) +
                    _c['g_n'][3][1] * sinh(a_[1][1] * beta[i] * x)
                ) * sin(beta[i] * y)
            )

        return result

    return coeff_sigma_y_eq


def _wrapper_tau(H, a_, d, alpha, beta):
    M = len(alpha)

    _c = dict()
    _c['h_0'] = {
        3: H[2, 2] * H[6, 6] * (a_[0][0]**2 - a_[0][1]**2) / d[0]
    }

    _c['h_n'] = [
        [
            H[6, 6]**2 * (
                H[1, 1] + 2 * H[1, 2] * a_[0][0]**2 + H[2, 2] * a_[0][0]**4
            ) / (d[0] * a_[0][0]),
            - H[6, 6]**2 * (
                H[1, 1] + 2 * H[1, 2] * a_[0][1]**2 + H[2, 2] * a_[0][1]**4
            ) / (d[0] * a_[0][1])
        ],
        [
            H[6, 6] * (H[1, 2]**2 - H[1, 1]*H[2, 2]) / d[0],
            - H[6, 6] * (H[1, 2]**2 - H[1, 1]*H[2, 2]) / d[0]
        ],
        [
            - H[6, 6] * (H[1, 1] + H[1, 2]*a_[0][0]**2) / (d[0] * a_[0][0]),
            H[6, 6] * (H[1, 1] + H[1, 2]*a_[0][1]**2) / (d[0] * a_[0][1])
        ],
        [
            H[6, 6] * (H[1, 2] + H[2, 2] * a_[0][0]**2) / d[0],
            - H[6, 6] * (H[1, 2] + H[2, 2] * a_[0][1]**2) / d[0]
        ]
    ]

    _c['g_0'] = {
        3: H[1, 1] * H[6, 6] * (a_[1][0]**2 - a_[1][1]**2) / d[1]
    }

    _c['g_n'] = [
        [
            H[6, 6] * (H[1, 2]**2 - H[1, 1] * H[2, 2]) / d[1],
            - H[6, 6] * (H[1, 2]**2 - H[1, 1] * H[2, 2]) / d[1]
        ],
        [
            H[6, 6]**2 * (
                H[2, 2] + 2 * H[1, 2] * a_[1][0]**2 + H[1, 1] * a_[1][0]**4
            ) / (d[1] * a_[1][0]),
            - H[6, 6]**2 * (
                H[2, 2] + 2 * H[1, 2] * a_[1][1]**2 + H[1, 1] * a_[1][1]**4
            ) / (d[1] * a_[1][1])
        ],
        [
            - H[6, 6] * (H[2, 2] + H[1, 2] * a_[1][0]**2) / (d[1] * a_[1][0]),
            H[6, 6] * (H[2, 2] + H[1, 2] * a_[1][1]**2) / (d[1] * a_[1][1])
        ],
        [
            H[6, 6] * (H[1, 2] + H[1, 1] * a_[1][0]**2) / d[1],
            - H[6, 6] * (H[1, 2] + H[1, 1] * a_[1][1]**2) / d[1]
        ]
    ]

    def coeff_tau_eq(x, y):
        result = [0, _c['h_0'][3]]

        for i in range(M):
            result.append(
                (
                    _c['h_n'][0][0] * sinh(a_[0][0] * alpha[i] * y) +
                    _c['h_n'][0][1] * sinh(a_[0][1] * alpha[i] * y)
                ) * cos(alpha[i] * x) * alpha[i]
            )
            result.append(
                (
                    _c['h_n'][1][0] * cosh(a_[0][0] * alpha[i] * y) +
                    _c['h_n'][1][1] * cosh(a_[0][1] * alpha[i] * y)
                ) * cos(alpha[i] * x) * alpha[i]
            )
            result.append(
                (
                    _c['h_n'][2][0] * sinh(a_[0][0] * alpha[i] * y) +
                    _c['h_n'][2][1] * sinh(a_[0][1] * alpha[i] * y)
                ) * cos(alpha[i] * x)
            )
            result.append(
                (
                    _c['h_n'][3][0] * cosh(a_[0][0] * alpha[i] * y) +
                    _c['h_n'][3][1] * cosh(a_[0][1] * alpha[i] * y)
                ) * cos(alpha[i] * x)
            )

        result.extend([0, _c['g_0'][3]])

        for i in range(M):
            result.append(
                (
                    _c['g_n'][0][0] * cosh(a_[1][0] * beta[i] * x) +
                    _c['g_n'][0][1] * cosh(a_[1][1] * beta[i] * x)
                ) * cos(beta[i] * y) * beta[i]
            )
            result.append(
                (
                    _c['g_n'][1][0] * sinh(a_[1][0] * beta[i] * x) +
                    _c['g_n'][1][1] * sinh(a_[1][1] * beta[i] * x)
                ) * cos(beta[i] * y) * beta[i]
            )
            result.append(
                (
                    _c['g_n'][2][0] * sinh(a_[1][0] * beta[i] * x) +
                    _c['g_n'][2][1] * sinh(a_[1][1] * beta[i] * x)
                ) * cos(beta[i] * y)
            )
            result.append(
                (
                    _c['g_n'][3][0] * cosh(a_[1][0] * beta[i] * x) +
                    _c['g_n'][3][1] * cosh(a_[1][1] * beta[i] * x)
                ) * cos(beta[i] * y)
            )

        return result
    return coeff_tau_eq


def get_coeff_function(H, a_, d, alpha, beta, func_name):
    return globals()['_wrapper_' + func_name](H, a_, d, alpha, beta)


# TODO: ? don't use M and range(M) in coeff functions
# TODO: ? refactor
if __name__ == '__main__':
    H = {(i, j): i + 0.5 * j for i in (1, 2, 6) for j in (1, 2, 6)}
    M = 30
    a = [[1, 2], [3, 4]]
    d = [5, 6]
    alpha = [i * 0.1 + 3 for i in range(1, M + 1)]
    beta = [i + 2 for i in range(1, M + 1)]

    from timeit import timeit

    c_func = get_coeff_function(H, a, d, alpha, beta, 'tau')
    print(timeit(
        'c_func(3, 2)',
        number=100,
        globals=globals()
    ))

    print('-'*15)
    print('Expected len:', M*8 + 4)
    print('when x = 0:', len(c_func(0, 2)))
    print('when y = 0:', len(c_func(2, 0)))
    print('when x, y != 0:', len(c_func(3, 2)))

"""