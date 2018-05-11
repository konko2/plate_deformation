from math_utilites import sin, cos, sinh, cosh


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
