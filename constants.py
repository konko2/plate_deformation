from enum import Enum


UNICODE_SIGMA = '\u03C3'
UNICODE_TAU = '\u03C4'
UNICODE_NU = '\u03BD'


class Borders(Enum):
    UPPER = 1
    RIGHT = 2
    BOTTOM = 3
    LEFT = 4

    def is_horizontal(self):
        return self in (self.UPPER, self.BOTTOM)

    def is_vertical(self):
        return self in (self.LEFT, self.RIGHT)


class MaterialTypes(Enum):
    ORTHOTROPY = 1
    ISOTROPY = 2


class ConditionTypes(Enum):
    CINEMATIC = 1
    FORCE = 2
    MIXED_1 = 3
    MIXED_2 = 4


border_func_args = {
    Borders.UPPER: ('0', 'y'),
    Borders.RIGHT: ('x', 'a'),
    Borders.BOTTOM: ('h', 'y'),
    Borders.LEFT: ('x', '0')
}

cond_funcs = {
    ConditionTypes.CINEMATIC: {border: ('u', 'v') for border in Borders},
    ConditionTypes.FORCE: {
        border: (
            UNICODE_SIGMA + '_' + ('x' if border.is_horizontal() else 'y'),
            UNICODE_TAU + '_xy'
        )
        for border in Borders
    },
    ConditionTypes.MIXED_1: {
        border: (f'{UNICODE_SIGMA}_x', 'v') if border.is_horizontal() else (f'{UNICODE_SIGMA}_y', 'u')
        for border in Borders
    },
    ConditionTypes.MIXED_2: {
        border: (
            'u' if border.is_horizontal() else 'v',
            UNICODE_TAU + '_xy'
        )
        for border in Borders
    }
}


if __name__ == '__main__':
    print(UNICODE_NU, UNICODE_SIGMA, UNICODE_TAU)
    print(Borders.UPPER.is_horizontal(), Borders.UPPER.is_vertical())
