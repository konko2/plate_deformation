from tkinter import *
import re
import sympy as sp
from constants import Borders, ConditionTypes
from abc import ABC, abstractmethod


NUMBER_PATTERN = r'(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?'
UNICODE_SIGMA = '\u03C3'
UNICODE_TAU = '\u03C4'
UNICODE_NU = '\u03BD'


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


def convert_str_nums_to_mpf(string):
    matches = list(re.finditer(NUMBER_PATTERN, string))
    for match in reversed(matches):
        string = string[:match.start()] + f"mpf('{match[0]}')" + string[match.end():]

    return string


class ValidatingEntry(Entry, ABC):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.default_bg = self['bg']

        self.bind('<Return>', self.check)
        self.bind('<FocusOut>', self.check)

    @abstractmethod
    def check_func(self):
        pass

    def check(self, *_):
        if not self.get():
            self.config(bg=self.default_bg)
            return False

        if self.check_func():
            self.config(bg=self.default_bg)
            return True
        else:
            self.config(bg='red')
            return False


class IntEntry(ValidatingEntry):
    def check_func(self):
        try:
            int(self.get())
            return True
        except Exception:
            return False


class NaturalsEntry(ValidatingEntry):
    def check_func(self):
        try:
            if int(self.get()) > 0:
                return True
        except Exception:
            pass
        return False


class PositiveNumEntry(ValidatingEntry):
    def check_func(self):
        try:
            if sp.N(self.get()) > 0:
                return True
        except Exception:
            pass
        return False


class RealEntry(ValidatingEntry):
    def check_func(self):
        try:
            if sp.N(self.get()).is_real:
                return True
        except Exception:
            pass
        return False


class FuncEntry(ValidatingEntry):
    def __init__(self, *args, func_args, **kwargs):
        super().__init__(*args, **kwargs)

        self.func_args_names = func_args

    def check_func(self):
        try:
            func_args = set(sp.symbols(self.func_args_names))
            if not (sp.N(self.get()).free_symbols - func_args):
                return True
        except Exception:
            pass
        return False
