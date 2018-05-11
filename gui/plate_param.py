from tkinter import *
import sys
from math import floor
import sympy as sp
from constants import MaterialTypes, Borders, ConditionTypes
from itertools import cycle
from .utils import UNICODE_TAU, UNICODE_SIGMA, UNICODE_NU, cond_funcs, border_func_args, ValidatingEntry, \
    IntEntry, NaturalsEntry, PositiveNumEntry, RealEntry, FuncEntry


class MaterialTypeFrame(Frame):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.var = Variable(self, name='problem_type')

        for m_type in MaterialTypes:
            Radiobutton(
                self,
                text=m_type.name.title(),
                variable=self.var,
                value=m_type.value
            ).pack(side=LEFT, expand=True)

    def get(self):
        return next(m_type for m_type in MaterialTypes if m_type.value == self.var.get())


class PlateSizeFrame(Frame):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.canvas = Canvas(self, width=260, height=160)
        self.canvas.grid(row=0, column=0, columnspan=2)

        self.scheme_color = next(self.colors)
        self.create_scheme()

        Label(self, text="Plate size:").grid(row=1, column=0, columnspan=2)

        self.entries = dict()
        for row, name in ((2, 'h'), (3, 'a')):
            Label(self, text=(name + " =")).grid(row=row, column=0, sticky='e')
            entry = PositiveNumEntry(self)
            entry.grid(row=row, column=1, sticky='w')

            entry.bind('<Return>', lambda _: self.update_scheme(), add='+')
            entry.bind('<FocusOut>', lambda _: self.update_scheme(), add='+')

            self.entries[name] = entry

    def create_scheme(self, h=100, a=200):
        canvas = self.canvas

        origin_c = (30, 30)
        rectangle_c = origin_c + (origin_c[0] + a, origin_c[1] + h)
        x_axe_c = (origin_c[0], 10, origin_c[0], rectangle_c[3] + 30)
        y_axe_c = (10, origin_c[1], rectangle_c[2] + 30, origin_c[1])

        rectangle = canvas.create_rectangle(
            *rectangle_c,
            width=0,
            fill=self.scheme_color
        )
        canvas.create_line(*x_axe_c, arrow='last')
        canvas.create_line(*y_axe_c, arrow='last')

        canvas.create_text(*origin_c, anchor='se', text='0')
        canvas.create_text(x_axe_c[2] - 5, x_axe_c[3], anchor='se', text='x')
        canvas.create_text(y_axe_c[2], y_axe_c[3] - 5, anchor='se', text='y')
        canvas.create_text(x_axe_c[0], rectangle_c[3], anchor='e', text='h')
        canvas.create_text(rectangle_c[2], y_axe_c[1], anchor='s', text='a')

        canvas.tag_bind(
            rectangle,
            '<Button-1>',
            lambda _: canvas.itemconfig(rectangle, fill=self.get_next_color())
        )

    def update_scheme(self):
        try:
            # TODO: try to union with EntryExpr check
            h = sp.N(self.entries['h'].get())
            a = sp.N(self.entries['a'].get())

            if h > 0 and a > 0:
                koeff = min([100 / h, 200 / a])
                pix_h = max([2, floor(koeff * h)])
                pix_a = max([2, floor(koeff * a)])

                self.canvas.delete('all')
                self.create_scheme(pix_h, pix_a)
        except Exception:
            pass

    def get_next_color(self):
        self.scheme_color = next(self.colors)
        return self.scheme_color

    def get(self):
        return {name: self.entries[name].get() for name in self.entries}

    colors = cycle(['red', 'green', 'blue', 'cyan', 'yellow', 'magenta'])


class CharacteristicsFrame(Frame):
    def __init__(self, *args, var_names, **kwargs):
        super().__init__(*args, **kwargs)

        self.var_names = var_names

        Label(self, text="Material characteristics:").grid(row=0, columnspan=2)

        for row, var_name in enumerate(var_names, start=1):
            Label(self, text=(var_name + ' =')).grid(row=row, column=0, sticky='e')
            RealEntry(
                self,
                textvariable=Variable(self, name=var_name),
            ).grid(row=row, column=1, sticky='w')

    def get(self):
        return {name: self.getvar(name) for name in self.var_names}


class BorderConditionsFrame(Frame):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.spinboxes = dict()
        self.entries = dict()
        self.func_labels = dict()

        Label(self, text="Values of selected functions for each border:").pack()

        for border in Borders:
            frame = Frame(self)

            Label(
                frame,
                text=f'Condition type for {border.name.lower()} border:'
            ).pack(side=LEFT)

            spinbox = Spinbox(
                frame,
                values=[c.name.lower() for c in ConditionTypes],
                state='readonly',
                width=10
            )
            self.spinboxes[border] = spinbox
            spinbox.pack(side=LEFT)

            frame.pack()
            self.create_funcs_frame(border).pack(fill=X)

            def spinbox_command(border=border):
                color = 'black' if self.check_conditions() else 'red'
                for spinbox in self.spinboxes.values():
                    spinbox.config(fg=color)

                self.change_labels_text(border)

            spinbox.config(command=spinbox_command)

    def create_funcs_frame(self, border):
        frame = Frame(self)

        self.func_labels[border] = list()
        self.entries[border] = list()

        for row in range(2):
            label = Label(frame)
            entry = FuncEntry(frame, func_args=border_func_args[border])

            self.func_labels[border].append(label)
            self.entries[border].append(entry)

            label.grid(row=row, column=0, sticky='e')
            entry.grid(row=row, column=1, sticky='ew')
        self.change_labels_text(border)

        frame.grid_columnconfigure(1, weight=1)

        return frame

    def change_labels_text(self, border):
        condition_name = self.spinboxes[border].get()
        condition = ConditionTypes[condition_name.upper()]

        f_names = cond_funcs[condition][border]
        args = border_func_args[border]

        labels = self.func_labels[border]
        for label, f_name in zip(labels, f_names):
            label.config(text=f'{f_name}({args[0]}, {args[1]}) =')

    # TODO: fix
    def check_conditions(self):
        conditions = [ConditionTypes[spinbox.get().upper()] for spinbox in self.spinboxes.values()]
        if ConditionTypes.CINEMATIC in conditions:
            return True
        if ConditionTypes.MIXED_1 in conditions and ConditionTypes.MIXED_2 in conditions:
            return True
        return False

    def get(self):
        if not self.check_conditions():
            raise Exception

        result = dict()
        for border in Borders:
            spinbox_value = self.spinboxes[border].get()
            condition = ConditionTypes[spinbox_value.upper()]
            func_names = cond_funcs[condition][border]

            # TODO: check entries, raise exc
            entry_values = [entry.get() for entry in self.entries[border]]
            result[border] = {name: str_func for name, str_func in zip(func_names, entry_values)}

        return result


# TODO: refactor
class DigitsFrame(Frame):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        Label(self, text='Digits: ').grid(row=0, column=0, sticky='e')

        self.entry_digit_num = NaturalsEntry(self, text=64, width=5)
        self.entry_digit_num.grid(row=0, column=1, sticky='w')

        Label(self, text='Number of sum terms: ').grid(row=0, column=2, sticky='e')

        self.entry_terms_num = NaturalsEntry(self, text=10, width=5)
        self.entry_terms_num.grid(row=0, column=3, sticky='w')

    def get(self):
        return {
            'digits': self.entry_digit_num.get(),
            'M': self.entry_terms_num.get()
        }


# TODO: refactor
class PreferenceWin(Toplevel):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.computation_param_frame = DigitsFrame(self)

        self.conditions_frame = BorderConditionsFrame(self)
        self.size_frame = PlateSizeFrame(self)

        self.start_button = Button(self, text='Start!')
        self.type_frame = MaterialTypeFrame(self)

        names = ['E_x', 'E_y', 'G_xy', f'{UNICODE_NU}_xy']
        self.orthotropy_frame = CharacteristicsFrame(self, var_names=names)

        names = ['E', f'{UNICODE_NU}_xy']
        self.isotropy_frame = CharacteristicsFrame(self, var_names=names)

        self.type_frame.var.trace_variable('w', self.change_frame_of_characteristic)

        self.type_frame.pack(fill=X)

        self.start_button.pack(side=BOTTOM)
        self.computation_param_frame.pack(side=BOTTOM)

        self.conditions_frame.pack(side=LEFT, fill=X, expand=True, anchor='n')
        self.size_frame.pack()

        self.type_frame.var.set(MaterialTypes.ORTHOTROPY.value)

        self.resizable(width=True, height=False)

    def change_frame_of_characteristic(self, *_):
        problem_type = self.type_frame.get()
        if problem_type == MaterialTypes.ISOTROPY:
            self.orthotropy_frame.forget()
            self.isotropy_frame.pack()
        elif problem_type == MaterialTypes.ORTHOTROPY:
            self.isotropy_frame.forget()
            self.orthotropy_frame.pack()

    # TODO
    def get(self):
        result = {'conditions': self.conditions_frame.get()}
        characteristic_frame = self.isotropy_frame \
            if self.type_frame.get() == MaterialTypes.ISOTROPY else self.orthotropy_frame
        param_frames = [self.computation_param_frame, characteristic_frame, self.size_frame]
        for frame in param_frames:
            result.update(frame.get())

        return result


if __name__ == '__main__':
    tk = Tk()
    tk.iconify()

    win = PreferenceWin(tk)

    win.start_button.config(command=lambda *_: print(win.get()))

    mainloop()
