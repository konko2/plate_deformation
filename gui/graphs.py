from tkinter import *
from math_utilites import mp, lambdify, linspace
from constants import UNICODE_SIGMA, UNICODE_NU, UNICODE_TAU, Borders
from orthotropy import get_solution
from orthotropy.method import find_dots, Timer
from constants import Borders, border_func_args
import matplotlib.pyplot as plt
import matplotlib.backends.backend_tkagg as tkagg
from math_utilites import linspace
from copy import deepcopy
import mpmath
from .utils import convert_str_nums_to_mpf, UNICODE_SIGMA, UNICODE_NU, UNICODE_TAU, PositiveNumEntry


def create_lines_labels(var_name, sections_num):
    fixed_var_name = 'x' if var_name == 'y' else 'y'
    fixed_var_endpoint_name = 'h' if fixed_var_name == 'x' else 'a'

    fixed_var_vals_str = ['0', ]
    fixed_var_vals_str += [f'{i} / {sections_num + 1} * {fixed_var_endpoint_name}' for i in range(1, sections_num + 1)]
    fixed_var_vals_str.append(f'{fixed_var_endpoint_name}')

    labels = [f'{fixed_var_name} = ' + val_str for val_str in fixed_var_vals_str]
    return labels


# TODO: entry validation
def create_graph_frame(master, *args, func_name, var_name, func, h, a, **kwargs):
    root = Frame(master, *args, **kwargs)

    figure = plt.Figure()
    canvas = tkagg.FigureCanvasTkAgg(master=root, figure=figure)
    axes = figure.add_subplot(1, 1, 1)
    axes.set_xlabel(var_name)
    axes.set_ylabel(func_name)

    toolbar_frame = Frame(root)
    tkagg.NavigationToolbar2TkAgg(canvas=canvas, window=toolbar_frame)

    sections_entry = PositiveNumEntry(root)
    plot_button = Button(root, text='plot')

    canvas.get_tk_widget().pack(fill='both', expand=True)
    toolbar_frame.pack(fill='x')

    plot_button.pack(side='bottom')

    Label(root, text='Number of sections:').pack(side='left', expand=True, anchor='e')
    sections_entry.pack(side='right', expand=True, anchor='w')

    if var_name == 'x':
        endpoint, fixed_var_endpoint = h, a
        fixed_var_name = 'y'
    else:
        endpoint, fixed_var_endpoint = a, h
        fixed_var_name = 'x'

    timer = Timer()

    plot_dots = linspace(0, endpoint, 200)

    # TODO:
    # 1) don't use fixed_var_name
    # 2) set_xlabel and ylabel: maybe you can remember it out of this function?
    # 3) last_used_values: how to make it simpler?
    last_used_values = dict()
    def butt_command():
        axes.cla()
        axes.set_xlabel(var_name)
        axes.set_ylabel(func_name)

        print("Plotting graph...")
        timer.start()

        sections_num = int(sections_entry.get())
        fixed_var_vals = linspace(0, fixed_var_endpoint, sections_num + 2)
        labels = create_lines_labels(var_name, sections_num)
        values = dict()
        for label, fixed_var in zip(labels, fixed_var_vals):
            if fixed_var in last_used_values:
                values[fixed_var] = last_used_values[fixed_var]
            else:
                values[fixed_var] = [func(**{var_name: var, fixed_var_name: fixed_var}) for var in plot_dots]
            axes.plot(plot_dots, values[fixed_var], label=label)

        axes.legend(loc='best').draggable()

        canvas.draw()
        print(f'(done in {timer})')

        last_used_values.clear()
        last_used_values.update(values)

    plot_button.config(command=butt_command)

    return root


def create_graph_win(*args, data_from_pr_win, **kwargs):
    str_initial_data = deepcopy(data_from_pr_win)

    mp.dps = int(data_from_pr_win.pop('digits'))
    M = int(data_from_pr_win.pop('M'))

    conditions = data_from_pr_win.pop('conditions')

    data_from_pr_win = {
        k.replace(UNICODE_NU, 'nu'): lambdify([], convert_str_nums_to_mpf(v), modules=mpmath)()
        for k, v in data_from_pr_win.items()
    }

    for border in Borders:
        conditions[border] = {
            (
                name.replace(UNICODE_SIGMA, 'sigma').replace(f'{UNICODE_TAU}_xy', 'tau')
            ): lambdify(
                ['x' if border.is_horizontal() else 'y'],
                convert_str_nums_to_mpf(str_func),
                modules=mpmath
            )
            for name, str_func in conditions[border].items()
        }

    solution = get_solution(report=True, conditions=conditions, M=M, **data_from_pr_win)

    win = Toplevel(*args, **kwargs)
    label_text = 'Border conditions:\n'
    str_conditions = str_initial_data.pop('conditions')
    for border in Borders:
        args = '{}, {}'.format(*border_func_args[border])
        for func_name in str_conditions[border]:
            label_text += f'{func_name}({args}) = {str_conditions[border][func_name]}\n'

    label_text += '\n'.join([k + ' = ' + str(v) for k, v in str_initial_data.items()])

    left_frame = Frame(win)
    Label(left_frame, text=label_text).pack(side='top')

    visible_solution_index = IntVar()
    var_of_func = Variable()

    func_names_and_indicies = [
        ('u', 0),
        ('v', 1),
        ('sigma_x', 2),
        ('sigma_y', 3),
        ('tau', 4)
    ]

    choose_graph_frame = Frame(left_frame)

    func_buttons_frame = Frame(choose_graph_frame)
    graph_frames = dict()
    for func_name, index in func_names_and_indicies:
        Radiobutton(
            func_buttons_frame,
            text=func_name,
            variable=visible_solution_index,
            value=index
        ).pack(anchor='w')
        for var_name in ('x' ,'y'):
            graph_frames[func_name, var_name] = create_graph_frame(
                win,
                func_name=func_name,
                var_name=var_name,
                func=solution[index],
                h=data_from_pr_win['h'],
                a=data_from_pr_win['a']
            )

    var_buttons_frame = Frame(choose_graph_frame)
    for var_name in ('x', 'y'):
        Radiobutton(
            var_buttons_frame,
            text=var_name,
            variable=var_of_func,
            value=var_name
        ).pack(anchor='w')


    def show_choosed_graph(*_):
        for graph_frame in graph_frames.values():
            graph_frame.forget()
        func_name = [couple for couple in func_names_and_indicies if couple[1] == visible_solution_index.get()][0][0]
        graph_frames[func_name, var_of_func.get()].pack(side='right', expand=True,fill='both')

    visible_solution_index.trace('w', show_choosed_graph)
    var_of_func.trace('w', show_choosed_graph)

    #visible_solution_index.set(0)
    #var_of_func.set('x')

    func_buttons_frame.pack(side='left')
    var_buttons_frame.pack(side='right')

    choose_graph_frame.pack(side='bottom')
    left_frame.pack(side='left')
