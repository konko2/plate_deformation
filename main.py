from gui import PreferenceWin, create_graph_win
from tkinter import Tk

tk = Tk()
tk.iconify()

pr_win = PreferenceWin(tk)
pr_win.start_button.config(command=lambda *_: create_graph_win(tk, data_from_pr_win=pr_win.get()))

tk.mainloop()
