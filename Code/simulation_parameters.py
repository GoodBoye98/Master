import tkinter as tk

def rgb_to_hex(r, g, b):
  return ('#{:X}{:X}{:X}').format(r, g, b)


init_vals = {
    "a_0": "0.6",
    "a_1": "0.38",
    "K_w": "0.38",
    "K_l": "1.89",
    "A": "192.2",
    "B": "3.85",
    "C": "13.2",
    "L": "5.0",
    "N": "1024",
    "s": f"{1/40:.3f}",
    "time": "0.0 3.33; 6.67 10.0"}
var_index = {
    "a_0": 0,
    "a_1": 1,
    "K_w": 2,
    "K_l": 3,
    "A": 4,
    "B": 5,
    "C": 6,
    "L": 7,
    "N": 8,
    "s": 9,
    "time": 10}


with open("Simulation Parameters/simulation.cfg", "r") as file:
    lines = [line.rstrip() for line in file]

for line in lines:
    line = line.split(' ', 1)
    if line[0] == "N":
        init_vals["N"] = line[1]
    elif line[0] == "s":
        init_vals["s"] = line[1]
    elif line[0] == "time":
        init_vals["time"] = line[1]


good = "✔"
good_c = "green"
bad = "✘"
bad_c = "red"

bM = rgb_to_hex(20, 20, 20)
fM = rgb_to_hex(220, 220, 220)

master = tk.Tk()
master.configure(bg=bM)

entries = []
labels = []
checks = []

i = 0
for key, value in init_vals.items():
    labels.append(tk.Label(master, text=key, fg=fM, bg=bM))
    labels[-1].grid(row=i)
    checks.append(tk.Label(master, text="✔", fg="green", bg=bM))
    checks[-1].grid(row=i, column=2)
    entries.append(tk.Entry(master, fg=fM, bg=bM))
    entries[-1].grid(row=i, column=1, padx=10, pady=3)
    entries[-1].insert(0, value)
    i += 1

r = len(init_vals)

spectral_check = tk.Label(master, text="✔", fg="green", bg=bM)
spectral_check.grid(row=r, column=2)
finitediff_check = tk.Label(master, text="✔", fg="green", bg=bM)
finitediff_check.grid(row=r+1, column=2)
midpoint_check = tk.Label(master, text="✔", fg="green", bg=bM)
midpoint_check.grid(row=r+2, column=2)

tk.Label(master, text="Method", fg=fM, bg=bM).grid(row=r)
spectral_var = tk.IntVar(value=0)
spectral = tk.Checkbutton(master, text="spectral", variable=spectral_var, fg=fM, bg=bM, selectcolor=bM)
spectral.grid(row=r, column=1)
finitediff_var = tk.IntVar(value=0)
finitediff = tk.Checkbutton(master, text="finite diff", variable=finitediff_var, fg=fM, bg=bM, selectcolor=bM)
finitediff.grid(row=r + 1, column=1)
midpoint_var = tk.IntVar(value=0)
midpoint = tk.Checkbutton(master, text="midpoint", variable=midpoint_var, fg=fM, bg=bM, selectcolor=bM)
midpoint.grid(row=r + 2, column=1)

for line in lines:
    line = line.split(' ', 1)
    if line[0] == "run":
        if "spectral" in line[1]:
            spectral_var.set(1)
        if "finitediff" in line[1]:
            finitediff_var.set(1)
        if "midpoint" in line[1]:
            midpoint_var.set(1)



def writeFile():
    with open("Simulation Parameters/simulation.cfg", "w") as file:
        for i in range(len(entries)):
            file.write(f"{labels[i]['text']} {entries[i].get()}\n")
        str = "spectral " if spectral_var.get() else ""
        str += "finitediff " if finitediff_var.get() else ""
        str += "midpoint " if midpoint_var.get() else ""
        if str:
            file.write(f"run {str}")
        else:
            spectral_var.set(1)
            file.write(f"run spectral")
        file.truncate()

def reset():
    i = 0
    for key, value in init_vals.items():
        entries[i].delete(0, tk.END)
        entries[i].insert(0, value)
        i += 1

def check():

    with open("Simulation Parameters/simulation.cfg", "r") as file:
        lines = [line.rstrip() for line in file]

    for line in lines:
        line = line.split(' ', 1)
        if line[0] not in init_vals:
            if line[0] == "run":
                if "finitediff" in line[1]:
                    if finitediff_var.get():
                        finitediff_check.config(text=good, fg=good_c)
                    else:
                        finitediff_check.config(text=bad, fg=bad_c)
                else:
                    if finitediff_var.get():
                        finitediff_check.config(text=bad, fg=bad_c)
                    else:
                        finitediff_check.config(text=good, fg=good_c)

                if "spectral" in line[1]:
                    if spectral_var.get():
                        spectral_check.config(text=good, fg=good_c)
                    else:
                        spectral_check.config(text=bad, fg=bad_c)
                else:
                    if spectral_var.get():
                        spectral_check.config(text=bad, fg=bad_c)
                    else:
                        spectral_check.config(text=good, fg=good_c)

                if "midpoint" in line[1]:
                    if midpoint_var.get():
                        midpoint_check.config(text=good, fg=good_c)
                    else:
                        midpoint_check.config(text=bad, fg=bad_c)
                else:
                    if midpoint_var.get():
                        midpoint_check.config(text=bad, fg=bad_c)
                    else:
                        midpoint_check.config(text=good, fg=good_c)

        elif entries[var_index[line[0]]].get() != line[1]:
            checks[var_index[line[0]]].config(text=bad, fg=bad_c)
        else:
            checks[var_index[line[0]]].config(text=good, fg=good_c)

    master.after(200, check)

r = len(init_vals) + 3
tk.Button(master, text='Quit', command=master.quit, bg=bM, fg=fM).grid(row=r, column=0, sticky=tk.W, pady=4, padx=10)
tk.Button(master, text='Update', command=writeFile, bg=bM, fg=fM).grid(row=r, column=1, sticky=tk.W, pady=4, padx=10)
tk.Button(master, text='Reset', command=reset, bg=bM, fg=fM).grid(row=r, column=2, sticky=tk.W, pady=4, padx=10)

master.after(0, check)
master.mainloop()
