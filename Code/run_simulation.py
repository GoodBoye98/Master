import os
import sys
import time
import types
import string
import shutil
import random
import numpy as n
import pandas as pd
import subprocess
import matplotlib.pyplot as plt
import scipy.interpolate as interpolate
from matplotlib.widgets import Slider, Button
from matplotlib.pyplot import figure
from os import listdir, remove, mkdir, rmdir
from os.path import isfile, isdir, join
from multipledispatch import dispatch


############ GENERAL PARAMETERS ############

# Albedo and Q
a_0 = 0.6
a_1 = 0.38
Q = 500

# Heat diffusion
K_w = 0.38
# K_l = 0.38
K_l = 1.89

# Heat conductivity
C_w = 13.2
C_l = 10.2

# Continental edges
C_0 = -1.0
C_1 = 1.0

# Model parameters
A = 192.2
B = 3.85
C = 13.2

# Simulation domain
L = 5.0
M = 50000
N = 1024

# Extra
# dt = C / K_w * 0.5
# dt = 2.68
s = 1 / 40
dt = 1.45
nSavePoints = 76

times = None
run = []
T_vals = n.array([])
S_vals = n.array([])
perturbation = None

valid_runs = ["spectral", "finitediff", "midpoint"]

with open("Simulation Parameters/simulation.cfg", "r") as file:
    lines = [line.rstrip() for line in file]

with open("Simulation Parameters/simulation_particular.cfg", "r") as file:
    lines_2 = [line.rstrip() for line in file]

lines = lines + lines_2

for arg in lines:
    line = arg.split(" ")
    if line[0] == "Q":
        Q = float(line[1])
    elif line[0] == "K_w":
        K_w = float(line[1])
    elif line[0] == "K_l":
        K_l = float(line[1])
    elif line[0] == "C_0":
        C_0 = float(line[1])
    elif line[0] == "C_1":
        C_1 = float(line[1])
    elif line[0] == "L":
        L = float(line[1])
    elif line[0] == "s":
        s = float(line[1])
    elif line[0] == "N":
        N = int(line[1])
    elif line[0] == "time":
        times = []
        times_row = []
        for i in range(1, len(line)):
            if ';' in line[i]:
                times_row.append(float(line[i][:-1]))
                times.append(times_row)
                times_row = []
            else:
                times_row.append(float(line[i]))
        times.append(times_row)
    elif line[0] == "run":
        for i in range(1, len(line)):
            if line[i] in valid_runs:
                run.append(line[i])
    elif line[0] == "T":
        T_vals = n.array(line[1:], dtype=float)
    elif line[0] == "S":
        S_vals = n.array(line[1:], dtype=float)
    elif line[0] == "P":
        perturbation = interpolate.interp1d(n.linspace(-10, 10, 2001), n.array(line[1:], dtype=float))

times = n.array(times)

# print(f"Q = {Q}")
# print(f"K_w = {K_w}")
# print(f"K_l = {K_l}")
# print(f"C_0 = {C_0}")
# print(f"C_1 = {C_1}")
# print(f"L = {L}")
# print(f"N = {N}")
# print(f"time = {time}")
# print(f"times = {times}")
# print(f"run = {run}")
# print(f"T = {T_vals}")
# print(f"S = {S_vals}")

delete_leftovers = False
delete = True


lss = {'spectral': '-.',
       'finitediff': '--',
       'midpoint': '-.',
       'analytic': '-'}
lws = {'spectral': 1.75,
       'finitediff': 1.75,
       'midpoint': 1.75,
       'analytic': 1.75}
cs = {'spectral': 'Orange',
      'finitediff': 'Blue',
      'midpoint': 'Orange',
      'analytic': 'Black'}
labels = {'spectral': 'Spectral',
          'finitediff': 'Finite difference',
          'midpoint': 'Midpoint',
          'analytic': 'Analytic'}

############ ARTIFICIAL SOURCE DEFINITION ############

def Source_1(x, t):
    return n.cos(x * n.pi / 10) * n.cos(t)

def Source_2(x, t):
    m = 0.38 / 1.89
    idx = n.argwhere(n.abs(x) > 1)
    ret = - m * (x * x - 25) * n.cos(t)
    ret[idx] = (1 + 24 * m - x[idx] * x[idx]) * n.cos(t)
    return ret

def Source_3(x, t):
    return n.exp(- x ** 2) / (t + 1)

def Source_4(x, t):
    return n.cos(x + t)

############ CONFIGURATION CLASS ############

class Config:
    def __init__(self):
        self.a_0 = a_0
        self.a_1 = a_1
        self.Q = Q

        # Heat diffusion
        self.K_w = K_w
        self.K_l = K_l

        # Heat conductivity
        self.C_w = C_w
        self.C_l = C_l

        # Continental edges
        self.C_0 = C_0
        self.C_1 = C_1

        # Model parameters
        self.A = A
        self.B = B
        self.C = C

        # Simulation domain
        self.L = L
        self.M = M
        self.N = N

        # Extra
        self.s = s
        self.dt = dt
        self.nSavePoints = nSavePoints

        self.T = n.zeros(N)
        self.S = n.zeros(N)

        # Need to be defined
        self.name = None
        self.mode = None
        self.artificial_source = 0

    def write(self):
        with open(self.name, 'w') as f:
            f.write(f"a_0 {self.a_0}\n")
            f.write(f"a_1 {self.a_1}\n")
            f.write(f"Q {self.Q}\n")

            f.write(f"K_w {self.K_w}\n")
            f.write(f"K_l {self.K_l}\n")

            f.write(f"C_w {self.C_w}\n")
            f.write(f"C_l {self.C_l}\n")

            f.write(f"C_0 {self.C_0}\n")
            f.write(f"C_1 {self.C_1}\n")

            f.write(f"A {self.A}\n")
            f.write(f"B {self.B}\n")
            f.write(f"C {self.C}\n")

            f.write(f"L {self.L}\n")
            f.write(f"M {self.M}\n")
            f.write(f"s {self.s}\n")
            f.write(f"dt {self.dt}\n")
            f.write(f"mode {self.mode}\n")
            if self.artificial_source:
                f.write(f"artificial source {self.artificial_source}\n")
            f.write(f"savePoints {self.nSavePoints}\n")
            f.write('T(x) ' + ' '.join([f'{num:.15f}' for num in self.T]) + '\n')
            f.write('S(x) ' + ' '.join([f'{num:.15f}' for num in self.S]) + '\n')

############ PLOTTING CLASS ############

class PlotHandler:
    def __init__(self):
        pass

    @dispatch(list, types.FunctionType, list)
    def plotArtificialSource(self, modes, source, plot):
        # Create array of times used in plots
        times = n.array(plot)
        t_max = n.max(times)

        # Running simulation
        results = self.runArtificialSource(modes, source, t_max)

        times /= t_max
        fig, ax = plt.subplots(times.shape[0], times.shape[1], figsize=(8, 6))
        # fig.suptitle(f"Artificial source {source.__name__[-1]}")
        fig_2, ax_2 = plt.subplots(times.shape[0], times.shape[1], figsize=(8, 6))
        # fig_2.suptitle(f"Artificial source {source.__name__[-1]}")
        for row in range(times.shape[0]):
            for col in range(times.shape[1]):
                x = n.linspace(-L, L, N)
                t = list(results.values())[0][1]
                ax[row, col].plot(x, source(x, t[int(times[row, col] * (t.shape[0] - 1))]), ls=lss['analytic'], lw=lws['analytic'], c=cs['analytic'], label=labels['analytic'])
                for mode, params in sorted(results.items()):
                    i = int(times[row, col] * (params[2].shape[0] - 1))
                    ax[row, col].plot(params[0], params[2][i],
                        ls=lss[mode], lw=lws[mode], c=cs[mode], label=labels[mode])
                ax[row, col].set_title(f"t = {params[1][i]:.2f}")
                ax[row, col].set_ylim(self.plot_lims(params[2]))

                ax[row, col].set_xlabel("x")
                ax[row, col].set_ylabel("T(x, t)")

                ax[row, col].legend()


                L_2 = 1.5
                x_2 = n.linspace(-L_2, L_2, N)
                t = list(results.values())[0][1]
                ax_2[row, col].plot(x_2, source(x_2, t[int(times[row, col] * (t.shape[0] - 1))]), ls=lss['analytic'], lw=lws['analytic'], c=cs['analytic'], label=labels['analytic'])
                for mode, params in sorted(results.items()):
                    i = int(times[row, col] * (params[2].shape[0] - 1))
                    idxx = n.argwhere(abs(params[0]) < L_2)
                    ax_2[row, col].plot(params[0].flatten()[idxx], params[2][i].flatten()[idxx],
                        ls=lss[mode], lw=lws[mode], c=cs[mode], label=labels[mode])

                ax_2[row, col].set_title(f"t = {params[1][i]:.2f}")
                # ax_2[row, col].set_ylim(self.plot_lims(params[2]))

                ax_2[row, col].set_xlabel("x")
                ax_2[row, col].set_ylabel("T(x, t)")

                ax_2[row, col].legend()

        fig.tight_layout()
        fig_2.tight_layout()
        plt.show()

    @dispatch(list, types.FunctionType, float)
    def plotArtificialSource(self, modes, source, plot):

        # Running simulation
        results = self.runArtificialSource(modes, source, plot)
        times = n.array(plot)
        fig, ax = plt.subplots(figsize=(12, 8))
        fig.suptitle(f"Artificial source {source.__name__[-1]}")

        ax.set_xlabel("x")
        ax.set_ylabel("T(x, t)")

        x = n.linspace(-L, L, N)
        analytic, = plt.plot(x, source(x, 0), ls='--', lw=1.25, c='Black', label='Analytic')
        for mode, params in sorted(results.items()):
            if mode == 'spectral':
                spectral, = plt.plot(params[0], params[2][0],
                    ls=lss[mode], lw=lws[mode], c=cs[mode], label=labels[mode])
            if mode == 'finitediff':
                finitediff, = plt.plot(params[0], params[2][0],
                    ls=lss[mode], lw=lws[mode], c=cs[mode], label=labels[mode])
            if mode == 'midpoint':
                midpoint, = plt.plot(params[0], params[2][0],
                    ls=lss[mode], lw=lws[mode], c=cs[mode], label=labels[mode])


        ax.set_title(f"t = {0.0:.2f}")
        ax.set_ylim(self.plot_lims(params[2]))
        ax.legend()

        # Make a slider to control the time.
        ax_time = plt.axes([0.15, 0.02, 0.55, 0.03])
        time_slider = Slider(
            ax=ax_time,
            label='Time [t]',
            valmin=0,
            valmax=1,
            valinit=0)

        # The function to be called anytime a slider's value changes
        def update(val):
            for mode, params in sorted(results.items()):
                i = int(val * (params[2].shape[0] - 1))
                if mode == 'spectral':
                    spectral.set_ydata(params[2][i])
                if mode == 'finitediff':
                    finitediff.set_ydata(params[2][i])
                if mode == 'midpoint':
                    midpoint.set_ydata(params[2][i])
            analytic.set_ydata(source(x, params[1][i]))
            ax.set_title(f"t = {params[1][i]:.2f}")
            fig.canvas.draw_idle()

        time_slider.on_changed(update)

        # Create a `matplotlib.widgets.Button` to reset the sliders to initial values.
        resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
        button = Button(resetax, 'Reset', hovercolor='0.975')

        def reset(event):
            time_slider.reset()
        button.on_clicked(reset)

        # plt.tight_layout()
        plt.show()

    @dispatch(list, types.FunctionType, types.FunctionType, float, list, list)
    def plotSimulation(self, modes, T_0, S_x, Q, analytic, plot):
        # Create array of times used in plots
        times = n.array(plot)
        t_max = n.max(times)

        # Running simulation
        results = self.runSimulation(modes, T_0, S_x, Q, t_max)

        times /= t_max
        fig, ax = plt.subplots(times.shape[0], times.shape[1], figsize=(8, 6))
        # fig.suptitle(f"Artificial source {source.__name__[-1]}")
        for row in range(times.shape[0]):
            for col in range(times.shape[1]):
                x = n.linspace(-L, L, N)
                t = list(results.values())[0][1]
                i = int(times[row, col] * (analytic[1].shape[0] - 1))
                ax[row, col].plot(analytic[0], analytic[1][i], c='Black', label='Analytic')
                for mode, params in sorted(results.items()):
                    i = int(times[row, col] * (params[2].shape[0] - 1))
                    ax[row, col].plot(params[0], params[2][i],
                        ls=lss[mode], lw=lws[mode], c=cs[mode], label=labels[mode])
                ax[row, col].set_title(f"t = {params[1][i]:.2f}")
                ax[row, col].set_ylim(self.plot_lims(params[2]))

                ax[row, col].set_xlabel("x")
                ax[row, col].set_ylabel("T(x, t)")

                ax[row, col].legend()
        plt.tight_layout()
        plt.show()

    @dispatch(list, interpolate.interpolate.interp1d, interpolate.interpolate.interp1d, float, n.ndarray)
    def plotSimulation(self, modes, T_0, S_x, Q, plot):
        # Create array of times used in plots
        times = n.array(plot)
        t_max = n.max(times)

        # Running simulation
        results = self.runSimulation(modes, T_0, S_x, Q, t_max)

        times /= t_max
        fig, ax = plt.subplots(times.shape[0], times.shape[1], figsize=(8, 6))
        # fig.suptitle(f"Artificial source {source.__name__[-1]}")
        fig.suptitle(f"Q = {Q:.1f}")
        for row in range(times.shape[0]):
            for col in range(times.shape[1]):
                x = n.linspace(-L, L, N)
                t = list(results.values())[0][1]
                ax[row, col].plot(x, T_0(x) / 10, lw=1.75, c='Black', ls='--', label='Stationary')

                # data = pd.read_csv("WE/1 - Ice only.csv", header=None).values
                # data = data[n.argwhere(abs(data[:, 0] - Q) < 0.25)[0][0]][4:]
                # T_1 = interpolate.interp1d(n.linspace(-10, 10, 2001), data)
                # data = pd.read_csv("WE/2 - Partial ice.csv", header=None).values
                # data = data[n.argwhere(abs(data[:, 0] - Q) < 0.25)[-1][0]][4:]
                # T_2 = interpolate.interp1d(n.linspace(-10, 10, 2001), data)
                # ax[row, col].plot(x, T_1(x), lw=1.75, ls=':', c='Black', label='Stable solutions')
                # ax[row, col].plot(x, T_2(x), lw=1.75, ls=':', c='Black')

                for mode, params in sorted(results.items()):
                    i = int(times[row, col] * (params[2].shape[0] - 1))
                    ax[row, col].plot(params[0], params[2][i] / 10,
                        ls=lss[mode], lw=lws[mode], c=cs[mode], label=labels[mode])
                ax[row, col].set_title(f"t = {params[1][i]:.2f}")
                ax[row, col].set_ylim(self.plot_lims(params[2] / 10))

                ax[row, col].set_xlabel("x")
                ax[row, col].set_ylabel("T(x, t)")

                ax[row, col].legend()
        plt.tight_layout()
        plt.show()

    def runArtificialSource(self, modes, source, t_max):

        # Using discretization
        x = n.linspace(-L, L, N)
        spectral_x = n.linspace(-L - 1, L + 1, N, endpoint=False)

        dx = 2 * L / (N - 1)
        M = int(t_max / (dt * dx ** 2))

        # Running simulation for 20 units of time
        results = {}
        dirs = []

        for mode in modes:
            if mode == 'spectral':
                # Creating working directory
                dir = f"s_dir_{self.randomString(10)}";
                dirs.append(dir)
                if isdir(dir):
                    files = [f for f in listdir(dir) if isfile(join(dir, f))]
                    for file in files:
                        remove(dir + "/" + file)
                else:
                    mkdir(dir)

                # Calculatin new length to match other simulations
                dxx = 2 * (L + 1) / (N - 1)
                M_new = int(t_max / (dt * dxx ** 2))

                # Configuring simulation
                configuration = Config()
                configuration.name = dir + "/simulationConfig.cfg"
                configuration.mode = 'spectral'
                configuration.artificial_source = source.__name__[-1]
                configuration.T = source(spectral_x, 0)
                configuration.L = L + 1
                configuration.M = M_new
                if (source.__name__[-1] in "3 4"):
                    configuration.K_l = K_w
                configuration.write()

                # Running simulation in new directory
                shutil.copy('numerics.exe', dir + '/numerics.exe')
                subprocess.Popen([f"./numerics.exe"], cwd=dir)

            elif mode == 'finitediff':
                # Creating working directory
                dir = f"f_dir_{self.randomString(10)}";
                dirs.append(dir)
                if isdir(dir):
                    files = [f for f in listdir(dir) if isfile(join(dir, f))]
                    for file in files:
                        remove(dir + "/" + file)
                else:
                    mkdir(dir)

                # Configuring simulation
                configuration = Config()
                configuration.name = dir + "/simulationConfig.cfg"
                configuration.mode = 'finitediff'
                configuration.artificial_source = source.__name__[-1]
                configuration.T = source(x, 0)
                configuration.M = M
                configuration.write()

                # Running simulation in new directory
                shutil.copy('numerics.exe', dir + '/numerics.exe')
                subprocess.Popen([f"./numerics.exe"], cwd=dir)
            elif mode == 'midpoint':
              # Creating working directory
              dir = f"m_dir_{self.randomString(10)}";
              dirs.append(dir)
              if isdir(dir):
                  files = [f for f in listdir(dir) if isfile(join(dir, f))]
                  for file in files:
                      remove(dir + "/" + file)
              else:
                  mkdir(dir)

              # Configuring simulation
              configuration = Config()
              configuration.name = dir + "/simulationConfig.cfg"
              configuration.mode = 'midpoint'
              configuration.artificial_source = source.__name__[-1]
              configuration.T = source(x, 0)
              configuration.M = M
              configuration.write()

              # Running simulation in new directory
              shutil.copy('numerics.exe', dir + '/numerics.exe')
              subprocess.Popen([f"./numerics.exe"], cwd=dir)

        directories = [x[1] for x in os.walk("./")][0]
        while len(results) < len(modes):
            for mode in modes:
                for dir in dirs:
                    files = [f for f in listdir(dir) if isfile(join(dir, f))]
                    result = [file for file in files if mode + ".end" in file]

                    if result:
                        data = pd.read_csv(dir + "/" + result[0].replace(".end", ".sim"), header=None).values
                        t = data[:, 0]
                        simulation = data[:, 1:]
                        if mode == 'spectral':
                            x_idx = n.argwhere(n.abs(spectral_x) < L)
                            results[mode] = [spectral_x[x_idx].flatten(), t, simulation[:, x_idx]]
                        else:
                            results[mode] = [x, t, simulation]

                        if delete:
                            remove(dir + "/" + result[0])
                            remove(dir + "/" + result[0].replace(".end", ".sim"))
                    time.sleep(0.1)

        if delete:
            for dir in dirs:
                files = [f for f in listdir(dir) if isfile(join(dir, f))]
                for file in files:
                    remove(dir + "/" + file)
                rmdir(dir)

        return results

    def runSimulation(self, modes, T_0, S_x, Q, t_max):
        # Using discretization
        x = n.linspace(-L, L, N)
        spectral_x = n.linspace(-L - 1, L + 1, N, endpoint=False)

        # Computing number of iterations
        dx = 2 * L / (N - 1)
        M = int(t_max / (dt * dx ** 2))
        spectral_dx = 2 * (L + 1) / (N - 1)
        spectral_M = int(t_max / (dt * spectral_dx ** 2))

        # Running simulation for 20 units of time
        results = {}
        dirs = []

        for mode in modes:
            if mode == 'spectral':
                # Creating working directory
                dir = f"s_dir_{self.randomString(10)}";
                dirs.append(dir)
                if isdir(dir):
                    files = [f for f in listdir(dir) if isfile(join(dir, f))]
                    for file in files:
                        remove(dir + "/" + file)
                else:
                    mkdir(dir)

                # Configuring simulation
                configuration = Config()
                configuration.name = dir + "/simulationConfig.cfg"
                configuration.mode = 'spectral'
                configuration.Q = Q
                configuration.T = T_0(spectral_x)
                if perturbation is not None:
                    configuration.T += perturbation(spectral_x)
                configuration.S = S_x(spectral_x)
                configuration.L = L + 1
                configuration.M = spectral_M
                configuration.write()

                # Running simulation in new directory
                shutil.copy('numerics.exe', dir + '/numerics.exe')
                subprocess.Popen([f"./numerics.exe"], cwd=dir)

            elif mode == 'finitediff':
                # Creating working directory
                dir = f"f_dir_{self.randomString(10)}";
                dirs.append(dir)
                if isdir(dir):
                    files = [f for f in listdir(dir) if isfile(join(dir, f))]
                    for file in files:
                        remove(dir + "/" + file)
                else:
                    mkdir(dir)

                # Configuring simulation
                configuration = Config()
                configuration.name = dir + "/simulationConfig.cfg"
                configuration.mode = 'finitediff'
                configuration.Q = Q
                configuration.T = T_0(x)
                if perturbation is not None:
                    configuration.T += perturbation(x)
                configuration.S = S_x(x)
                configuration.M = M
                configuration.write()

                # Running simulation in new directory
                shutil.copy('numerics.exe', dir + '/numerics.exe')
                subprocess.Popen([f"./numerics.exe"], cwd=dir)
            elif mode == 'midpoint':
                # Creating working directory
                dir = f"m_dir_{self.randomString(10)}";
                dirs.append(dir)
                if isdir(dir):
                  files = [f for f in listdir(dir) if isfile(join(dir, f))]
                  for file in files:
                      remove(dir + "/" + file)
                else:
                  mkdir(dir)

                # Configuring simulation
                configuration = Config()
                configuration.name = dir + "/simulationConfig.cfg"
                configuration.mode = 'midpoint'
                configuration.Q = Q
                configuration.T = T_0(x)
                if perturbation is not None:
                    configuration.T += perturbation(x)
                configuration.S = S_x(x)
                configuration.M = M
                configuration.write()

                # Running simulation in new directory
                shutil.copy('numerics.exe', dir + '/numerics.exe')
                subprocess.Popen([f"./numerics.exe"], cwd=dir)

        while len(results) < len(modes):
            for mode in modes:
                for dir in dirs:
                    files = [f for f in listdir(dir) if isfile(join(dir, f))]
                    result = [file for file in files if mode + ".end" in file]

                    if result:
                        data = pd.read_csv(dir + "/" + result[0].replace(".end", ".sim"), header=None).values
                        t = data[:, 0]
                        simulation = data[:, 1:]
                        if mode == 'spectral':
                            x_idx = n.argwhere(n.abs(spectral_x) < L)
                            results[mode] = [spectral_x[x_idx].flatten(), t, simulation[:, x_idx]]
                        else:
                            results[mode] = [x, t, simulation]

                        if delete:
                            remove(dir + "/" + result[0])
                            remove(dir + "/" + result[0].replace(".end", ".sim"))
                    time.sleep(0.1)

        if delete:
            for dir in dirs:
                files = [f for f in listdir(dir) if isfile(join(dir, f))]
                for file in files:
                    remove(dir + "/" + file)
                rmdir(dir)

        return results

    def randomString(self, length):
        return ''.join(random.choices(string.ascii_lowercase + string.digits, k=length))

    def plot_lims(self, data):
        y_min = n.min(data)
        y_max = n.max(data)

        diff = y_max - y_min

        return (y_min - diff * 0.05, y_max + diff * 0.05)


def main():

    # Deletes directories with simulations
    if delete_leftovers:
        dirs = os.listdir()
        for dir in dirs:
            if "s_dir_" in dir or "f_dir_" in dir or "m_dir_" in dir:
                files = [f for f in listdir(dir) if isfile(join(dir, f))]
                for file in files:
                    remove(dir + "/" + file)
                rmdir(dir)

    T_init = interpolate.interp1d(n.linspace(-10, 10, 2001), T_vals)
    S_x = interpolate.interp1d(n.linspace(-10, 10, 2001), S_vals)

    handler = PlotHandler()

    # if times.shape == (1, 1):
    #     handler.plotSimulation(run, T_init, S_x, Q, times[0][0])
    # if times is not None:
    #     handler.plotSimulation(run, T_init, S_x, Q, times)

    handler.plotArtificialSource(
        ["spectral"], Source_2, [[0, 2.5], [5.0, 7.5]]
    )
    #
    # # handler.plotArtificialSource(
    # #     ["spectral"], Source_t, [[0.0, 3.0], [6.0, 10.0]]
    # # )
    #
    # analytic_x = n.linspace(-5, 5, 101)
    # analytic_T = pd.read_csv("Analytic_Sol.csv", header=None).values
    # analytic_T = analytic_T.reshape((101, 101))
    #
    # handler.plotSimulation(['spectral'], T_0, S_gaussian, 340.0,
    #     [analytic_x, analytic_T], [[0.0, 3.0], [6.0, 10.0]])

    # handler.plotSimulation(['spectral'], T_0, S_gaussian, 340.0, [[0.0, 3.0], [6.0, 10.0]])


if __name__ == "__main__":
    main()
