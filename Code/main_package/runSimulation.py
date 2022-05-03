import sys
import time
import shutil
import subprocess
import numpy as n
import pandas as pd
import scipy.integrate as integrate
import scipy.interpolate as interpolate
import matplotlib.pyplot as plt
from os.path import isfile, isdir, join
from os import listdir, remove, mkdir, rmdir, getcwd

try:
    from Miscellaneous import randomString
except ModuleNotFoundError:
    from main_package.Miscellaneous import randomString


class runSimulation:
    def __init__(self, mode, data, dir, S):
        self.mode = mode
        self.data = data
        self.dir = dir
        self.S = eval(S)

        self.simRunning = False

        if self.mode == 'spectral':
            self.simDir = f'{self.dir}\\s_dir_{randomString(10)}'
        elif self.mode == 'finitediff':
            self.simDir = f'{self.dir}\\f_dir_{randomString(10)}'
        else:
            raise ValueError(f'mode must be either "spectral" or "finitediff"')
        self.simDir = self.simDir if self.dir else self.simDir[1:]

        # Copy numerics to new directory
        mkdir(self.simDir)
        shutil.copy('C++\\numerics.exe', self.simDir + '\\numerics.exe')

        self.Q = self.data[0]
        self.C_0 = self.data[2]
        self.C_1 = self.data[3]
        self.T = self.data[4:]
        self.T = interpolate.interp1d(n.linspace(-6, 6, self.T.shape[0]), self.T)

        self._N = 2048
        self._L = 6.0
        self._s = 0.025
        self.maxT = 5
        self.savePoints = 1000

        self.sim = None
        self.result = None
        self.times = None

    @property
    def N(self):
        return self._N
    @N.setter
    def N(self, val):
        if val <= 0:
            raise ValueError('N set out of bounds, must be integer over 0')
        else:
            self._N = val
            
    @property
    def L(self):
        return self._L
    @L.setter
    def L(self, val):
        if val <= 0.0 or val > 6.0:
            raise ValueError('L set out of bounds, must be in half open interval (0.0, 6.0]')
        else:
            self._L = val
            
    @property
    def s(self):
        return self._s
    @s.setter
    def s(self, val):
        if val <= 0:
            raise ValueError('s set out of bounds, must be float over 0.0')
        else:
            self._s = val

    def startSimulation(self, perturbation=0):
        # Starting parameters
        if self.mode == "spectral":
            self.x = n.linspace(-self.L, self.L, self._N, endpoint=False)
            dx = self.x[1] - self.x[0]
            dt = 1.45
            M = int(self.maxT / (dt * dx ** 2))
        if self.mode == "finitediff":
            self.x = n.linspace(-self.L, self.L, self._N, endpoint=True)
            dx = self.x[1] - self.x[0]
            dt = 3.48
            M = int(self.maxT / (dt * dx ** 2))

        T = self.T(self.x) + perturbation
        S = self.S(self.x)
        F = self._discretizedGaussian(dx)

        self.T_vals = T - perturbation

        # Writes parameters to file
        with open(f'{self.simDir}/simulationConfig.cfg', 'w') as f:
            f.write(f'mode {self.mode}\n')
            f.write(f'savePoints {self.savePoints}\n')
            f.write(f'Q {self.Q}\n')
            f.write(f'C_0 {self.C_0}\n')
            f.write(f'C_1 {self.C_1}\n')
            f.write(f'N {self._N}\n')
            f.write(f'M {M}\n')
            f.write(f'L {self._L}\n')
            f.write(f's {self._s}\n')
            f.write(f'dt {dt}\n')
            f.write('T(x) ' + ' '.join([f'{num:.15f}' for num in T]) + '\n')
            f.write('S(x) ' + ' '.join([f'{num:.15f}' for num in S]) + '\n')
            f.write('F(x) ' + ' '.join([f'{num:.15f}' for num in F]) + '\n')

        # Runs simulation
        self.sim = subprocess.Popen([f".\\numerics.exe"], cwd=self.simDir, shell=True)

    def runArtificialSource(self, num):
        # Starting parameters
        if self.mode == "spectral":
            self.x = n.linspace(-self.L, self.L, self._N, endpoint=False)
            dx = self.x[1] - self.x[0]
            dt = 1.45
            M = int(self.maxT / (dt * dx ** 2))
        if self.mode == "finitediff":
            self.x = n.linspace(-self.L, self.L, self._N, endpoint=True)
            dx = self.x[1] - self.x[0]
            dt = 3.48
            M = int(self.maxT / (dt * dx ** 2))

        T = self.T(self.x)
        F = self._discretizedGaussian(dx)

        # Writes parameters to file
        with open(f'{self.simDir}/simulationConfig.cfg', 'w') as f:
            f.write(f'mode {self.mode}\n')
            f.write(f'artificial source {num}\n')
            f.write(f'C_0 {self.C_0}\n')
            f.write(f'C_1 {self.C_1}\n')
            f.write(f'N {self._N}\n')
            f.write(f'M {M}\n')
            f.write(f'L {self._L}\n')
            f.write(f's {self._s}\n')
            f.write(f'dt {dt}\n')
            f.write('T(x) ' + ' '.join([f'{num:.15f}' for num in T]) + '\n')
            f.write('F(x) ' + ' '.join([f'{num:.15f}' for num in F]) + '\n')

        # Runs simulation
        self.sim = subprocess.Popen([f".\\numerics.exe"], cwd=self.simDir, shell=True)

    def _discretizedGaussian(self, dx):
        num = int(n.ceil(4 * self.s / dx))
        x = n.linspace(-dx * num, dx * num, 2 * num + 1)
        weights = []
        for i in range(x.shape[0]):
            weights.append(integrate.quad(lambda x: 1 / (self.s * n.sqrt(n.pi)) * n.exp(- (x / self.s) ** 2), x[i] - 1/2*dx, x[i] + 1/2*dx)[0])
        weights = n.array(weights)
        weights /= n.sum(weights)
        return weights
    
    def getProgress(self):
        try:
            # Open progress-file and read content
            s = ''
            with open(f'{self.simDir}/progress.txt', 'r') as f:
                s = f.read()

            # If content is in order, return progress
            if len(s) and s[-1] == 'f':
                return float(s.replace('f','')) / 10
            # Otherwise return False
            return False
        except FileNotFoundError as e:
            return 0

    def getResultFirst(self, rows=[]):
        if self.result is None:
            files = [f for f in listdir(self.simDir) if isfile(join(self.simDir, f))]
            if 'result.end' in files:
                if rows:
                    rows = n.array(rows) % 2
                    skipRows = [i for i in range(2) if i not in rows]
                    self.result = pd.read_csv(f'{self.simDir}\\result.sim', header=None, skiprows=skipRows).values
                else:
                    self.result = pd.read_csv(f'{self.simDir}\\result.sim', header=None).values
                return True
            else:
                return False
        else:
            return False

    def getResult(self):
        return self.x, self.result

    def cleanup(self):
        shutil.rmtree(self.simDir)


def main():
    analytic_x = n.linspace(-5, 5, 101)
    analytic_T = pd.read_csv("Analytic_Sol.csv", header=None).values
    analytic_T = analytic_T.reshape((101, 101))


    # Artificial Sources to choose from
    def a_source_1(x, t):
        return n.cos(x * n.pi / 10) * n.cos(t)
    def a_source_2(x, t):
        return n.cos(t) * n.where(n.abs(x) < 1, -0.38 / 1.89 * (x ** 2 - 25), 1 + 24 * 0.38 / 1.89 - x ** 2)
    def a_source_3(x, t):
        return n.exp(- x ** 2) * 1 / (t + 1)
    def a_source_4(x, t):
        return n.cos(x + t)
    def analytic(x, t):
        idx = int(100 * t / 10)
        interpolating_func = interpolate.interp1d(analytic_x, analytic_T[idx], fill_value='extrapolate')
        return interpolating_func(x)

    # Pick
    a_source = a_source_2
    plotTimes = n.array([[0.0, 3.33], [6.67, 10.0]])
    mode = "finitediff"
    zoomed = True

    ## Running artificial source ##
    init = n.array([0, 0, 0, 0])
    x = n.linspace(-6, 6, 12001)
    data = a_source(x, 0)

    data = n.concatenate((init, data))
    obj = runSimulation(mode, data, "", "None")
    obj.N = 4001
    obj.L = 5
    obj.s = 0.005
    obj.maxT = 10.0
    obj.runArtificialSource(a_source.__name__[-1])
    # obj.startSimulation()

    # Waiting for C++ to finish simulation
    while True:
        a = obj.getProgress()
        if a:
            sys.stdout.write("\b" * 20 + str(a))
            sys.stdout.flush()
        if obj.getResultFirst():
            break
        time.sleep(0.1)
    
    # Load data
    x, results = obj.getResult()
    obj.cleanup()

    # Definit visual parameters
    lss = {'spectral': 'dashdot',
           'finitediff': (0, (3, 1, 1, 1, 1, 1)),
           'analytic': 'solid'}
    lws = {'spectral': 1.75,
           'finitediff': 1.75,
           'analytic': 1.75}
    cs = {'spectral': 'darkorange',
          'finitediff': 'deepskyblue',
          'analytic': 'Black'}
    labels = {'spectral': 'Spectral',
              'finitediff': 'Finite difference',
              'analytic': 'Analytic'}

    y_max = n.max(results[:, 1:])
    y_min = n.min(results[:, 1:])
    y_dif = y_max - y_min

    # Plotting
    fig, ax = plt.subplots(2, 2, figsize=(8, 6))
    idx = n.argwhere(n.abs(x) < 5)
    for i in range(2):
        for j in range(2):
            t_idx = int(plotTimes[i, j] / n.max(plotTimes) * (results.shape[0] - 1))
            # t = plotTimes[i, j]
            t = results[t_idx, 0]
            ax[i, j].set_ylim((y_min - 0.05 * y_dif, y_max + 0.05 * y_dif))
            ax[i, j].set_title(fr"$t={t:.2f}$")
            ax[i, j].plot(x[idx], a_source(x, t)[idx], ls=lss['analytic'], lw=lws['analytic'], c=cs['analytic'], label=labels['analytic'])
            ax[i, j].plot(x[idx], results[t_idx, 1:][idx], ls=lss[mode], lw=lws[mode], c=cs[mode], label=labels[mode])
            ax[i, j].legend()
    plt.tight_layout()

    if zoomed:
        fig_2, ax_2 = plt.subplots(2, 2, figsize=(8, 6))
        idx = n.argwhere(n.abs(x) < 1.5)
        for i in range(2):
            for j in range(2):
                t_idx = int(plotTimes[i, j] / n.max(plotTimes) * (results.shape[0] - 1))
                # t = plotTimes[i, j]
                t = results[t_idx, 0]
                ax_2[i, j].set_title(fr"$t={t:.2f}$")
                ax_2[i, j].plot(x[idx], a_source(x, t)[idx], ls=lss['analytic'], lw=lws['analytic'], c=cs['analytic'], label=labels['analytic'])
                ax_2[i, j].plot(x[idx], results[t_idx, 1:][idx], ls=lss[mode], lw=lws[mode], c=cs[mode], label=labels[mode])
                ax_2[i, j].legend()
        plt.tight_layout()

    plt.show()


if __name__ == "__main__":
    main()
    