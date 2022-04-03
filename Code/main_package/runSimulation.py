import shutil
import subprocess
import numpy as n
import pandas as pd
import scipy.integrate as integrate
import scipy.interpolate as interpolate
from Miscellaneous import randomString
from os.path import isfile, isdir, join
from os import listdir, remove, mkdir, rmdir, getcwd

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

        # Copy numerics to new directory
        mkdir(self.simDir)
        shutil.copy('C++\\numerics.exe', self.simDir + '\\numerics.exe')

        self.Q = self.data[0]
        self.C_0 = self.data[2]
        self.C_1 = self.data[3]
        self.T = self.data[4:]
        self.T = interpolate.interp1d(n.linspace(-6, 6, self.T.shape[0]), self.T)

        self._N = 1024
        self._L = 5.0
        self._s = 0.025

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


    def startSimulation(self):
        # Starting parameters
        if self.mode == "spectral":
            self.x = n.linspace(-self.L, self.L, self.N, endpoint=False)
            dx = self.x[1] - self.x[0]
            dt = 1.45
            M = int(5 / (dt * dx ** 2))
        if self.mode == "finitediff":
            self.x = n.linspace(-self.L, self.L, self.N, endpoint=True)
            dx = self.x[1] - self.x[0]
            dt = 3.48
            M = int(5 / (dt * dx ** 2))

        T = self.T(self.x)
        S = self.S(self.x)
        F = self._discretizedGaussian(dx)

        # Writes parameters to file
        with open(f'{self.simDir}/simulationConfig.cfg', 'w') as f:
            f.write(f'mode {self.mode}\n')
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


    def _discretizedGaussian(self, dx):
        num = int(n.ceil(5 * self.s / dx))
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
                return float(s[:-1]) / 10
            # Otherwise return False
            return False
        except FileNotFoundError as e:
            return 0


    def getResultFirst(self):
        if self.result is None:
            files = [f for f in listdir(self.simDir) if isfile(join(self.simDir, f))]
            if 'result.end' in files:
                self.result = pd.read_csv(f'{self.simDir}\\result.sim', header=None).values
                return True
            else:
                return False
        else:
            return False


    def getResult(self):
        return self.x, self.result

