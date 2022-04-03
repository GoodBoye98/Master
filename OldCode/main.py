import time
import subprocess
import numpy as n
import pandas as pd
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import scipy.interpolate as interpolate
from os import listdir, remove, mkdir, rmdir
from os.path import isfile, isdir, join
from matplotlib.widgets import Slider
from Supplementary.make_video import makeVideo


solutionFamilies = []


class solutionFamily:
    def __init__(self, dir, filename):
        """
        Load a family of solutions from filename and directory.
        """
        data = pd.read_csv(dir + '/' + filename, header=None, usecols=[0, 1], delimiter=',').values
        self.name = filename[filename.find('-') + 2:].replace('.csv', '')
        self.dir = dir
        self.solutions = []
        self.C_0 = -1
        self.C_1 = 1

        # Fix offset and s(x)
        if dir[0] == 'G':
            if 'W' in dir:
                self.C_0 = 0
                self.C_1 = 0
            val = float(dir[dir.find('G')+1:]) if 'W' not in dir else 10
            S = lambda x : 4 / (n.sqrt(n.pi * val)) * n.exp(- x * x / val)
        else:
            if 'W' not in dir:
                self.C_0 += float(dir[2:])
                self.C_1 += float(dir[2:])
            else:
                self.C_0 = 0
                self.C_1 = 0
            S = lambda x : n.exp(-n.abs(x) / 2)

        for i in range(data.shape[0]):
            sol = solution(data[i])
            sol.C_0 = self.C_0
            sol.C_1 = self.C_1
            sol.name = self.name
            sol.S = S
            self.solutions.append(sol)

    def getSolCoordinates(self, measure='T_0'):
        """
        Get coordinates to all solutions in this family
        with a chosen measure.
        T_0: temperature at x=0.
        T_center: temerature at the center of the continent.
            (Same as T_0 in the case of no or centered continent)
        T_max: max temperature on the line or continent if there is one
        """
        points = n.zeros((len(self.solutions), 2))
        for i, sol in enumerate(self.solutions):
            points[i] = sol.getPoint(measure)
        return points

    def addToDisplay(self, ax, measure='T_0'):
        """
        Add all solutions in this family to plot.
        T_0: temperature at x=0.
        T_center: temerature at the center of the continent.
            (Same as T_0 in the case of no or centered continent)
        T_max: max temperature on the line or continent if there is one
        """
        points = self.getSolCoordinates(measure)
        if "IceWaterSnowLandSnowWaterIce" in self.name:
            points = points[n.argsort(points[:, 0])]
        else:
            points = points[n.argsort(points[:, 1])]
        if self.name == 'Partial ice':
            if self.dir == "E W":
                up = -0.5
                low = -0.52
            elif self.dir == "G W":
                up = -0.4
                low = -0.42
            idx = n.argwhere(points[:, 1] > low)
            ax.plot(points[idx,0], points[idx,1], label=self.name)
            idx = n.argwhere(points[:, 1] < up)
            ax.plot(points[idx,0], points[idx,1], label=self.name + " unstable", color='Black')
        elif self.dir == 'G 6.0':
            ax.plot(points[:,0], points[:,1])
        else:
            ax.plot(points[:,0], points[:,1], label=self.name)



class solution:
    def __init__(self, solInfo):
        """
        Helper class for a single solution in a family.
        """
        self.Q_point = solInfo[0]
        self.T_0 = solInfo[1]
        self.solution = interpolate.interp1d(n.linspace(-6, 6, solInfo[2:].shape[0]), solInfo[2:])
        # self.solution = solInfo[4:]
        self.L = 5
        self.C_0 = -1
        self.C_1 = 1

    def getPoint(self, measure):
        """
        Get coordinates of single solution based upon chosen merit.
        T_0: temperature at x=0.
        T_center: temerature at the center of the continent.
            (Same as T_0 in the case of no or centered continent)
        T_max: max temperature on the line or continent if there is one
        """
        if measure == 'T_0':
            return n.array([self.Q_point, self.T_0])
        if measure == 'T_max':
            return n.array([self.Q_point, self.T_max])
        if measure == 'T_center':
            return n.array([self.Q_point, self.T_center])

    def runSimulation(self):
        pert = 0
        x = n.linspace(-10, 10, 2001)
        valsT = 10 * self.solution(x)
        valsS = self.S(x)
        valsP = pert * n.ones(2001)

        K_l = 1.89 if self.C_0 != self.C_1 else 0.38

        # Writes information particular to this simulation
        with open("Simulation Parameters/simulation_particular.cfg", "w") as file:
            file.write(f'Q {self.Q_point}\n')
            file.write(f'C_0 {self.C_0}\n')
            file.write(f'C_1 {self.C_1}\n')
            file.write('T ' + ' '.join([f'{num:.10f}' for num in valsT]) + '\n')
            file.write('S ' + ' '.join([f'{num:.10f}' for num in valsS]) + '\n')
            file.write('P ' + ' '.join([f'{num:.10f}' for num in valsP]) + '\n')
            file.truncate()

        # Sends information to script that runs simulation
        subprocess.Popen(
        ['python',
        'run_simulation.py'], shell=True)


def stabilityTest(sol, ε=-0.15):
    Q = sol.Q_point
    A = 192.2
    B = 3.85
    C = 13.2
    K_w = 0.38
    K_l = 1.89
    a_0 = 0.6
    a_1 = 0.38

    N = 1001
    x = n.linspace(-5, 5, N)
    dx = x[1]-x[0]
    T = sol.solution(x)
    S = sol.S(x)
    ε = ε * n.ones(N)

    a = n.zeros(N)
    K = n.zeros(N)
    for i in range(N):
        if sol.C_0 < x[i] < sol.C_1:
            a[i] = (a_1 - a_0) * 0.5 / dx if T[i] * T[i+1] < 0 or T[i] * T[i-1] < 0 else 0
            K[i] = K_l
        else:
            if i == 0 or i == N-1:
                a[i] = 0
            else:
                a[i] = (a_1 - a_0) * 0.5 / dx if (T[i]+1) * (T[i+1]+1) < 0 or (T[i]+1) * (T[i-1]+1) < 0 else 0
            K[i] = K_w
    a = smooth(a, dx, 0.025)
    K = smooth(K, dx, 0.025)

    # Comuting N'(T)
    # a = (a[2:] - a[:-2]) / (2 * dx)
    F = -Q * S * a

    plt.ion()

    fig, ax = plt.subplots(1, 1)

    # ax.plot(x, T)
    vals, = ax.plot(x, ε, label='ε')
    ax.set_ylim((-0.5, 0.5))
    ax.legend()
    plt.tight_layout()
    fig.canvas.draw()
    fig.canvas.flush_events()

    dt = dx ** 2
    dε = n.zeros(N)
    while True:
        for j in range(N):
            if j == 0:
                ε_j1 = 3 * ε[j] - 3 * ε[j+1] + ε[j+2]
                dε[j] = dt / C * (  1 / dx ** 2 * (K[j] * (ε[j+1] - ε[j]) - K[j-0] * (ε[j] - ε_j1) + 
                                                   K[j+1] * (ε[j+1] - ε[j]) - K[j] * (ε[j] - ε_j1))
                                  + (F[j] - B) * ε[j])
            elif j == N - 1:
                ε_j1 = 3 * ε[j] - 3 * ε[j-1] + ε[j-2]
                dε[j] = dt / C * (  1 / dx ** 2 * (K[j] * (ε_j1 - ε[j]) - K[j-1] * (ε[j] - ε[j-1]) + 
                                                   K[j+0] * (ε_j1 - ε[j]) - K[j] * (ε[j] - ε[j-1]))
                                  + (F[j] - B) * ε[j])
            else:
                dε[j] = dt / C * (  1 / dx ** 2 * (K[j] * (ε[j+1] - ε[j]) - K[j-1] * (ε[j] - ε[j-1]) + 
                                                   K[j+1] * (ε[j+1] - ε[j]) - K[j] * (ε[j] - ε[j-1]))
                                  + (F[j] - B) * ε[j])
        
        ε += dε

        if i % 100 == 0:
            print(n.sum(ε * dx))
            vals.set_ydata(ε)

            fig.canvas.draw()
            fig.canvas.flush_events()
      

def gauss(x, s):
    return 1 / (s * n.sqrt(n.pi)) * n.exp(- (x / s) ** 2)


def convolutionWeights(dx, s):
    num = int(n.ceil(5 * s / dx))
    x = n.linspace(-dx * num, dx * num, 2 * num + 1)
    weights = []
    for i in range(x.shape[0]):
        weights.append(integrate.quad(lambda x: gauss(x, s), x[i] - 1/2*dx, x[i] + 1/2*dx)[0])
    weights = n.array(weights)
    weights /= n.sum(weights)
    return weights


def smooth(input, dx, s):

    conv_weights = convolutionWeights(dx, s)
    length = conv_weights.shape[0]

    pad_l = n.ones(length // 2) * input[0]
    pad_r = n.ones((length - 1) // 2) * input[-1]
    input = n.concatenate((pad_l, input, pad_r))
    input = n.convolve(input, conv_weights, 'valid')

    return input


def plotSolutions(dirName):

    # Finding what S(x) to put in title of plot
    if dirName[0] == 'G':
        val = float(dirName[dirName.find('G')+1:]) if 'W' not in dirName else 10
        sFunction = rf"$\frac{{4}}{{\sqrt{{ {val:.1f} \pi}}}} e^{{-\frac{{x^2}}{{ {val:.1f} }}}}$"
    else:
        sFunction = r"$e^{-\frac{|x|}{2}}$"

    fig, ax = plt.subplots(1, 1)
    # fig.suptitle(f"Continent at [{cnt[0]}, {cnt[1]}] and S(x)={s_function}")
    fig.suptitle(f"Bifurcation diagram with S(x)={sFunction}")
    ax.set_xlabel("Q")
    ax.set_ylabel("T(0)")

    point = plt.scatter([n.nan], [n.nan], s=5, facecolors='none', edgecolors='Black')

    for family in solutionFamilies:
        family.addToDisplay(ax, measure='T_0')

    xlim = ax.get_xlim()
    xL = xlim[1] - xlim[0]
    ylim = ax.get_ylim()
    yL = ylim[1] - ylim[0]

    allSols = [s for sols in solutionFamilies for s in sols.solutions]
    solsPos = n.array([sol.getPoint('T_0') for sol in allSols])
    solsPos[:, 0] = (solsPos[:, 0] - xlim[0]) / xL
    solsPos[:, 1] = (solsPos[:, 1] - ylim[0]) / yL

    def onclick(event):
        """
        Event handeler for clicking on a datapoint
        """
        if event.button == 3:
            # Finding family closest to point
            mousePos = n.array([event.xdata, event.ydata])
            mousePos[0] = (mousePos[0] - xlim[0]) / xL
            mousePos[1] = (mousePos[1] - ylim[0]) / yL

            # Computing distances to mouse pos
            dists = n.sum((solsPos - mousePos) ** 2, axis=1)

            # The closes solution to click if it's close enough
            closest = allSols[n.argmin(dists)]
            # stabilityTest(closest)
            closest.runSimulation()


    def onmove(event):
        """
        Event handeler for clicking on a datapoint
        """
        # Finding family closest to point
        try:
            mousePos = n.array([event.xdata, event.ydata])
            mousePos[0] = (mousePos[0] - xlim[0]) / xL
            mousePos[1] = (mousePos[1] - ylim[0]) / yL

            # Computing distances to mouse pos
            dists = n.sum((solsPos - mousePos) ** 2, axis=1)

            # The closes solution to click if it's close enough
            closest = allSols[n.argmin(dists)]
            point.set_offsets([closest.Q_point, closest.T_0])
            fig.canvas.draw_idle()
        except TypeError:
            pass

    subprocess.Popen(
    ['python',
    'simulation_parameters.py'], shell=True)

    cid = fig.canvas.mpl_connect('button_press_event', onclick)
    # cid = fig.canvas.mpl_connect('scroll_event', onmove)
    cid = fig.canvas.mpl_connect('motion_notify_event', onmove)
    ax.legend()
    plt.show()


def main():

    # Choosing what folder to gather data from
    dirName = 'G W'

    # Loading solution data from folder
    filenames = [f for f in listdir(dirName) if isfile(join(dirName, f))]

    for solName in filenames:
        with open(dirName + '/' + solName, 'r+') as file:
            content = file.read()
            if len(content) == 0:
                continue
            content = content.replace('{', '').replace('}', '').replace('"', '').replace('*^', 'e')
            file.seek(0)
            file.write(content)
            file.truncate()

        newFamily = solutionFamily(dirName, solName)
        solutionFamilies.append(newFamily)

    # Plotting data from folder
    plotSolutions(dirName)


def createVideo():
    # Create a video of data in specified folder
    videoMaker = makeVideo("E 0.00")
    videoMaker.generate("Videos", "E 0.00")


if __name__ == '__main__':
    main()
    # createVideo()
