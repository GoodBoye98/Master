import time
import subprocess
import numpy as n
import pandas as pd
import matplotlib.pyplot as plt
import scipy.interpolate as interpolate
from os import listdir, remove, mkdir, rmdir
from os.path import isfile, isdir, join
from matplotlib.widgets import Slider
from Supplementary.make_video import makeVideo


solutionFamilies = []


class solFamily:
    def __init__(self, dir, cedges, filename):
        """
        Load a family of solutions from filename and directory.
        """
        data = pd.read_csv(dir + '/' + filename, header=None, delimiter=',').values
        try:
            self.offset = float(dir)
        except ValueError:
            self.offset = 0
        self.name = filename[filename.find('-') + 2:].replace('.csv', '')
        self.dir = dir
        self.solutions = []
        self.C_0 = cedges[0]
        self.C_1 = cedges[1]
        for i in range(data.shape[0]):
            sol = solution(data[i])
            sol.C_0 = self.C_0
            sol.C_1 = self.C_1
            sol.name = self.name
            if 'G' in dir:
                val = float(dir[2:]) if dir[0] == 'G' else 10
                sol.S = lambda x : 4 / (n.sqrt(n.pi * val)) * n.exp(- x * x / val)
            else:
                sol.S = lambda x : n.exp(-n.abs(x) / 2)
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
            if self.dir == "WE":
                up = -0.5
                low = -0.52
            elif self.dir == "WG":
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

    def minDist(self, point, measure='T_0'):
        """
        Get minimum distance from a point to any of solution in
        this family.
        T_0: temperature at x=0.
        T_center: temerature at the center of the continent.
            (Same as T_0 in the case of no or centered continent)
        T_max: max temperature on the line or continent if there is one
        """
        points = self.getSolCoordinates(measure)
        dists = n.sum((points - point) ** 2, axis=1)
        return n.sqrt(n.min(dists))


class solution:
    def __init__(self, solInfo):
        """
        Helper class for a single solution in a family.
        """
        self.Q_point = solInfo[0]
        self.T_max = solInfo[1]
        self.T_0 = solInfo[2]
        self.T_center = solInfo[3]
        self.solution = interpolate.interp1d(n.linspace(-10, 10, 2001), solInfo[4:])
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

    def runFiniteDifference(self):
        pert = 0
        N = 1024
        x = n.linspace(-10, 10, 2001)
        valsT = 10 * self.solution(x)
        valsS = self.S(x)
        valsP = pert * n.cos(0.5 * x * n.pi / self.L)

        K_l = 1.89 if self.C_0 != self.C_1 else 0.38

        # Adds information particular to this simulation
        with open("Simulation Parameters/simulation_particular.cfg", "a") as file:
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

        # if (self.C_0 == self.C_1):
        #     vals3 = 0.38 * n.ones(1001)
        #
        # with open('init_fdPlot.cfg', 'w') as f:
        #     f.write(f"name {self.name}\n")
        #     f.write(f"Q {self.Q_point}\n")
        #     f.write(f"L {self.L}\n")
        #     f.write(f"C_0 {self.C_0}\n")
        #     f.write(f"C_1 {self.C_1}\n")
        #     f.write(f"M {500001}\n")
        #     f.write(f"N {1001}\n")
        #     f.write(f"s {100}\n")
        #     f.write('T(x) ' + ' '.join([f'{num:.10f}' for num in vals]) + '\n')
        #     f.write('S(x) ' + ' '.join([f'{num:.10f}' for num in vals2]) + '\n')
        #
        #
        # subprocess.Popen(["python3", "fdPlot.py"])

    def perturbationStabiliy(self, epsilon):

        with open("Simulation Parameters/simulation.cfg") as file:
            lines = [line.rstrip() for line in file]

        a_0 = 0.6
        a_1 = 0.38
        K_w = 0.38
        K_l = 1.89
        A = 192.2
        B = 3.85
        dx = 20 / 2001
        s = 0.025  # default s value
        for line in lines:
            line = line.split(' ', 1)
            if line[0] == 's':
                s = float(line[1])

        j = 0
        val = 1.0
        while (val > 1e-10):
            val = 1.0 / (s * n.sqrt(n.pi)) * n.exp(- (j * dx / s) ** 2)
            j += 1
        j -= 1

        conv_weights = n.zeros(2 * j + 1)
        for i in range(j + 1):
            if not i:
                conv_weights[j] = 1.0 / (s * n.sqrt(n.pi))
            else:
                conv_weights[j + i] = 1.0 / (s * n.sqrt(n.pi)) * n.exp(- (i * dx / s) ** 2)
                conv_weights[j - i] = 1.0 / (s * n.sqrt(n.pi)) * n.exp(- (i * dx / s) ** 2)
        conv_weights /= n.sum(conv_weights)

        x = n.linspace(-10, 10, 2001)
        T = self.solution(x)

        B_vec = n.ones(2001) * B
        A_vec = n.where(T < -1, a_0, a_1)
        A_vec = n.where((x >= self.C_0) & (x <= self.C_1) & (T < 0), a_0, A_vec)
        A_vec = n.where((x >= self.C_0) & (x <= self.C_1) & (T >= 0), a_1, A_vec)

        padding = n.ones((conv_weights.shape[0] - 1) // 2) * a_0
        A_vec = n.concatenate((padding, n.convolve(A_vec, conv_weights, "valid"), padding))

        N_vec = n.zeros(2001)
        N_vec[1:-1] = (A_vec[2:] - A_vec[:-2]) / (2 * dx)
        N_vec[0] = N_vec[1]
        N_vec[-1] = N_vec[-2]
        N_vec = N_vec * self.Q_point * self.S(x)

        delta_epsilon = - B_vec * epsilon + epsilon * N_vec

        return delta_epsilon






def main():
    sol_name = 'G 5.0'
    if sol_name[0] == 'W':
        cnt = [0, 0]
        s_function = r"$e^{-\frac{|x|}{2}}$" if sol_name[1] == 'E' else r"$\frac{4}{\sqrt{10 \pi}}e^{-\frac{x^2}{10}}$"
    elif sol_name[0] == 'G':
        cnt = [-1, 1]
        s_function = rf"$\frac{{4}}{{\sqrt{{ {float(sol_name[2:]):.1f} \pi}}}} e^{{-\frac{{x^2}}{{ {float(sol_name[2:]):.1f} }}}}$"
    elif sol_name[0] == 'E':
        cnt = [-1 + float(sol_name[2:]), 1 + float(sol_name[2:])]
        s_function = r"$e^{-\frac{|x|}{2}}$"


    filenames = [f for f in listdir(sol_name) if isfile(join(sol_name, f))]

    for filename in filenames:
        with open(sol_name + '/' + filename, 'r+') as file:
            content = file.read()
            if len(content) == 0:
                continue
            content = content.replace('{', '').replace('}', '').replace('"', '').replace('*^', 'e')
            file.seek(0)
            file.write(content)
            file.truncate()

        newFamily = solFamily(sol_name, cnt, filename)
        solutionFamilies.append(newFamily)

    fig, ax = plt.subplots(1, 1)
    # fig.suptitle(f"Continent at [{cnt[0]}, {cnt[1]}] and S(x)={s_function}")
    fig.suptitle(f"Bifurcation diagram with S(x)={s_function}")
    ax.set_xlabel("Q")
    ax.set_ylabel("T(0)")

    point = plt.scatter([n.nan], [n.nan], s=5, facecolors='none', edgecolors='Black')

    for family in solutionFamilies:
        family.addToDisplay(ax, measure='T_0')

    xlim = ax.get_xlim()
    xL = xlim[1] - xlim[0]
    ylim = ax.get_ylim()
    yL = ylim[1] - ylim[0]

    def onclick(event):
        """
        Event handeler for clicking on a datapoint
        """
        if event.button == 3:
            # Finding family closest to point
            dists = []
            mousePos = n.array([event.xdata, event.ydata])
            mousePos[0] = (mousePos[0] - xlim[0]) / xL
            mousePos[1] = (mousePos[1] - ylim[0]) / yL
            for solutionFamily in solutionFamilies:
                points = solutionFamily.getSolCoordinates()
                points[:, 0] = (points[:, 0] - xlim[0]) / xL
                points[:, 1] = (points[:, 1] - ylim[0]) / yL
                ds = n.sum((points - mousePos) ** 2, axis=1)
                dists.append(n.min(ds))
            closestFamily = solutionFamilies[n.argmin(dists)]

            # Locating closest solution in family
            dists = []
            solutions = closestFamily.getSolCoordinates()
            solutions[:, 0] = (solutions[:, 0] - xlim[0]) / xL
            solutions[:, 1] = (solutions[:, 1] - ylim[0]) / yL
            for i in range(solutions.shape[0]):
                dists.append(n.sum((solutions[i] - mousePos) ** 2))

            # The closes solution to click if it's close enough
            closest = closestFamily.solutions[n.argmin(dists)]
            closest.runFiniteDifference()
        if event.button == 1:
            # Finding family closest to point
            dists = []
            mousePos = n.array([event.xdata, event.ydata])
            mousePos[0] = (mousePos[0] - xlim[0]) / xL
            mousePos[1] = (mousePos[1] - ylim[0]) / yL
            for solutionFamily in solutionFamilies:
                points = solutionFamily.getSolCoordinates()
                points[:, 0] = (points[:, 0] - xlim[0]) / xL
                points[:, 1] = (points[:, 1] - ylim[0]) / yL
                ds = n.sum((points - mousePos) ** 2, axis=1)
                dists.append(n.min(ds))
            closestFamily = solutionFamilies[n.argmin(dists)]

            # Locating closest solution in family
            dists = []
            solutions = closestFamily.getSolCoordinates()
            solutions[:, 0] = (solutions[:, 0] - xlim[0]) / xL
            solutions[:, 1] = (solutions[:, 1] - ylim[0]) / yL
            for i in range(solutions.shape[0]):
                dists.append(n.sum((solutions[i] - mousePos) ** 2))

            # The closes solution to click if it's close enough
            closest = closestFamily.solutions[n.argmin(dists)]
            pert_1 = closest.perturbationStabiliy(0.05)
            pert_2 = closest.perturbationStabiliy(-0.05)

            fig, ax = plt.subplots(1, 1)
            ax.plot(pert_1, label=r'$\epsilon=0.05$')
            ax.plot(pert_2, label=r'$\epsilon=-0.05$')
            ax.legend()
            plt.show()




    def onmove(event):
        """
        Event handeler for clicking on a datapoint
        """
        # Finding family closest to point
        try:
            dists = []
            mousePos = n.array([event.xdata, event.ydata])
            mousePos[0] = (mousePos[0] - xlim[0]) / xL
            mousePos[1] = (mousePos[1] - ylim[0]) / yL
            for solutionFamily in solutionFamilies:
                points = solutionFamily.getSolCoordinates()
                points[:, 0] = (points[:, 0] - xlim[0]) / xL
                points[:, 1] = (points[:, 1] - ylim[0]) / yL
                ds = n.sum((points - mousePos) ** 2, axis=1)
                dists.append(n.min(ds))
            closestFamily = solutionFamilies[n.argmin(dists)]

            # Locating closest solution in family
            dists = []
            solutions = closestFamily.getSolCoordinates()
            solutions[:, 0] = (solutions[:, 0] - xlim[0]) / xL
            solutions[:, 1] = (solutions[:, 1] - ylim[0]) / yL
            for i in range(solutions.shape[0]):
                dists.append(n.sum((solutions[i] - mousePos) ** 2))

            # The closes solution to click if it's close enough
            closest = closestFamily.solutions[n.argmin(dists)]
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


def createVideo():
    videoMaker = makeVideo("G 6.0")
    videoMaker.generate("Video G6.0")


if __name__ == '__main__':
    main()
