import string
import subprocess
import numpy as n
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from os import listdir, remove, mkdir, rmdir
from os.path import isfile, isdir, join


def consecutive(data, stepsize=1):
    return n.split(data, n.where(n.diff(data) != stepsize)[0] + 1)


class makeVideo:
    def __init__(self, dir):
        files = [f for f in listdir(dir) if isfile(join(dir, f))]
        data = [pd.read_csv(f"{dir}/{file}", header=None).values for file in files]
        order = n.zeros(len(files))

        for i in range(len(files)):
            if "_" in files[i]:
                num = int(''.join([s for s in files[i][1:] if s in string.digits]))
                if "_u_" in files[i]:
                    order[i] = 10000 + num
                    data[i] = data[i][n.argsort(data[i][:, 0])]
                    if num % 2:
                        data[i] = data[i][::-1]
                else:
                    order[i] = 20000 - num
                    data[i] = data[i][n.argsort(data[i][:, 0])]
                    if num % 2:
                        data[i] = data[i][::-1]
            elif "IceSnowIce" in files[i]:
                order[i] = 0
                data[i] = data[i][n.argsort(data[i][:, 0])]
            elif "IceWaterSnowWaterIce" in files[i]:
                order[i] = 1
                data[i] = data[i][n.argsort(data[i][:, 1])]
            elif "IceWaterLandWaterIce" in files[i]:
                order[i] = 30000
                data[i] = data[i][n.argsort(data[i][:, 0])]

        order = n.argsort(order)
        newData = []
        for i in order:
            newData.append(data[i])

        self.data = n.vstack(newData)
        # Load and normalize q-values
        self.Q_stored = self.data[:, 0].copy()
        self.Q = self.data[:, 0]
        self.Q -= n.min(self.Q)
        self.Q /= n.max(self.Q)
        # Load and normalize t_values
        self.T_0 = self.data[:, 1]
        self.T_0 -= n.min(self.T_0)
        self.T_0 /= n.max(self.T_0)
        # Load temp profiles
        self.T = self.data[:, 4:]

        # Create list of normalized points for solutions
        self.solPoints = n.column_stack((self.Q, self.T_0))


    def generate(self, savePath):
        prev_p = n.array([0, 0])
        x = n.linspace(-10, 10, 2001)
        idx = n.argwhere(n.abs(x) < 5)
        l_0, l_1 = n.argwhere(x[idx] == -1)[0][0], n.argwhere(x[idx] == 1)[0][0]
        if not isdir(savePath):
            mkdir(savePath)
        for i in range(self.solPoints.shape[0]):
            plt.cla()
            plt.title(f"Q={self.Q_stored[i]:.2f}")
            plt.xlabel("x")
            plt.ylabel("T(0)")

            ice = n.argwhere(self.T[i][idx] < -1)[:, 0]
            ice = ice[(ice < l_0) | (ice > l_1)]
            water = n.argwhere(self.T[i][idx] >= -1)[:, 0]
            water = water[(water < l_0) | (water > l_1)]

            snow = n.argwhere(self.T[i][idx] < 0)[:, 0]
            snow = snow[(snow >= l_0) & (snow <= l_1)]
            land = n.argwhere(self.T[i][idx] >= 0)[:, 0]
            land = land[(land >= l_0) & (land <= l_1)]

            if ice.shape[0]:
                plt.plot(0, 0, c="lightblue", path_effects=[pe.Stroke(linewidth=3, foreground='black'), pe.Normal()], label='Ice')
                segments = consecutive(ice)
                for seg in segments:
                    plt.plot(x[idx[seg]], self.T[i][idx[seg]], c="lightblue", path_effects=[pe.Stroke(linewidth=3, foreground='black'), pe.Normal()])

            if water.shape[0]:
                plt.plot(0, 0, c="blue", path_effects=[pe.Stroke(linewidth=3, foreground='black'), pe.Normal()], label='Water')
                segments = consecutive(water)
                for seg in segments:
                    plt.plot(x[idx[seg]], self.T[i][idx[seg]], c="blue", path_effects=[pe.Stroke(linewidth=3, foreground='black'), pe.Normal()])

            if snow.shape[0]:
                plt.plot(0, 0, c="snow", path_effects=[pe.Stroke(linewidth=3, foreground='black'), pe.Normal()], label='Snow')
                segments = consecutive(snow)
                for seg in segments:
                    plt.plot(x[idx[seg]], self.T[i][idx[seg]], c="snow", path_effects=[pe.Stroke(linewidth=3, foreground='black'), pe.Normal()])

            if land.shape[0]:
                plt.plot(0, 0, c="brown", path_effects=[pe.Stroke(linewidth=3, foreground='black'), pe.Normal()], label='Land')
                segments = consecutive(land)
                for seg in segments:
                    plt.plot(x[idx[seg]], self.T[i][idx[seg]], c="brown", path_effects=[pe.Stroke(linewidth=3, foreground='black'), pe.Normal()])

            plt.legend()
            plt.savefig(savePath + "/file%02d.png" % i)
            print(f"{100 * i / self.data.shape[0]:.2f}%")

        subprocess.run([
            'ffmpeg', '-framerate', '60', '-i', 'file%02d.png', '-r', '30', '-pix_fmt', 'yuv420p',
            'video_name.mp4'], cwd=savePath)

        files = [f for f in listdir(savePath) if isfile(join(savePath, f))]
        for file in files:
            if file != "video_name.mp4":
                remove(f"{savePath}/{file}")
