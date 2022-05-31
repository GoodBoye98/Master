import string
import subprocess
import numpy as n
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from main_package import solutionData, colorStationarySolution
from tqdm import tqdm
from os import listdir, remove, mkdir, rmdir
from os.path import isfile, isdir, join


def consecutive(data, stepsize=1):
    return n.split(data, n.where(n.diff(data) != stepsize)[0] + 1)


class makeVideo:
    def __init__(self, mainDir, subDir):

        self.solData = solutionData(mainDir, subDir)

        IceSnowIce                      = self.solData.getSolutionFamily('1 - IceSnowIce')
        IceWaterSnowWaterIce            = self.solData.getSolutionFamily('2 - IceWaterSnowWaterIce')
        IceWater_Snow_WaterIce          = self.solData.getSolutionFamily('3B - IceWater_Snow_WaterIce')
        IceWater_Land_WaterIce          = self.solData.getSolutionFamily('3A - IceWater_Land_WaterIce')
        IceWaterLandWaterIce            = self.solData.getSolutionFamily('4 - IceWaterLandWaterIce')


        self.data = {
            '1 - IceSnowIce': IceSnowIce,
            '2 - IceWaterSnowWaterIce': IceWaterSnowWaterIce,
            '3B - IceWater_Snow_WaterIce': IceWater_Snow_WaterIce,
            '3A - IceWater_Land_WaterIce': IceWater_Land_WaterIce,
            '4 - IceWaterLandWaterIce': IceWaterLandWaterIce
        }


    def generate(self, savePath, saveName):
        if not isdir(savePath):
            mkdir(savePath)

        fig, ax = plt.subplots(1, 2, figsize=(10, 5))

        # Find limits and stuff
        T = self.solData.getSol(list(self.data.values())[-1].shape[0] - 1, specialFamily=list(self.data.keys())[-1])
        self.solData.plot(fig, ax[0])
        colorStationarySolution(fig, ax[1], n.concatenate((n.array([0, 0, -1, 1]), T)), scale=10)
        ax0_xlim = ax[0].get_xlim()
        ax0_ylim = ax[0].get_ylim()
        ax1_xlim = ax[1].get_xlim()
        ax1_ylim = ax[1].get_ylim()
        plt.tight_layout()


        # Start making images
        j = 0
        N = n.sum([k.shape[0] for k in self.data.values()])
        with tqdm(total = N) as pb:

            for family, data in self.data.items():

                for i in range(data.shape[0]):

                    ax[0].cla()
                    ax[1].cla()

                    # Load solution data
                    T = self.solData.getSol(i, specialFamily=family)[::-1]

                    # Plot bifurcation diagram
                    self.solData.plot(fig, ax[0])
                    ax[0].scatter(*data[i], s=5, facecolors='none', edgecolors='Black')
                    handles, labels = ax[0].get_legend_handles_labels()
                    by_label = dict(zip(labels, handles))
                    ax[0].legend(by_label.values(), by_label.keys(), prop={'size': 8})
                    ax[0].set_xlim(ax0_xlim)
                    ax[0].set_ylim(ax0_ylim)

                    # Plot solution
                    colorStationarySolution(fig, ax[1], n.concatenate((n.array([0, 0, -1, 1]), T)), scale=10)
                    ax[1].legend()
                    ax[1].set_xlim(ax1_xlim)
                    ax[1].set_ylim(ax1_ylim)

                    # Save figure
                    plt.savefig(savePath + "/file%02d.png" % j)
                    j += 1
                    pb.update(1)

        subprocess.run([
            'ffmpeg', '-framerate', '60', '-i', 'file%02d.png', '-r', '30', '-pix_fmt', 'yuv420p',
            f'{saveName}.mp4'], cwd=savePath)

        files = [f for f in listdir(savePath) if isfile(join(savePath, f))]
        for file in files:
            if '.png' in file:
                remove(f"{savePath}/{file}")


def main():
    genVid = makeVideo('Gaussian', 'Std 6')
    genVid.generate('img_dir', 'video')


if __name__ == '__main__':
    main()