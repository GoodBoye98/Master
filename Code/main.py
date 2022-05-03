import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc
from main_package import solutionData, mergeFiles, sortFiles, processFiles, highlightSolution


def main():
    # rc('text', usetex=True)
    # rc('font', size=14)
    # rc('legend', fontsize=13)
    # rc('text.latex', preamble=r'\usepackage{cmbright}')

    sType = "Gaussian"
    dType = "Water"

    # If it is the first time running the files
    preprocessing = True
    familyNames = [
        # "1 - IceSnowIce", 
        # "2 - IceWaterSnowWaterIce", 
        # "3 - IceWaterSnowLandSnowWaterIce", 
        # "3 - IceWaterLandSnowLandWaterIce", 
        # "4 - IceWaterLandWaterIce", 
        # "2 - IceWaterSnowIce", 
        # "3 - IceWaterSnowWaterIce", 
        # "3A - IceWater_Snow_WaterIce", 
        # "4 - IceWaterSnowLandSnowWaterIce", 
        # "4# - IceWaterSnowLandSnowWaterIce", 
        # "4B - IceWaterSnowLandSnowIce", 
        # "4B# - IceWaterSnowLandSnowIce", 
        # "5A - IceWaterLandSnowWaterIce", 
        # "5B - IceWaterSnowLandWaterIce", 
        # "5C - IceWaterLandSnowIce", 
        # "5D - IceSnowLandWaterIce", 
        # "6 - IceWaterLandWaterIce"
        ]
    sortFrom = 'bot'

    if preprocessing:
        for family in familyNames:
            # processFiles(f"{sType}/{dType}", family)
            # mergeFiles(f"{sType}/{dType}", family)
            sortFiles(f"{sType}/{dType}", family, sortFrom)

    sData = solutionData(sType, dType)

    fig, ax = plt.subplots(1, 1, figsize=(8, 6))
    fig.canvas.manager.window.wm_geometry("+%d+%d" % (0, 0))

    highlightSolution(fig, ax, sData)

    sData.plot(fig, ax, "plot")
    ax.set_xlabel('Q')
    ax.set_ylabel('T(0)')

    # fig.suptitle('')
    # ax.scatter([485.152], [0.80328], s=20, facecolors='none', edgecolors='red')
    # ax.scatter([476.729], [0.8841], s=20, facecolors='none', edgecolors='black')
    # legend = ax.legend()
    # legend.remove()
    # fig.suptitle("Bifurcation diagram with " + r"$S(x)=e^{-\frac{|x|}{2}}$")
    # fig.suptitle("Bifurcation diagram with " + r"$S(x)=\frac{4}{\sqrt{10 \pi}} e^{-\frac{x^2}{10}}$")
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys())
    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    main()