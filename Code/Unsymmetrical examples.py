import numpy as n
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from main_package import solutionData, colorStationarySolution

def main():
    eData = solutionData("Exponential", "Continent Offset 0.000")

    # Loading exponential data
    y_1 = eData.getSol(19, specialFamily="5A - IceWaterLandSnowWaterIce")
    x_1 = n.linspace(-6, 6, y_1.shape[0])

    # Loading gaussian data
    y_2 = eData.getSol(19, specialFamily="5B - IceWaterSnowLandWaterIce")
    x_2 = n.linspace(-6, 6, y_2.shape[0])

    fig, ax = plt.subplots(1, 1, figsize=(6, 3))

    ax.plot(x_1, y_1, label=r'$T_{5}(x)$',  zorder=1, c='coral')
    ax.plot(x_2, y_2, label=r'$T_{5C}(x)$', zorder=1, c='dodgerblue')

    ax.plot(x_1, y_1[::-1], label=r'$T_{5}(-x)$',  zorder=2, ls=(0, (1, 4)), c='coral')
    ax.plot(x_2, y_2[::-1], label=r'$T_{5C}(-x)$', zorder=2, ls=(0, (1, 4)), c='dodgerblue')

    ax.fill_between(x_1, y_1, y_2, color='grey', alpha=0.5, zorder=0)

    ax.legend()
    ax.set_ylabel(r'$T(x)$')
    ax.set_xlabel('x')
    ax.set_xlim((-2, 2))
    ax.set_ylim((-3.2, 1.2))

    plt.grid(True)
    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    main()