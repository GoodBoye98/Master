import numpy as n
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from main_package import solutionData, colorStationarySolution
from matplotlib.gridspec import GridSpec


def main():
    eData = solutionData("Exponential", "Continent Offset 0.000")
    gData = solutionData("Gaussian", "Std 6")

    # Loading exponential data
    y_exp = eData.getSol(216, specialFamily="3 - IceWaterSnowLandSnowWaterIce")
    x_exp = n.linspace(-6, 6, y_exp.shape[0])

    # Loading gaussian data
    y_gss = gData.getSol(0, specialFamily="3B - IceWater_Snow_WaterIce")
    x_gss = n.linspace(-6, 6, y_gss.shape[0])

    fig = plt.figure(constrained_layout=True, figsize=(6, 6))

    gs = GridSpec(3, 2, figure=fig)
    ax1 = fig.add_subplot(gs[0:2, :])
    ax2 = fig.add_subplot(gs[2, 0])
    ax3 = fig.add_subplot(gs[2, 1])

    ax1.plot(x_exp, y_exp, label='Laplacian', zorder=0)
    ax1.plot(x_gss, y_gss, label='Gaussian', zorder=0)

    ax1.hlines(0, -3, 3, color='black', ls='--', alpha=0.4)
    # ax1.scatter([-1, 1], [0, 0], facecolors='none', s=20, edgecolors='r')

    ax1.legend()
    ax1.set_ylim((-3.7, 0.5))
    ax1.set_xlim((-2.7, 2.7))

    ax1.set_ylabel(r'$T(x)$')
    ax1.set_xlabel('x')

    eData.plot(fig, ax2,'plot', linewidth=1.0)
    gData.plot(fig, ax3,'plot', linewidth=1.0)
    ax2.get_legend().remove()
    ax2.set_title('Laplacian')
    ax2.scatter([595.6000134571465], [0], s=20, facecolor='none', edgecolor='r')
    ax2.set_ylabel(r'$T(0)$')
    ax2.set_xlabel('Q')

    ax3.get_legend().remove()
    ax3.set_title('Gaussian')
    ax3.scatter([522.61928], [-0.060929], s=20, facecolor='none', edgecolor='r')
    ax3.set_ylabel(r'$T(0)$')
    ax3.set_xlabel('Q')

    fig.suptitle('')

    plt.show()

    # # ax_2 = ax.twinx()
    # # x = n.linspace(-6, 6, 200)
    # # exp = n.exp(- n.abs(x) / 2)
    # # gss = 4 / n.sqrt(10 * n.pi) * n.exp(- x ** 2 / 10)
    # # ax_2.plot(x, exp, color='black', ls='--', alpha=0.5)
    # # ax_2.plot(x, gss, color='black', ls='--', alpha=0.5)
    # # ax_2.set_ylabel('S(x)')

    # plt.tight_layout()
    # plt.show()

if __name__ == '__main__':
    main()