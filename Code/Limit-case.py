import numpy as n
import matplotlib.pyplot as plt
import matplotlib.ticker as tkr     # has classes for tick-locating and -formatting
from main_package import solutionData, colorStationarySolution


def main():
    eData = n.genfromtxt("Gaussian/Std 6/LimitingCase.csv", delimiter=',')
    print(eData)

    # Loading exponential data
    y_1 = eData[1:]
    x_1 = n.linspace(-6, 6, y_1.shape[0])

    fig, ax = plt.subplots(1, 1, figsize=(6, 4))

    # colorStationarySolution(fig, ax, n.concatenate(([0], [0], [-1], [1], 10 * y_1)))
    idx = n.argwhere(n.abs(x_1) < 3.8)
    ax.plot(x_1[idx], y_1[idx],  zorder=1, c='coral')

    ax.legend()
    ax.set_ylabel(r'$T(x)$')
    ax.set_xlabel('x')
    # ax.set_xlim((-4, 4))
    ax.set_ylim((-5.0, 1.0))

    # plt.grid(True)
    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    main()