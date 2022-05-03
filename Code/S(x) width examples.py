import numpy as n
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from main_package import solutionData, colorStationarySolution

def main():
    eData = solutionData("Exponential", "Water")
    gData = solutionData("Gaussian", "Water")

    # Loading exponential data
    y_exp = eData.getSol(120)
    x_exp = n.linspace(-6, 6, y_exp.shape[0])

    # Loading gaussian data
    y_gss = gData.getSol(0)
    x_gss = n.linspace(-6, 6, y_gss.shape[0])

    fig, ax = plt.subplots(1, 1, figsize=(6, 2))
    # fig.suptitle('Q=400', x=0.54)

    ax.plot(x_exp, y_exp, label='Laplacian', zorder=0)
    ax.plot(x_gss, y_gss, label='Gaussian', zorder=0)

    max_exp = n.max(y_exp)
    max_gss = n.max(y_gss)
    ax.hlines([max_exp, max_gss], -3, 0, color='black', ls=':', alpha=0.9)
    ax.text(-3.8, max_exp- 0.15, f"T={max_exp:.2f}", ha='center', size='small')
    ax.text(-3.8, max_gss- 0.15, f"T={max_gss:.2f}", ha='center', size='small')

    ax.scatter([0, 0], [max_exp, max_gss], c='black', s=5, marker='x', zorder=1)

    ax.legend()
    ax.set_ylim((-5.2, -0.5))

    ax.set_ylabel(r'$T(x)$')
    ax.set_xlabel('x')

    # ax_2 = ax.twinx()
    # x = n.linspace(-6, 6, 200)
    # exp = n.exp(- n.abs(x) / 2)
    # gss = 4 / n.sqrt(10 * n.pi) * n.exp(- x ** 2 / 10)
    # ax_2.plot(x, exp, color='black', ls='--', alpha=0.5)
    # ax_2.plot(x, gss, color='black', ls='--', alpha=0.5)
    # ax_2.set_ylabel('S(x)')

    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    main()