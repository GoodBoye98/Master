import numpy as n
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from main_package import solutionData, colorStationarySolution

def main():
    
    data = pd.read_csv("Gaussian/LimitingCase.csv", header=None).values[0]

    y = data[2:]
    x = n.linspace(-6, 6, y.shape[0])

    fig, ax = plt.subplots(1, 1, figsize=(6, 2))

    sep_x = [-7, -1, -1, 1, 1, 7]
    sep_y = [-10.3, -10.3, 0, 0, -10.3, -10.3]
    top_y = [20] * len(sep_x)
    bot_y = [-55] * len(sep_x)

    colorStationarySolution(fig, ax, 10 * n.concatenate(([0], [0], [-0.1], [0.1], y)))
    # ax.plot(sep_x, sep_y, color=(0, 0, 0))
    # ax.fill_between(sep_x, sep_y, bot_y, color='lightskyblue', alpha=0.1)
    # ax.fill_between(sep_x, sep_y, top_y, color='blue', alpha=0.1)
    ax.vlines([-1, 1], -55, 20, color='black', alpha=0.1, ls=':')
    # ax.text(-6, -3, r'$a(T)=a_1$')
    # ax.text(-6, -25, r'$a(T)=a_0$')

    # ax.text(-2.5, -46, r'$T_1$', ha='center')
    # ax.text(-0, -20, r'$T_2$', ha='center')
    # ax.text( 2.5, -46, r'$T_3$', ha='center')

    ax.set_xlim((-6.5, 6.5))
    ax.set_ylim((-53, 4))
    
    plt.xticks((-6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6), ('6', '5', '4', '3', '2', '1', '0', '1', '2', '3', '4', '5', '6'))
    plt.yticks((0, -10, -20, -30, -40, -50), ('0', '-1', '-2', '-3', '-4', '-5'))
    plt.ylabel(r'$T(x)$')
    plt.xlabel(r'$x$')
    plt.tight_layout()
    plt.legend()
    plt.show()

if __name__ == '__main__':
    main()