import numpy as n
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from main_package import solutionData, colorStationarySolution

def main():

    sType = "Exponential"
    dType = "Water"

    sData = solutionData(sType, dType)

    y = sData.getSol(150)
    x = n.linspace(-6, 6, y.shape[0])

    fig, ax = plt.subplots(1, 1, figsize=(6, 2))

    colorStationarySolution(fig, ax, 10 * n.concatenate(([0], [0], [0], [0], y)))
    ax.hlines(-10.3, -7, 7, color=(0, 0, 0))
    ax.fill_between([-7, 7], -10, -55, color='lightskyblue', alpha=0.1)
    ax.fill_between([-7, 7], 20, -10, color='blue', alpha=0.1)
    ax.text(-6, -7, r'$a(T)=a_1$')
    ax.text(-6, -25, r'$a(T)=a_0$')

    ax.set_xlim((-6.5, 6.5))
    ax.set_ylim((-51, 0))

    ax.text(0, -26, r'$T$', ha='center')
    ax.set_xticks([0])
    ax.set_xticklabels(["0"])

    plt.yticks((0, -10, -20, -30, -40), ('0', '-1', '-2', '-3', '-4'))
    plt.ylabel(r'$T(x)$')
    plt.xlabel('x')
    plt.tight_layout()
    plt.show()


def main_2():

    sType = "Exponential"
    dType = "Water"

    sData = solutionData(sType, dType)

    y = sData.getSol(570)
    x = n.linspace(-6, 6, y.shape[0])

    fig, ax = plt.subplots(1, 1, figsize=(6, 2))

    colorStationarySolution(fig, ax, 10 * n.concatenate(([0], [0], [0], [0], y)))
    ax.vlines((-0.797, 0.797), -48, 16, color=(0, 0, 0, 0.2), ls=':')
    ax.hlines(-10.3, -7, 7, color=(0, 0, 0))
    ax.fill_between([-7, 7], -10, -55, color='lightskyblue', alpha=0.1)
    ax.fill_between([-7, 7], 20, -10, color='blue', alpha=0.1)
    ax.text(-6, -0, r'$a(T)=a_1$')
    ax.text(-6, -25, r'$a(T)=a_0$')

    ax.text(-2.5, -46, r'$T_1$', ha='center')
    ax.text(-0, -2, r'$T_2$', ha='center')
    ax.text( 2.5, -46, r'$T_3$', ha='center')

    ax.set_xlim((-6.5, 6.5))
    ax.set_ylim((-51, 20))
    
    plt.xticks((-0.797, 0.797), (r'$A_0$', r'$A_1$'))
    plt.yticks((20, 10, 0, -10, -20, -30, -40), ('2', '1', '0', '-1', '-2', '-3', '-4'))
    plt.ylabel(r'$T(x)$')
    plt.xlabel('x')
    plt.tight_layout()
    plt.show()


def main_3():

    sType = "Exponential"
    dType = "Continent Offset 0.000"

    sData = solutionData(sType, dType)

    y = sData.getSol(250)
    x = n.linspace(-6, 6, y.shape[0])

    fig, ax = plt.subplots(1, 1, figsize=(6, 2))

    sep_x = [-7, -1, -1, 1, 1, 7]
    sep_y = [-10.3, -10.3, 0, 0, -10.3, -10.3]
    top_y = [20] * len(sep_x)
    bot_y = [-55] * len(sep_x)

    colorStationarySolution(fig, ax, 10 * n.concatenate(([0], [0], [-0.1], [0.1], y)))
    ax.plot(sep_x, sep_y, color=(0, 0, 0))
    ax.fill_between(sep_x, sep_y, bot_y, color='lightskyblue', alpha=0.1)
    ax.fill_between(sep_x, sep_y, top_y, color='blue', alpha=0.1)
    ax.vlines([-1, 1], -55, 20, color='black', alpha=0.1, ls=':')
    ax.text(-6, -3, r'$a(T)=a_1$')
    ax.text(-6, -25, r'$a(T)=a_0$')

    ax.text(-2.5, -46, r'$T_1$', ha='center')
    ax.text(-0, -20, r'$T_2$', ha='center')
    ax.text( 2.5, -46, r'$T_3$', ha='center')

    ax.set_xlim((-6.5, 6.5))
    ax.set_ylim((-51, 7))
    
    plt.xticks((-1, 1), ('-1', '1'))
    plt.yticks((0, -10, -20, -30, -40, -50), ('0', '-1', '-2', '-3', '-4', '-5'))
    plt.ylabel(r'$T(x)$')
    plt.xlabel('x')
    plt.tight_layout()
    plt.show()


def main_4():

    sType = "Gaussian"
    dType = "Std 6"

    sData = solutionData(sType, dType)

    y = sData.getSol(350, specialFamily="3A - IceWater_Land_WaterIce")
    x = n.linspace(-6, 6, y.shape[0])

    fig, ax = plt.subplots(1, 1, figsize=(6, 2))

    sep_x = [-7, -1, -1, 1, 1, 7]
    sep_y = [-10.3, -10.3, 0, 0, -10.3, -10.3]
    top_y = [20] * len(sep_x)
    bot_y = [-55] * len(sep_x)

    colorStationarySolution(fig, ax, 10 * n.concatenate(([0], [0], [-0.1], [0.1], y)))
    ax.plot(sep_x, sep_y, color=(0, 0, 0))
    ax.fill_between(sep_x, sep_y, bot_y, color='lightskyblue', alpha=0.1)
    ax.fill_between(sep_x, sep_y, top_y, color='blue', alpha=0.1)
    ax.vlines([-1, 1], -55, 20, color='black', alpha=0.1, ls=':')
    ax.text(-2.4, -8, r'$a(T)=a_1$')
    ax.text(-2.4, -16, r'$a(T)=a_0$')

    ax.text(-0.91, 0.04, r'$G_{-N}$', ha='center')
    ax.text(-0.87, -0.06, r'$G_{-N+1}$', ha='center')
    ax.text(-0.15, -0.09, r'$G_{-1}$', ha='center')
    ax.text( 0.0, 0.03, r'$G_0$', ha='center')
    ax.text( 0.15, -0.09, r'$G_1$', ha='center')
    ax.text( 0.87, -0.06, r'$G_{N-1}$', ha='center')
    ax.text( 0.91, 0.04, r'$G_{N}$', ha='center')

    k = n.linspace(-1, 1, 15)
    for i in range(15):
        ax.plot([k[i], 0], [-0.4, -10], color='black', alpha=0.6, lw=0.75)

    ax.set_xlim((-1.1, 1.1))
    ax.set_ylim((-0.1, 0.1))
    
    plt.xticks((-1, 0, 1), ('-1', '0', '1'))
    plt.yticks((0.1, 0.05, 0, -0.05, -0.1), ('0.01', '0.005', '0', '-0.005', '-0.1'))
    plt.ylabel(r'$T(x)$')
    plt.xlabel('x')
    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    # main()
    # main_2()
    # main_3()
    main_4()