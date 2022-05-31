import matplotlib.pyplot as plt
from main_package import solutionData, highlightSolution


def main():
    sType = "Gaussian"
    dType = "Std 6"

    sData = solutionData(sType, dType)

    fig, ax = plt.subplots(1, 1, figsize=(8, 6))

    highlightSolution(fig, ax, sData)

    sData.plot(fig, ax, "plot")
    ax.set_xlabel('Q')
    ax.set_ylabel('T(0)')

    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys())
    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    main()