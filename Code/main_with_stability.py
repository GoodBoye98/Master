import numpy as n
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib.gridspec import GridSpec
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm
from main_package import solutionData, mergeFiles, sortFiles, processFiles, highlightSolution


def main():

    sType = "Exponential"
    dType = "Continent Offset 0.030"

    sData = solutionData(sType, dType)

    fig = plt.figure(figsize=(8, 8))

    gs = GridSpec(3, 2, figure=fig)
    ax1 = fig.add_subplot(gs[0:2, :])
    ax2 = fig.add_subplot(gs[2, 0])
    ax3 = fig.add_subplot(gs[2, 1])

    sData.plot(fig, ax1, "plot")
    ax1.set_xlabel('Q')
    ax1.set_ylabel('T(0)')
    
    handles, labels = ax1.get_legend_handles_labels()
    by_label = dict(sorted(dict(zip(labels, handles)).items()))
    ax1.legend(by_label.values(), by_label.keys())

    data = n.load(f'Stability/{sType}-{dType}-stability-N4001.npy')

    sData = solutionData(sType, dType)
    sols = sData.getSols()

    idx = data[:, 0]
    pert = data[:, 1]
    error = data[:, 2]

    idx_u = idx[pert > 0].astype(int)
    idx_d = idx[pert < 0].astype(int)

    error_u = error[pert > 0]
    error_d = error[pert < 0]

    set_lims = False
    size = 0.1
    x_lims = (470, 495)
    y_lims = (0.6, 1.2)
    if set_lims:
        size = 1.0
        ax1.set_xlim(x_lims)
        ax2.set_xlim(x_lims)
        ax3.set_xlim(x_lims)

        ax1.set_ylim(y_lims)
        ax2.set_ylim(y_lims)
        ax3.set_ylim(y_lims)

    scatter_0 = ax2.scatter(sols[idx_u, 0], sols[idx_u, 1], 
        s=size, c=error_u[idx_u], cmap='rainbow', vmin=0, vmax=0.2)
    ax2.set_title('Upward perturbation')
    ax2.set_xlabel('Q')
    ax2.set_ylabel('T(0)')

    bar = fig.colorbar(scatter_0, ax=ax2)
    bar.ax.set_title('Error')

    scatter_1 = ax3.scatter(sols[idx_d, 0], 
        sols[idx_d, 1], s=size, c=error_d[idx_d], cmap='rainbow', vmin=0, vmax=0.2)
    ax3.set_title('Downward perturbation')
    ax3.set_xlabel('Q')
    ax3.set_ylabel('T(0)')

    bar = fig.colorbar(scatter_1, ax=ax3)
    bar.ax.set_title('Error')

    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    main()
