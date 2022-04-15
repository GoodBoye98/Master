import matplotlib.pyplot as plt
from main_package import solutionData, mergeFiles, sortFiles, processFiles, highlightSolution


def main():
    sType = "Exponential"
    dType = "Continent Offset 0.000"

    # If it is the first time running the files
    preprocessing = True
    familyNames = [
        # "1 - IceSnowIce", 
        # "2 - IceWaterSnowIce", 
        # "3 - IceWaterSnowWaterIce", 
        # "4 - IceWaterSnowLandSnowWaterIce", 
        # "4B - IceWaterSnowLandSnowIce", 
        # "5 - IceWaterLandSnowWaterIce", 
        # "5B - IceWaterLandSnowIce", 
        # "5C - IceWaterSnowLandWaterIce", 
        # "6 - IceWaterLandWaterIce"
        ]
    sortFrom = 'right'

    if preprocessing:
        for family in familyNames:
            processFiles(f"{sType}/{dType}", family)
            mergeFiles(f"{sType}/{dType}", family)
            sortFiles(f"{sType}/{dType}", family, sortFrom)

    sData = solutionData(sType, dType)

    fig, ax = plt.subplots(1, 2, figsize=(10, 5))

    highlightSolution(fig, ax[0], sData)

    sData.plot(fig, ax[0], "scatter")
    sData.plot(fig, ax[1], "scatter")
    ax[0].set_xlabel('Q')
    ax[0].set_ylabel('T(0)')
    ax[1].set_xlabel('Q')
    ax[1].set_ylabel('T(0)')
    ax[1].set_xlim((465, 515))
    ax[1].set_ylim((0.3, 1.4))

    fig.suptitle("Bifurcation diagram. Continent at (-0.950, 1.050)")
    # ax[0].get_legend().remove()
    # ax[1].get_legend().remove()
    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    main()