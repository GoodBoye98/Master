import matplotlib.pyplot as plt
from main_package import solutionData, mergeFiles, sortFiles, processFiles, highlightSolution


def main():
    sType = "Gaussian"
    dType = "Water"

    # If it is the first time running the files
    preprocessing = False
    familyNames = ["2 - IceWaterIce"]
    sortStart = 'bot'

    if preprocessing:
        for family in familyNames:
            processFiles(f"{sType}/{dType}", family)
            mergeFiles(f"{sType}/{dType}", family)
            sortFiles(f"{sType}/{dType}", family, sortStart)

    sData = solutionData(sType, dType)

    fig, ax = plt.subplots(1, 1, figsize=(8, 6))

    highlightSolution(fig, ax, sData)

    sData.plot(fig, ax, "plot")
    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    main()