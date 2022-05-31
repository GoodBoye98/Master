import matplotlib
import matplotlib.pyplot as plt
from main_package import solutionData, mergeFiles, sortFiles, processFiles, invertFiles, highlightSolution


def main():

    sType = "Gaussian"
    dType = "Std 6"

    # If it is the first time running the files
    preprocessing = True
    familyNames = [
        "1 - IceSnowIce", 
        # "2 - IceWaterSnowWaterIce", 
        # "3 - IceWaterSnowLandSnowWaterIce", 
        # "3 - IceWaterLandSnowLandWaterIce", 
        # "4 - IceWaterLandWaterIce", 
        # "2 - IceWaterSnowIce", 
        # "3 - IceWaterSnowWaterIce", 
        # "3A - IceWater_Land_WaterIce", 
        # "4 - IceWaterSnowLandSnowWaterIce", 
        # "4# - IceWaterSnowLandSnowWaterIce", 
        # "4B - IceWaterSnowLandSnowIce", 
        # "4B# - IceWaterSnowLandSnowIce", 
        # "5A - IceWaterLandSnowWaterIce", 
        # "5B - IceWaterSnowLandWaterIce", 
        # "5C - IceWaterLandSnowIce", 
        # "5D - IceSnowLandWaterIce", 
        # "6 - IceWaterLandWaterIce"
        ]
    sortFrom = 'right'

    if preprocessing:
        for family in familyNames:
            # processFiles(f"{sType}/{dType}", family)
            # mergeFiles(f"{sType}/{dType}", family)
            # sortFiles(f"{sType}/{dType}", family, sortFrom)
            invertFiles(f"{sType}/{dType}", family)

if __name__ == '__main__':
    main()