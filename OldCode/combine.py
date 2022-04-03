import sys
import numpy as n
import pandas as pd
import matplotlib.pyplot as plt
from plot import solFamily, solution
from os import listdir, remove, mkdir, rmdir
from os.path import isfile, isdir, join

if __name__ == '__main__':

    if len(sys.argv) == 3:
        dir = sys.argv[1]
        filename = sys.argv[2]

        if not isdir(dir):
            raise Exception(f"Directory {dir} does not exist")

        files = [f for f in listdir(dir) if isfile(join(dir, f))]

        files_correct = [f for f in files if filename in f]

        if len(files_correct) <= 1:
            raise Exception(f"No extra files matching '{filename}' was found")

        print("Files found:")
        for file in files_correct:
            print(f"\t{file}")

        all_string = ""

        for file in files_correct:
            if filename in file:
                with open(f"{dir}/{file}", "r") as f:
                    all_string += f.read()
                remove(f"{dir}/{file}")

        with open(f"{dir}/{filename}.csv", "w") as f:
            f.write(all_string)
    else:
        print("Invalid arguments")
