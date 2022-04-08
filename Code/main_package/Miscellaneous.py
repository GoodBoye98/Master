import random
import string
import numpy as n
import pandas as pd
from os import listdir, remove
from os.path import isfile, join
from argparse import ArgumentError


def processFiles(dir, base_name):
    files = [f for f in listdir(dir) if isfile(join(dir, f))]
    correctFiles = [f for f in files if base_name in f]

    for file in correctFiles:
        with open(f"{dir}/{file}", "r+") as f:
            s = f.read()
            if len(s) == 0:
                continue
            s = s.replace('"', '').replace('{', '').replace('}', '').replace('*^', 'e')
            f.seek(0)
            f.write(s)
            f.truncate()


def mergeFiles(dir, base_name):
    files = [f for f in listdir(dir) if isfile(join(dir, f))]
    correctFiles = [f for f in files if base_name in f]

    if len(correctFiles) <= 1:
        print(f"Warning: No files matching name: '{base_name}' were found in directory: '{dir}'")
        return

    allString = ""
    for file in correctFiles:
        with open(f"{dir}/{file}", "r") as f:
            allString += f.read()
        
    for file in correctFiles:
        remove(f"{dir}/{file}")

    with open(f"{dir}/{base_name}.csv", "w") as f:
        f.write(allString)
    
    
    print(f"Sucessfully merged files {correctFiles} into '{base_name}.csv' and removed lefovers.")


def sortFiles(dir, base_name, how='bot', usecols=[0, 1]):
    files = [f for f in listdir(dir) if isfile(join(dir, f))]
    correctFiles = [f for f in files if base_name in f]

    for file in correctFiles:
        # Full data matrix
        dataFull = pd.read_csv(f"{dir}/{file}", header=None, delimiter=',').values

        if dataFull.dtype != float:
            print(f"File {file} can not be sorted as it has dtype '{dataFull.dtype}' and not float.\n" + 
                  f"Consider running processFiles() first.")
            continue

        # Data to use in sorting
        data = dataFull[:, usecols]
        dataCopy = data.copy()

        # What element is first in sorting
        if how == 'bot':
            newIdx = n.argmin(data[:, 1])
        elif how == 'top':
            newIdx = n.argmax(data[:, 1])
        elif how == 'left':
            newIdx = n.argmin(data[:, 0])
        elif how == 'right':
            newIdx = n.argmax(data[:, 0])
        else:
            raise ArgumentError(f"how='{how}' is not a valid argument. Try 'top', 'bot', 'right' or 'left'.")

        scaling = n.array([n.max(data[:, 0]) - n.min(data[:, 0]),
                           n.max(data[:, 1]) - n.min(data[:, 1])]) ** 2

        idxSorted = [newIdx]
        point = data[newIdx]
        data = n.delete(data, newIdx, 0)
        while data.shape[0]:
            # Find next point in diminished array
            d = (data[:, 0] - point[0]) ** 2 / scaling[0] + (data[:, 1] - point[1]) ** 2 / scaling[1]
            newIdx = n.argmin(d)
            point = data[newIdx].copy()
            data = n.delete(data, newIdx, 0)

            # Find index in full array and appending to idxSorted
            d = (dataCopy[:, 0] - point[0]) ** 2 / scaling[0] + (dataCopy[:, 1] - point[1]) ** 2 / scaling[1]
            newIdx = n.argmin(d)
            idxSorted.append(newIdx)
        
        # Array of indexes to sort after
        idxSorted = n.array(idxSorted)
        n.savetxt(f"{dir}/{file}", dataFull[idxSorted], delimiter=",", fmt='%.14f')

        print(f"{file} has been sorted according to '{how}'.")


def randomString(length=10):
    return ''.join(random.choices(string.ascii_lowercase + string.digits, k=length))
