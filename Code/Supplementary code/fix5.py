import numpy as n
import pandas as pd
import matplotlib.pyplot as plt
from plot import solFamily, solution

if __name__ == '__main__':
    offset = '0.05'
    cnt = [-1 + float(offset), 1 + float(offset)]

    thing = pd.read_csv('0.05/5 - IceWaterLandSnowWaterIce.csv', header=None).values
    arr = thing[:, 3]
    idxSort = n.argsort(arr)
    thing = thing[idxSort]
    upper = thing[n.where(thing[:, 3] > 0.92)]
    lower = thing[n.where(thing[:, 3] < 0.92)]
    upper = upper[n.argsort(upper[:, 0])]
    lower = lower[n.argsort(lower[:, 0])]

    i = 0
    while i < lower.shape[0]:
        p1 = n.array([lower[i, 0], lower[i, 3]])

        delete = []

        for j in range(i, lower.shape[0] - 1):
            p2 = n.array([lower[j, 0], lower[j, 3]])

            dist = n.sum((p2 - p1) ** 2)
            if 0 < dist < 0.3 ** 2:
                delete.append(j)
        delete = n.array(delete, dtype=int)
        lower = n.delete(lower, delete, axis=0)

        i += 1
    i = 0
    while i < upper.shape[0]:
        p1 = n.array([upper[i, 0], upper[i, 3]])

        delete = []

        for j in range(i, upper.shape[0] - 1):
            p2 = n.array([upper[j, 0], upper[j, 3]])

            dist = n.sum((p2 - p1) ** 2)
            if 0 < dist < 0.3 ** 2:
                delete.append(j)
        delete = n.array(delete, dtype=int)
        upper = n.delete(upper, delete, axis=0)

        i += 1


    matrix = n.concatenate((upper, lower))
    plt.scatter(matrix[:, 0], matrix[:, 3])
    plt.show()


    # n.savetxt('0.05/test.csv', matrix, delimiter=', ', fmt='%.10f')
