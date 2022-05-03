import numpy as n
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm
from main_package import solutionData

def main():
    eData = solutionData('Gaussian', 'Test')

    fig, ax = plt.subplots(1, 1, figsize=(8, 6))

    eData.plot(fig, ax, 'plot')
    plt.tight_layout()
    plt.show()

    # higher = '3A - IceWater_Snow_WaterIce_u_'
    # lower = '3A - IceWater_Snow_WaterIce_'

    # h_tot = n.zeros((0, 12003))
    # l_tot = n.zeros((0, 12003))

    # for i in tqdm(range(1, 20)):
    #     h = pd.read_csv(f'Gaussian/Test/{higher}{i}.csv', header=None, delimiter=',').values
    #     Q = h[:, 0]
    #     h = h[n.argsort(Q)]
    #     if not i % 2:
    #         h = h[::-1]
    #     h_tot = n.vstack((h_tot, h))
        

    #     l = pd.read_csv(f'Gaussian/Test/{lower}{i}.csv', header=None, delimiter=',').values
    #     Q = l[:, 0]
    #     l = l[n.argsort(Q)]
    #     if i % 2:
    #         l = l[::-1]
    #     l_tot = n.vstack((l_tot, l))
    
    # n.savetxt(f'Gaussian/Test/{higher}{i}.csv', h_tot, delimiter=',', fmt='%.14f')
    # n.savetxt(f'Gaussian/Test/{lower}{i}.csv', l_tot, delimiter=',', fmt='%.14f')


    

        

if __name__ == '__main__':
    main()