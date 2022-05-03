import time
import shutil
import numpy as n
import matplotlib.pyplot as plt
from tqdm import tqdm
from matplotlib import cm
from os.path import isdir
from os import mkdir, rmdir
from multiprocessing import Pool
from main_package import solutionData, runSimulation


def secToString(t):
    days = int(t // 86400)
    t -= days * 86400

    hours = int(t // 3600)
    t -= hours * 3600

    minutes = int(t // 60)
    t -= minutes * 60

    seconds = int(t // 1)
    t -= seconds * 1

    millis = int(t * 1000 + 0.5)

    return f'{days:>2}d{hours:>2}h{minutes:>2}m{seconds:>2}s{millis:>3}ms'.replace(' ', '0')


def testStability(args):
    idx, sType, dType, C_0, C_1, pert = args

    sData = solutionData(sType, dType)
    data = sData(idx, run=False)
    data = n.concatenate(([data[0]], [0], [C_0], [C_1], 10 * data[1:]))

    k = 'u' if pert > 0 else 'd'
    dirName = f'simulations/tempDir_stability_{k}_{str(idx):>4}'.replace(' ', '0')

    mkdir(dirName)

    simRunner = runSimulation('finitediff', data, dirName, sData.S)

    simRunner.N = 4001
    simRunner.L = 6
    simRunner.s = 0.001
    simRunner.maxT = 50
    simRunner.savePoints = 2
    simRunner.startSimulation(perturbation=pert)

    while True:
        if simRunner.getResultFirst(rows=[-1]):
            break
        time.sleep(0.1)
    x, T_sim = simRunner.getResult()
    dx = x[1] - x[0]
    error = n.sum((simRunner.T_vals - T_sim[0, 1:]) ** 2) * dx

    return [idx, pert, error]


def run():
    if isdir('simulations'):
        shutil.rmtree('simulations')
    mkdir('simulations')

    sType = 'Gaussian'
    dType = 'Water'
    continentOffset = 0.0
    C_0 = 0
    C_1 = 0

    sData = solutionData(sType, dType)

    sols = sData.getSols()
    idx = n.arange(0, sols.shape[0], dtype=int)
    print(f'Max index : {idx[-1]}')

    args = []
    for i in idx:
        args.append([i, sType, dType, C_0, C_1, 0.1])
        args.append([i, sType, dType, C_0, C_1, -0.1])


    t_0 = time.time()

    res = []

    with Pool(16) as p:
        for x in tqdm(p.imap_unordered(testStability, args), total=idx.shape[0] * 2):
            res.append(x)
        
    t_1 = time.time()

    shutil.rmtree('simulations')
    print(secToString(t_1 - t_0))

    res = n.array(res)
    n.save(f'Stability/{sType}-{dType}-stability-N4001.npy', res)


def plot(arg):

    sType = 'Exponential'
    dType = 'Water'

    data = n.load(f'Stability/{sType}-{dType}-stability-N2001.npy')

    sData = solutionData(sType, dType)
    sols = sData.getSols()

    idx = data[:, 0]
    pert = data[:, 1]
    error = data[:, 2]

    idx_u = idx[pert > 0].astype(int)
    idx_d = idx[pert < 0].astype(int)

    error_u = error[pert > 0]
    error_d = error[pert < 0]

    fig, ax = plt.subplots(1, 1, figsize=(6, 4))

    if arg == '+':

        scatter = ax.scatter(sols[idx_u, 0], sols[idx_u, 1], 
            s=1, c=error_u, cmap='rainbow', vmin=0, vmax=0.2)
        ax.set_title('Upward perturbation')
        ax.set_xlabel('Q')
        ax.set_ylabel('T(0)')
    
    if arg == '-':

        scatter = ax.scatter(sols[idx_d, 0], 
            sols[idx_d, 1], s=1, c=error_d, cmap='rainbow', vmin=0, vmax=0.2)
        ax.set_title('Downward perturbation')
        ax.set_xlabel('Q')
        ax.set_ylabel('T(0)')

    bar = fig.colorbar(scatter)
    bar.ax.set_title('Error')
    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    # run()
    plot('+')
    plot('-')