import sys
import shutil
import numpy as n
import tkinter as tk
import matplotlib.pyplot as plt
from Miscellaneous import randomString
from plotFunctions import colorStationarySolution
from runSimulation import runSimulation
from tkinter import font, ttk
from subprocess import run, DEVNULL, Popen
from os.path import isfile, isdir, join
from os import listdir, remove, mkdir, rmdir


def rgb_to_hex(r, g, b):
    return ('#{:X}{:X}{:X}').format(r, g, b)


def main():
    # Simulation directory
    if len(sys.argv) > 1:
        dir = sys.argv[1]
        name = sys.argv[2]

        # Loading data
        data = n.genfromtxt(f'{dir}/solData.txt')
        with open(f'{dir}/S(x).txt', 'r') as f:
            S = f.read()
        Q = data[0]
        T = data[4:]
    else:
        dir = f'tempDir_{randomString()}'
        data = n.array([340.0, 0.0, 0.0, 0.0] + [0.0] * 12000)
        S = f"lambda x : n.exp(-n.abs(x) / 2)"
        Q = 340
        T = data[4:]
        name = 'test'

        mkdir(dir)


    # Opening simluation parameter GUI
    master = tk.Tk()
    
    # Initializing variables
    fProgress = None
    sProgress = None
    fBut = None
    sBut = None
    fTimes = None
    sTimes = None

    plotActive = False

    fSim_u = runSimulation('finitediff', data, dir, S)
    fSim_d = runSimulation('finitediff', data, dir, S)
    sSim_u = runSimulation('spectral', data, dir, S)
    sSim_d = runSimulation('spectral', data, dir, S)

    plt.ion()

    fig, ax = plt.subplots(2, 2, figsize=(8, 6))
    fig.canvas.manager.window.wm_geometry("+%d+%d" % (800, 0))
    fig.suptitle(f'Q = {Q:.2f}\n{name}')

    fig_err, ax_err_u = plt.subplots(1, 1, figsize=(16, 3))
    ax_err_d = ax_err_u.twinx()
    ax_err_u.tick_params(axis='y', colors='red')
    ax_err_d.tick_params(axis='y', colors='blue')
    fig_err.canvas.manager.window.wm_geometry("+%d+%d" % (0, 650))

    x = n.linspace(-6, 6, T.shape[0])
    # ax.plot(x, T, c='Black', ls='-', label='Stationary')
    for i in range(2):
        for j in range(2):
            # colorStationarySolution(fig, ax[i, j], data)
            ax[i, j].plot(x, data[4:], c='Black', ls='-', label='Stationary')

    ax[0, 0].set_xlabel('x (earth radiuses)')
    ax[0, 0].set_ylabel('T (℃)')
    ax[0, 0].legend()

    ax[0, 1].set_visible(False)
    ax[1, 0].set_visible(False)
    ax[1, 1].set_visible(False)
    # # ax[0, 1].plot(x, T, c='Black', ls='-', label='Stationary')
    # colorStationarySolution(fig, ax[0, 1], data)
    ax[0, 1].set_xlabel('x (earth radiuses)')
    ax[0, 1].set_ylabel('T (℃)')
    # ax[0, 1].legend()
    # # ax[1, 0].plot(x, T, c='Black', ls='-', label='Stationary')
    # colorStationarySolution(fig, ax[1, 0], data)
    ax[1, 0].set_xlabel('x (earth radiuses)')
    ax[1, 0].set_ylabel('T (℃)')
    # ax[1, 0].legend()
    # # ax[1, 1].plot(x, T, c='Black', ls='-', label='Stationary')
    # colorStationarySolution(fig, ax[1, 1], data)
    ax[1, 1].set_xlabel('x (earth radiuses)')
    ax[1, 1].set_ylabel('T (℃)')
    # ax[1, 1].legend()
    plt.tight_layout()


    def startFSimulation():
        global fProgress, fTimes
        try:
            fSim_u.times = n.array([s.split(' ') for s in tEntry.get().split('; ')], dtype=float)

            fSim_u.N = int(fEnt1.get())
            fSim_u.L = float(fEnt2.get())
            fSim_u.s = float(fEnt3.get())
            fSim_u.maxT = n.max(fSim_u.times)
            fSim_u.startSimulation(perturbation=float(sEnt4.get()))

            fSim_d.times = n.array([s.split(' ') for s in tEntry.get().split('; ')], dtype=float)

            fSim_d.N = int(fEnt1.get())
            fSim_d.L = float(fEnt2.get())
            fSim_d.s = float(fEnt3.get())
            fSim_d.maxT = n.max(fSim_d.times)
            fSim_d.startSimulation(perturbation=-float(sEnt4.get()))


            fBut.destroy()
            
            par4 = tk.Label(master, text='Progress:', fg=foregroundColor, bg=backgroundColor)
            par4.grid(row=7, column=0)

            fProgress = ttk.Progressbar(master, orient='horizontal', mode='determinate')
            fProgress.grid(row=7, column=2, pady=5)
        except NameError:
            pass
        except ValueError as e:
            print(e)

    def startSSimulation():
        global sProgress, sTimes
        try:
            sSim_u.times = n.array([s.split(' ') for s in tEntry.get().split('; ')], dtype=float)

            sSim_u.N = int(sEnt1.get())
            sSim_u.L = float(sEnt2.get())
            sSim_u.s = float(sEnt3.get())
            sSim_u.maxT = n.max(sSim_u.times)
            sSim_u.startSimulation(perturbation=float(sEnt4.get()))

            sSim_d.times = n.array([s.split(' ') for s in tEntry.get().split('; ')], dtype=float)

            sSim_d.N = int(sEnt1.get())
            sSim_d.L = float(sEnt2.get())
            sSim_d.s = float(sEnt3.get())
            sSim_d.maxT = n.max(sSim_d.times)
            sSim_d.startSimulation(perturbation=-float(sEnt4.get()))


            sBut.destroy()
            
            par4 = tk.Label(master, text='Progress:', fg=foregroundColor, bg=backgroundColor)
            par4.grid(row=7, column=0)

            sProgress = ttk.Progressbar(master, orient='horizontal', mode='determinate')
            sProgress.grid(row=7, column=1, pady=5)
        except NameError:
            pass
        except ValueError as e:
            print(e)

    def update():
        # Checking for simulation progress
        global fProgress, sProgress
        try:
            if fProgress is not None:
                val_u = fSim_u.getProgress()
                val_d = fSim_d.getProgress()
                if isinstance(val_u, float) and isinstance(val_d, float):
                    fProgress['value'] = (val_u + val_d) / 2
        except NameError:
            pass

        try:
            if sProgress is not None:
                val_u = sSim_u.getProgress()
                val_d = sSim_d.getProgress()
                if isinstance(val_u, float) and isinstance(val_d, float):
                    sProgress['value'] = (val_u + val_d) / 2
        except NameError:
            pass


        # Plotting solutions if they're finished
        if fSim_u.getResultFirst():
            ax[0, 1].set_visible(True)
            ax[1, 0].set_visible(True)
            ax[1, 1].set_visible(True)
            x, data = fSim_u.getResult()
            dx = x[1] - x[0]
            idx = int(fSim_u.times[0, 0] / n.max(fSim_u.times) * 1000)
            ax[0, 0].set_title(f't={fSim_u.times[0, 0]}')
            ax[0, 0].plot(x, data[idx, 1:], c='coral', ls='--', label='FD+')
            ax[0, 0].legend()
            idx = int(fSim_u.times[0, 1] / n.max(fSim_u.times) * 1000)
            ax[0, 1].set_title(f't={fSim_u.times[0, 1]}')
            ax[0, 1].plot(x, data[idx, 1:], c='coral', ls='--', label='FD+')
            ax[0, 1].legend()
            idx = int(fSim_u.times[1, 0] / n.max(fSim_u.times) * 1000)
            ax[1, 0].set_title(f't={fSim_u.times[1, 0]}')
            ax[1, 0].plot(x, data[idx, 1:], c='coral', ls='--', label='FD+')
            ax[1, 0].legend()
            idx = int(fSim_u.times[1, 1] / n.max(fSim_u.times) * 1000)
            ax[1, 1].set_title(f't={fSim_u.times[1, 1]}')
            ax[1, 1].plot(x, data[idx, 1:], c='coral', ls='--', label='FD+')
            ax[1, 1].legend()
            plt.tight_layout()
            fig.canvas.draw()

            err = n.zeros(1000)
            for i in range(1000):
                err[i] = n.sum((data[i, 1:] - fSim_u.T_vals) ** 2)
            t = n.linspace(0, n.max(fSim_u.times), 1000)
            ax_err_u.plot(t, err * dx, c='coral', label='FD+ error')
            ax_err_u.legend(loc=(0.9, 0.8))
            fig_err.canvas.draw()

        if fSim_d.getResultFirst():
            ax[0, 1].set_visible(True)
            ax[1, 0].set_visible(True)
            ax[1, 1].set_visible(True)
            x, data = fSim_d.getResult()
            dx = x[1] - x[0]
            idx = int(fSim_d.times[0, 0] / n.max(fSim_d.times) * 1000)
            ax[0, 0].set_title(f't={fSim_d.times[0, 0]}')
            ax[0, 0].plot(x, data[0, 1:], c='deepskyblue', ls=':', label='FD-')
            ax[0, 0].legend()
            idx = int(fSim_d.times[0, 1] / n.max(fSim_d.times) * 1000)
            ax[0, 1].set_title(f't={fSim_d.times[0, 1]}')
            ax[0, 1].plot(x, data[333, 1:], c='deepskyblue', ls=':', label='FD-')
            ax[0, 1].legend()
            idx = int(fSim_d.times[1, 0] / n.max(fSim_d.times) * 1000)
            ax[1, 0].set_title(f't={fSim_d.times[1, 0]}')
            ax[1, 0].plot(x, data[667, 1:], c='deepskyblue', ls=':', label='FD-')
            ax[1, 0].legend()
            idx = int(fSim_d.times[1, 1] / n.max(fSim_d.times) * 1000)
            ax[1, 1].set_title(f't={fSim_d.times[1, 1]}')
            ax[1, 1].plot(x, data[1000, 1:], c='deepskyblue', ls=':', label='FD-')
            ax[1, 1].legend()
            plt.tight_layout()
            fig.canvas.draw()

            err = n.zeros(1000)
            for i in range(1000):
                err[i] = n.sum((data[i, 1:] - fSim_d.T_vals) ** 2)
            t = n.linspace(0, n.max(fSim_u.times), 1000)
            ax_err_d.plot(t, err * dx, c='deepskyblue', label='FD- error')
            ax_err_d.legend(loc=(0.9, 0.7))
            fig_err.canvas.draw()
        
        if sSim_u.getResultFirst():
            ax[0, 1].set_visible(True)
            ax[1, 0].set_visible(True)
            ax[1, 1].set_visible(True)
            x, data = sSim_u.getResult()
            idx = int(fSim_d.times[0, 0] / n.max(fSim_d.times) * 1000)
            ax[0, 0].set_title(f't={sSim_u.times[0, 0]}')
            ax[0, 0].plot(x, data[0, 1:], c='sandybrown', ls='--', label='S+')
            ax[0, 0].legend()
            idx = int(fSim_d.times[0, 1] / n.max(fSim_d.times) * 1000)
            ax[0, 1].set_title(f't={sSim_u.times[0, 1]}')
            ax[0, 1].plot(x, data[333, 1:], c='sandybrown', ls='--', label='S+')
            ax[0, 1].legend()
            idx = int(fSim_d.times[1, 0] / n.max(fSim_d.times) * 1000)
            ax[1, 0].set_title(f't={sSim_u.times[1, 0]}')
            ax[1, 0].plot(x, data[667, 1:], c='sandybrown', ls='--', label='S+')
            ax[1, 0].legend()
            idx = int(fSim_d.times[1, 1] / n.max(fSim_d.times) * 1000)
            ax[1, 1].set_title(f't={sSim_u.times[1, 1]}')
            ax[1, 1].plot(x, data[1000, 1:], c='sandybrown', ls='--', label='S+')
            ax[1, 1].legend()
            plt.tight_layout()
            fig.canvas.draw()
        
        if sSim_d.getResultFirst():
            ax[0, 1].set_visible(True)
            ax[1, 0].set_visible(True)
            ax[1, 1].set_visible(True)
            x, data = sSim_d.getResult()
            idx = int(fSim_d.times[0, 0] / n.max(fSim_d.times) * 1000)
            ax[0, 0].set_title(f't={sSim_d.times[0, 0]}')
            ax[0, 0].plot(x, data[0, 1:], c='darkorange', ls='--', label='S-')
            ax[0, 0].legend()
            idx = int(fSim_d.times[0, 1] / n.max(fSim_d.times) * 1000)
            ax[0, 1].set_title(f't={sSim_d.times[0, 1]}')
            ax[0, 1].plot(x, data[333, 1:], c='darkorange', ls='--', label='S-')
            ax[0, 1].legend()
            idx = int(fSim_d.times[1, 0] / n.max(fSim_d.times) * 1000)
            ax[1, 0].set_title(f't={sSim_d.times[1, 0]}')
            ax[1, 0].plot(x, data[667, 1:], c='darkorange', ls='--', label='S-')
            ax[1, 0].legend()
            idx = int(fSim_d.times[1, 1] / n.max(fSim_d.times) * 1000)
            ax[1, 1].set_title(f't={sSim_d.times[1, 1]}')
            ax[1, 1].plot(x, data[1000, 1:], c='darkorange', ls='--', label='S-')
            ax[1, 1].legend()
            plt.tight_layout()
            fig.canvas.draw()

        if fSim_u.result is not None or sSim_u.result is not None or fSim_d.result is not None or sSim_d.result is not None:
            fig.canvas.flush_events()

        master.after(100, update)

    # Defining some colors
    backgroundColor = rgb_to_hex(20, 20, 20)
    foregroundColor = rgb_to_hex(220, 220, 220)

    # Setting up title
    master.configure(bg=backgroundColor)
    tk.Label(master, text=f'Simulation parameters for Q={Q:.1f} in solution "{name}"', 
        fg=foregroundColor, bg=backgroundColor, font=('Cambria Math', 12)).grid(row=0, columnspan=3)


    # First row
    sLab0 = tk.Label(master, text='Spectral:', fg=foregroundColor, bg=backgroundColor)
    sLab0.grid(row=1, column=1)

    fLab0 = tk.Label(master, text='Finite diff:', fg=foregroundColor, bg=backgroundColor)
    fLab0.grid(row=1, column=2)


    # Choosing parameter N
    par1 = tk.Label(master, text='N', fg=foregroundColor, bg=backgroundColor)
    par1.grid(row=2, column=0)

    sEnt1 = tk.Entry(master, fg=foregroundColor, bg=backgroundColor)
    sEnt1.insert(0, 1024)
    sEnt1.grid(row=2, column=1)

    fEnt1 = tk.Entry(master, fg=foregroundColor, bg=backgroundColor)
    fEnt1.insert(0, 1001)
    fEnt1.grid(row=2, column=2)


    # Choosing parameter L
    par2 = tk.Label(master, text='L', fg=foregroundColor, bg=backgroundColor)
    par2.grid(row=3, column=0)

    sEnt2 = tk.Entry(master, fg=foregroundColor, bg=backgroundColor)
    sEnt2.insert(0, 6)
    sEnt2.grid(row=3, column=1)

    fEnt2 = tk.Entry(master, fg=foregroundColor, bg=backgroundColor)
    fEnt2.insert(0, 6)
    fEnt2.grid(row=3, column=2)


    # Choosing parameter s
    par3 = tk.Label(master, text='s', fg=foregroundColor, bg=backgroundColor)
    par3.grid(row=4, column=0)

    sEnt3 = tk.Entry(master, fg=foregroundColor, bg=backgroundColor)
    sEnt3.insert(0, 0.001)
    sEnt3.grid(row=4, column=1)

    fEnt3 = tk.Entry(master, fg=foregroundColor, bg=backgroundColor)
    fEnt3.insert(0, 0.001)
    fEnt3.grid(row=4, column=2)


    # Choosing what perturbation
    par4 = tk.Label(master, text='pert.', fg=foregroundColor, bg=backgroundColor)
    par4.grid(row=5, column=0)

    sEnt4 = tk.Entry(master, fg=foregroundColor, bg=backgroundColor)
    sEnt4.insert(0, 0.1)
    sEnt4.grid(row=5, column=1, columnspan=3)


    # Choosing what times to plot
    par4 = tk.Label(master, text='t', fg=foregroundColor, bg=backgroundColor)
    par4.grid(row=6, column=0)

    tEntry = tk.Entry(master, fg=foregroundColor, bg=backgroundColor)
    tEntry.insert(0, "0.0 16.7; 33.3 50.0")
    tEntry.grid(row=6, column=1, columnspan=3)

    
    # Start simulation buttons
    sBut = tk.Button(master, text='Start', command=startSSimulation, fg=foregroundColor, bg=backgroundColor)
    sBut.grid(row=7, column=1)

    fBut = tk.Button(master, text='Start', command=startFSimulation, fg=foregroundColor, bg=backgroundColor)
    fBut.grid(row=7, column=2)

    # Cleanup code
    def cleanup():
        # Killing simulations if they're running and deletes directory
        if fSim_u.sim:
            run("TASKKILL /F /PID {pid} /T".format(pid=fSim_u.sim.pid), stdout=DEVNULL, stderr=DEVNULL)
        if fSim_d.sim:
            run("TASKKILL /F /PID {pid} /T".format(pid=fSim_d.sim.pid), stdout=DEVNULL, stderr=DEVNULL)
        if sSim_u.sim:
            run("TASKKILL /F /PID {pid} /T".format(pid=sSim_u.sim.pid), stdout=DEVNULL, stderr=DEVNULL)
        if sSim_d.sim:
            run("TASKKILL /F /PID {pid} /T".format(pid=sSim_d.sim.pid), stdout=DEVNULL, stderr=DEVNULL)
        try:
            shutil.rmtree(dir)
        except FileNotFoundError:
            pass
        master.quit()
        plt.close()

    
    def reset():
        # Killing simulations if they're running and deletes directory
        if fSim_u.sim:
            run("TASKKILL /F /PID {pid} /T".format(pid=fSim_u.sim.pid), stdout=DEVNULL, stderr=DEVNULL)
        if fSim_d.sim:
            run("TASKKILL /F /PID {pid} /T".format(pid=fSim_d.sim.pid), stdout=DEVNULL, stderr=DEVNULL)
        if sSim_u.sim:
            run("TASKKILL /F /PID {pid} /T".format(pid=sSim_u.sim.pid), stdout=DEVNULL, stderr=DEVNULL)
        if sSim_d.sim:
            run("TASKKILL /F /PID {pid} /T".format(pid=sSim_d.sim.pid), stdout=DEVNULL, stderr=DEVNULL)
        master.quit()
        plt.close()
        Popen(['python', 'simulationGUI.py', dir, name])

    resetBut = tk.Button(master, text='Reset', command=reset, fg=foregroundColor, bg=backgroundColor)
    resetBut.grid(row=8, column=1, columnspan=3)

    # Tk main loop
    master.after(0, update)
    master.protocol("WM_DELETE_WINDOW", cleanup)
    master.geometry('%dx%d+%d+%d' % (650, 300, 10, 10))
    master.mainloop()


if __name__ == '__main__':
    main()
