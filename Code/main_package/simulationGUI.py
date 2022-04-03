import sys
import shutil
import numpy as n
import tkinter as tk
import matplotlib.pyplot as plt
from Miscellaneous import randomString
from runSimulation import runSimulation
from tkinter import font, ttk
from subprocess import run, DEVNULL
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

    fSim = runSimulation('finitediff', data, dir, S)
    sSim = runSimulation('spectral', data, dir, S)

    plt.ion()
    fig, ax = plt.subplots(2, 2, figsize=(8, 6))
    fig.suptitle(f'Q = {Q}\n{name}')

    x = n.linspace(-6, 6, T.shape[0])
    ax[0, 0].plot(x, T, c='Black', ls='-', label='Stationary')
    ax[0, 0].set_xlabel('x (earth radiuses)')
    ax[0, 0].set_ylabel('T (℃)')
    ax[0, 0].legend()
    ax[0, 1].plot(x, T, c='Black', ls='-', label='Stationary')
    ax[0, 1].set_xlabel('x (earth radiuses)')
    ax[0, 1].set_ylabel('T (℃)')
    ax[0, 1].legend()
    ax[1, 0].plot(x, T, c='Black', ls='-', label='Stationary')
    ax[1, 0].set_xlabel('x (earth radiuses)')
    ax[1, 0].set_ylabel('T (℃)')
    ax[1, 0].legend()
    ax[1, 1].plot(x, T, c='Black', ls='-', label='Stationary')
    ax[1, 1].set_xlabel('x (earth radiuses)')
    ax[1, 1].set_ylabel('T (℃)')
    ax[1, 1].legend()
    plt.tight_layout()


    def startFSimulation():
        global fProgress, fTimes
        try:
            fSim.N = int(fEnt1.get())
            fSim.L = float(fEnt2.get())
            fSim.s = float(fEnt3.get())
            fSim.startSimulation()

            fSim.times = n.array([s.split(' ') for s in tEntry.get().split('; ')], dtype=float)
            print(fSim.times)

            fBut.destroy()
            
            par4 = tk.Label(master, text='Progress:', fg=foregroundColor, bg=backgroundColor)
            par4.grid(row=6, column=0)

            fProgress = ttk.Progressbar(master, orient='horizontal', mode='determinate')
            fProgress.grid(row=6, column=2, pady=5)
        except NameError:
            pass
        except ValueError as e:
            print(e)

    def startSSimulation():
        global sProgress, sTimes
        try:
            sSim.N = int(sEnt1.get())
            sSim.L = float(sEnt2.get())
            sSim.s = float(sEnt3.get())
            sSim.startSimulation()

            sSim.times = n.array([s.split(' ') for s in tEntry.get().split('; ')], dtype=float)
            print(sSim.times)

            sBut.destroy()
            
            par4 = tk.Label(master, text='Progress:', fg=foregroundColor, bg=backgroundColor)
            par4.grid(row=6, column=0)

            sProgress = ttk.Progressbar(master, orient='horizontal', mode='determinate')
            sProgress.grid(row=6, column=1, pady=5)
        except NameError:
            pass
        except ValueError as e:
            print(e)

    def update():
        # Checking for simulation progress
        global fProgress, sProgress
        try:
            if fProgress is not None:
                val = fSim.getProgress()
                if val is not None and val is not False:
                    fProgress['value'] = val
        except NameError:
            pass

        try:
            if sProgress is not None:
                val = sSim.getProgress()
                if val is not None and val is not False:
                    sProgress['value'] = val
        except NameError:
            pass


        # Plotting solutions if they're finished
        if fSim.getResultFirst():
            x, data = fSim.getResult()
            ax[0, 0].set_title(f't={fSim.times[0, 0]}')
            ax[0, 0].plot(x, data[0, 1:], c='blue', ls='-.', label='Finite difference')
            ax[0, 0].legend()
            ax[0, 1].set_title(f't={fSim.times[0, 1]}')
            ax[0, 1].plot(x, data[333, 1:], c='blue', ls='-.', label='Finite difference')
            ax[0, 1].legend()
            ax[1, 0].set_title(f't={fSim.times[1, 0]}')
            ax[1, 0].plot(x, data[667, 1:], c='blue', ls='-.', label='Finite difference')
            ax[1, 0].legend()
            ax[1, 1].set_title(f't={fSim.times[1, 1]}')
            ax[1, 1].plot(x, data[1000, 1:], c='blue', ls='-.', label='Finite difference')
            ax[1, 1].legend()
            plt.tight_layout()
            fig.canvas.draw()
        
        if sSim.getResultFirst():
            x, data = sSim.getResult()
            ax[0, 0].set_title(f't={fSim.times[0, 0]}')
            ax[0, 0].plot(x, data[0, 1:], c='orange', ls='--', label='Spectral')
            ax[0, 0].legend()
            ax[0, 1].set_title(f't={fSim.times[0, 1]}')
            ax[0, 1].plot(x, data[333, 1:], c='orange', ls='--', label='Spectral')
            ax[0, 1].legend()
            ax[1, 0].set_title(f't={fSim.times[1, 0]}')
            ax[1, 0].plot(x, data[667, 1:], c='orange', ls='--', label='Spectral')
            ax[1, 0].legend()
            ax[1, 1].set_title(f't={fSim.times[1, 1]}')
            ax[1, 1].plot(x, data[1000, 1:], c='orange', ls='--', label='Spectral')
            ax[1, 1].legend()
            plt.tight_layout()
            fig.canvas.draw()

        if fSim.result is not None or sSim.result is not None:
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
    sEnt1.insert(0, 4096)
    sEnt1.grid(row=2, column=1)

    fEnt1 = tk.Entry(master, fg=foregroundColor, bg=backgroundColor)
    fEnt1.insert(0, 4001)
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


    # Choosing what times to plot
    par4 = tk.Label(master, text='t', fg=foregroundColor, bg=backgroundColor)
    par4.grid(row=5, column=0)

    tEntry = tk.Entry(master, fg=foregroundColor, bg=backgroundColor)
    tEntry.insert(0, "0.0 1.0; 2.0 3.0")
    tEntry.grid(row=5, column=1, columnspan=3)

    
    # Start simulation buttons
    sBut = tk.Button(master, text='Start', command=startSSimulation, fg=foregroundColor, bg=backgroundColor)
    sBut.grid(row=6, column=1)

    fBut = tk.Button(master, text='Start', command=startFSimulation, fg=foregroundColor, bg=backgroundColor)
    fBut.grid(row=6, column=2)

    # Cleanup code
    def cleanup():
        # Killing simulations if they're running and deletes directory
        if fSim.sim:
            run("TASKKILL /F /PID {pid} /T".format(pid=fSim.sim.pid), stdout=DEVNULL, stderr=DEVNULL)
        if sSim.sim:
            run("TASKKILL /F /PID {pid} /T".format(pid=sSim.sim.pid), stdout=DEVNULL, stderr=DEVNULL)
        try:
            shutil.rmtree(dir)
        except FileNotFoundError:
            pass
        master.quit()
        plt.close()


    # Tk main loop
    master.after(0, update)
    master.protocol("WM_DELETE_WINDOW", cleanup)
    master.mainloop()


if __name__ == '__main__':
    main()