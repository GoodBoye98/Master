import subprocess
import numpy as n
import pandas as pd
from os import listdir, mkdir, rmdir
from os.path import isfile, join
from main_package.Miscellaneous import randomString
from main_package.solutionFamily import solutionFamily


class solutionData:
    def __init__(self, mainDir, subDir, usecols=[0, 1]):
        """
        Load all solutions from a directory.
        """
        self.name = f"{mainDir} + {subDir}"
        if mainDir == "Exponential":
            if subDir == "Water":
                self.C_0 = 0
                self.C_1 = 0
            elif "Continent Offset" in subDir:
                value = float(subDir.replace("Continent Offset ", ""))
                self.C_0 = -1 + value
                self.C_1 = 1 + value
            self.S = f"lambda x : n.exp(-n.abs(x) / 2)"
            self.S_latex = r"$S(x)=e^{-\frac{|x|}{2}}$"
        elif mainDir == "Gaussian":
            val = 10
            if subDir == "Water":
                self.C_0 = 0
                self.C_1 = 0
            elif "Std" in subDir:
                self.C_0 = -1
                self.C_1 = 1
                val = float(subDir.replace("Std ", ""))
            self.S = f"lambda x: 4 / (n.sqrt(n.pi * {val})) * n.exp(- x * x / {val})"
            self.S_latex = r"$S(x)=\frac{4}{\sqrt{V \pi}}e^{-\frac{x^2}{V}}$".replace('V', f"{val:.0f}")
        else:
            raise OSError(f"Directory '{mainDir}' does not exist")

        self.families = []
        self.params = n.array([])

        self.dir = f"{mainDir}/{subDir}"
        try:
            self.params = n.genfromtxt(self.dir+'/Measurements.txt', delimiter=', ', dtype=str)
        except UnboundLocalError:
            print('Warming: No Measurements.txt file found. Defaults to [Q, T_0]')
            self.params = n.array(['Q', 'T_0'], dtype=str)
        files = [f for f in listdir(self.dir) if isfile(join(self.dir, f))]
        for file in files:
            if file == "Measurements.txt":
                continue
            dataDir = f"{mainDir}/{subDir}/{file}"
            data = pd.read_csv(dataDir, header=None, usecols=usecols, delimiter=',').values
            sol = solutionFamily(data, name=file.replace('.csv', ''), params=self.params)
            self.families.append(sol)


    def __call__(self, index, run=True):
        for family in self.families:
            if index < family.Q.shape[0]:
                # Loads solution data
                loaddir = f"{self.dir}/{family.name}.csv"
                solution = pd.read_csv(loaddir, header=None, skiprows=index, nrows=1).values[0]
                
                if run is False:
                    return n.concatenate(([solution[0]], solution[self.params.shape[0]:]))

                # Saves solution data in new directory
                solData = f'{solution[0]} {solution[1]} {self.C_0} {self.C_1} ' + ' '.join(f'{10 * v:.14f}' for v in solution[self.params.shape[0]:])
                
                # Creates simulation directory if data was loaded succesfully
                newDir = f'tempDir_{randomString(15)}'
                mkdir(f'main_package/{newDir}')

                with open(f'main_package/{newDir}/solData.txt', 'w') as file:
                    file.write(solData)
                with open(f'main_package/{newDir}/S(x).txt', 'w') as file:
                    file.write(self.S)
                
                # Runs simulation
                subprocess.Popen(['python', 'simulationGUI.py', newDir, family.name], cwd='main_package', shell=True)
                return None
            else:
                index -= family.Q.shape[0]


    def plot(self, fig, ax, mode='plot', **kwargs):
        fig.suptitle(f"Bifurcation diagram\n{self.S_latex}")
        for family in self.families:
            family.plot(fig, ax, mode, **kwargs)
        ax.legend()


    def getSol(self, index, all=False, specialFamily=None):
        for family in self.families:
            if specialFamily is not None:
                if family.name == specialFamily:
                    loaddir = f"{self.dir}/{family.name}.csv"
                    solution = pd.read_csv(loaddir, header=None, skiprows=index, nrows=1).values[0]
                    return solution[self.params.shape[0]:] if not all else solution
                continue

            if index < family.Q.shape[0]:
                # Loads solution data
                loaddir = f"{self.dir}/{family.name}.csv"
                solution = pd.read_csv(loaddir, header=None, skiprows=index, nrows=1).values[0]

                # Returns solution
                return solution[self.params.shape[0]:] if not all else solution
            else:
                index -= family.Q.shape[0]
        return None


    def getSols(self, measure='0'):
        data = []
        for family in self.families:
            data.append(family.getSols(measure))
        return n.vstack(data)

    def getSolutionFamily(self, name):
        for family in self.families:
            if family.name == name:
                return family.getSols()

    def getSolutionName(self, index):
        for family in self.families:
            if index < family.Q.shape[0]:
                return family.name
            else:
                index -= family.Q.shape[0]
