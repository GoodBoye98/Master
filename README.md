# Master
Code written for master thesis "An Energy Balance Mode on an Infinite Line" (only for Windows)

To run the different programs you need to download data first.
For the Laplacian case, download the folder "Exponential" from ...
For the Gaussian case, download the folder "Gaussian" from shorturl.at/oxFN6
These folders are to be placed in the "Code" folder. In the same directory as the "main.py" file

In the Code folder there are a bunch of example scripts, they were used to plot various figures in the thesis.
To plot bifurcation and run the finite difference or spectral numerical schemes, run "main.py" from its own directory
In "main.py", change sType and dType to plot the different bifurcation diagrams, the possible combinations are

sType = "Exponential"
  - dType = "Water"
  - dType = "Continential Offset 0.000"
  - dType = "Continential Offset 0.001"
  - dType = "Continential Offset 0.003"
  - dType = "Continential Offset 0.005"
  - dType = "Continential Offset 0.010"
  - dType = "Continential Offset 0.015"
  - dType = "Continential Offset 0.020"
  - dType = "Continential Offset 0.030"
  - dType = "Continential Offset 0.040"
  - dType = "Continential Offset 0.050"

sType = "Gaussian"
  - dType = "Water"
  - dType = "Std 5"
  - dType = "Std 6"

dType = "Water" corresponds to the line without a continent, "Std 5" is the narrower gaussian and "Std 6" is the wider gaussian.

To plot a solution from the bifurcation diagram and run the numerical schemes for a perturbation, simply right click on a solution, this should open up a GUI for running the schemes. To do that, it is requiered to compile the c++ code in "Code/main_package/C++/numerics.cpp", which depends upon the library Eigen 3.4.0 (https://eigen.tuxfamily.org/index.php?title=Main_Page). The GUI gives the oppurtunity to change N, L, s, the perturbation and what t values to plot from the simulations.
