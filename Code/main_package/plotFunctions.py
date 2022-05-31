import numpy as n
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.widgets import Slider


def consecutive(data, stepsize=1):
    return n.split(data, n.where(n.diff(data) != stepsize)[0] + 1)


def highlightSolution(fig, ax, solutionData):

    # Loading data
    sols = solutionData.getSols()
    activeSol = None
    
    # Initializing point that shows active solution
    point = ax.scatter([n.nan], [n.nan], s=5, facecolors='none', edgecolors='Black')

    # Defining funciton to run on mouse-move event
    def onmove(event):
        global activeSol
        """
        Event handeler for clicking on a datapoint
        """
        # Finding family closest to point
        try:
            # Important data about size of solution domain
            x_l, x_h = ax.get_xlim()
            y_l, y_h = ax.get_ylim()

            xL = x_h - x_l  # Sixe of domain on x-axis
            yL = y_h - y_l  # Sixe of domain on y-axis

            # Transforming to normalized coordiantes
            mousePos = n.array([event.xdata, event.ydata])

            # Computing distances to mouse
            dists = ((sols[:, 0] - mousePos[0]) / xL) ** 2 + \
                    ((sols[:, 1] - mousePos[1]) / yL) ** 2

            # Finds solution closest to mouse
            closest = n.argmin(dists)

            # Checks if it is close enough
            if n.min(dists) < 0.025 ** 2:
                point.set_offsets(sols[closest])
                activeSol = closest
            else:
                point.set_offsets([n.nan, n.nan])
                activeSol = None

            # Draw the point
            fig.canvas.draw_idle()
        except TypeError:
            # Removes active solution and dot
            point.set_offsets([n.nan, n.nan])
            activeSol = None

    # Defining funciton to run on mouse-click event
    def onclick(event):
        global activeSol
        if event.button == 3 and activeSol is not None:
            # Runs simulation of highlighted solution
            try:
                solutionData(activeSol)
            except ValueError as e:
                print(f"Invalid data in {solutionData.getSolutionName(activeSol)}, run processFiles()")

    fig.canvas.mpl_connect('button_press_event', onclick)
    fig.canvas.mpl_connect('motion_notify_event', onmove)


def colorStationarySolution(fig, ax, data, scale=1):
    C_0 = data[2]
    C_1 = data[3]
    T = data[4:]

    x = n.linspace(-6, 6, T.shape[0])
    idx = n.arange(T.shape[0])
    l_0, l_1 = n.argmin(n.abs(x - C_0)), n.argmin(n.abs(x - C_1))

    ice = n.argwhere(T[idx] < -10 / scale)[:, 0]
    ice = ice[(ice <= l_0) | (ice >= l_1)]
    water = n.argwhere(T[idx] >= -10 / scale)[:, 0]
    water = water[(water <= l_0) | (water >= l_1)]

    snow = n.argwhere(T[idx] < 0)[:, 0]
    snow = snow[(snow >= l_0) & (snow <= l_1)]
    land = n.argwhere(T[idx] >= 0)[:, 0]
    land = land[(land >= l_0) & (land <= l_1)]

    if ice.shape[0]:
        ax.plot(n.nan, n.nan, c="lightblue", path_effects=[pe.Stroke(linewidth=3, foreground='black'), pe.Normal()], label='Ice')
        segments = consecutive(ice)
        for seg in segments:
            ax.plot(x[idx[seg]], T[idx[seg]], c="lightblue", path_effects=[pe.Stroke(linewidth=3, foreground='black'), pe.Normal()])

    if water.shape[0]:
        ax.plot(n.nan, n.nan, c="blue", path_effects=[pe.Stroke(linewidth=3, foreground='black'), pe.Normal()], label='Water')
        segments = consecutive(water)
        for seg in segments:
            ax.plot(x[idx[seg]], T[idx[seg]], c="blue", path_effects=[pe.Stroke(linewidth=3, foreground='black'), pe.Normal()])

    if snow.shape[0]:
        ax.plot(n.nan, n.nan, c="snow", path_effects=[pe.Stroke(linewidth=3, foreground='black'), pe.Normal()], label='Snow')
        segments = consecutive(snow)
        for seg in segments:
            ax.plot(x[idx[seg]], T[idx[seg]], c="snow", path_effects=[pe.Stroke(linewidth=3, foreground='black'), pe.Normal()])

    if land.shape[0]:
        ax.plot(n.nan, n.nan, c="brown", path_effects=[pe.Stroke(linewidth=3, foreground='black'), pe.Normal()], label='Land')
        segments = consecutive(land)
        for seg in segments:
            ax.plot(x[idx[seg]], T[idx[seg]], c="brown", path_effects=[pe.Stroke(linewidth=3, foreground='black'), pe.Normal()])
