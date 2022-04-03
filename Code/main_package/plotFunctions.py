import numpy as n
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider


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
            solutionData(activeSol)

    fig.canvas.mpl_connect('button_press_event', onclick)
    fig.canvas.mpl_connect('motion_notify_event', onmove)
