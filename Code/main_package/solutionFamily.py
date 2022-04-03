import numpy as n


class solutionFamily:
    def __init__(self, data, name, params):
        """
        Helper class for a single solution in a family.
        """
        self.Q = self.T_0 = self.T_c = self.T_m = None
        self.name = name
        for i, par in enumerate(params):
            if par == 'Q':
                self.Q = data[:, i]
            if par == 'T_0':
                self.T_0 = data[:, i]
            if par == 'T_center':
                self.T_c = data[:, i]
            if par == 'T_max':
                self.T_m = data[:, i]

    def getSols(self, measure='0'):
        """
        Get coordinates of all solutions in this family based upon chosen merit.
        T_0: temperature at x=0.
        T_center: temerature at the center of the continent.
            (Same as T_0 in the case of no or centered continent)
        T_max: max temperature on the line or continent if there is one
        """
        if measure.lower() in ['t_0', '0', 'default']:
            return n.column_stack([self.Q, self.T_0])
        if measure.lower() in ['t_m', 'max', 'maximum', 't_max']:
            return n.column_stack([self.Q, self.T_m])
        if measure.lower() in ['t_c', 'center', 'centre', 't_center', 't_centre']:
            return n.column_stack([self.Q, self.T_c])

    def plot(self, fig, ax, mode='plot'):
        if mode == 'plot':
            ax.plot(self.Q, self.T_0, label=self.name)
        elif mode == 'scatter':
            ax.scatter(self.Q, self.T_0, label=self.name, s=1)