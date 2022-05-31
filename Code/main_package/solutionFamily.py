import numpy as n


COLORS = {
        "1 - Ice": '#1f77b4',
        "2 - IceWaterIce": '#ff7f0e',

        "1 - IceSnowIce": '#1f77b4',
        "2 - IceWaterSnowWaterIce": '#ff7f0e',
        "3 - IceWaterSnowLandSnowWaterIce": '#2ca02c',
        "3 - IceWaterLandSnowLandWaterIce": '#2ca02c',
        "4 - IceWaterLandWaterIce": '#d62728',

#       "1 - IceSnowIce": '#1f77b4',
        "2 - IceWaterSnowIce": '#ff7f0e',
        "3 - IceWaterSnowWaterIce": '#8c564b',
        "4 - IceWaterSnowLandSnowWaterIce": '#2ca02c',
        "3A - IceWater_Land_WaterIce": '#2ca02c',
        "3B - IceWater_Snow_WaterIce": '#9467bd',
        "4B - IceWaterSnowLandSnowIce": '#e377c2',
        "5A - IceWaterLandSnowWaterIce": '#7f7f7f',
        "5B - IceWaterSnowLandWaterIce": '#17becf',
        "5C - IceWaterLandSnowIce": '#bcbd22',
        "5D - IceSnowLandWaterIce": '#e377c2',
        "6 - IceWaterLandWaterIce": '#d62728'
}



class solutionFamily:
    def __init__(self, data, name, params):
        """
        Helper class for a single solution in a family.
        """
        self.Q = self.T_0 = self.T_c = self.T_m = None
        # self.name = name if name != "6 - IceWaterLandWaterIce" else "4 - IceWaterLandWaterIce"  # Cause I'm stupid and fked the indexing of the stability simulation
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

    def plot(self, fig, ax, mode='plot', **kwargs):
        try:
            col = COLORS[self.name.replace('#', '')]
        except KeyError:
            col = 'red'

        if mode == 'plot':
            ax.plot(self.Q, self.T_0, label=self.name.replace('#', ''), zorder=1, color=col, **kwargs)
        elif mode == 'scatter':
            p = ax.plot([n.nan], [n.nan], label=self.name, color=col, **kwargs)
            ax.scatter(self.Q, self.T_0, c=p[0].get_color(), s=1, zorder=1, **kwargs)