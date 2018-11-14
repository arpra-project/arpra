
import numpy as np
import matplotlib.pyplot as plt

# SETUP
# %load_ext autoreload
# %autoreload 2
# from tools.arpra_mpfr_2d import arpra_mpfr_2d
# #####

def arpra_mpfr_2d (x, y, t, i_start, i_stop, path='./', ax_traj=None, ax_x=None, ax_y=None):

    with open(path + x, 'r') as xx_file, \
         open(path + y, 'r') as yy_file, \
         open(path + t, 'r') as tt_file:

        for i in range(i_start):
            xx_file.readline();
            yy_file.readline();
            tt_file.readline();

        xx = np.array([xx_file.readline() for i in range(i_start, i_stop)], dtype=np.float64)
        yy = np.array([yy_file.readline() for i in range(i_start, i_stop)], dtype=np.float64)
        tt = np.array([tt_file.readline() for i in range(i_start, i_stop)], dtype=np.float64)

    if ax_traj:
        # Plot (x, y) trajectory
        ax_traj.set_xlabel(x)
        ax_traj.set_ylabel(y)
        ax_traj.plot(xx, yy, color='b')#, marker='.')

    if ax_x:
        # Plot x through time
        ax_x.set_xlabel('time')
        ax_x.set_ylabel(x)
        ax_x.plot(tt, xx, color='b', label=x)

    if ax_y:
        # Plot y through time
        ax_y.set_xlabel('time')
        ax_y.set_ylabel(y)
        ax_y.plot(tt, yy, color='b', label=y)

    return
