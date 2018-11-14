
import numpy as np
import matplotlib.pyplot as plt

# SETUP
# %load_ext autoreload
# %autoreload 2
# from tools.arpra_mpfi_2d import arpra_mpfi_2d
# #####

def arpra_mpfi_2d (x, y, t, i_start, i_stop):

    with open(x + '.dat', 'r') as xx_file, \
         open(y + '.dat', 'r') as yy_file, \
         open(t + '.dat', 'r') as tt_file:

        for i in range(i_start):
            xx_file.readline();
            yy_file.readline();
            tt_file.readline();

        xx = np.array([xx_file.readline().split() for i in range(i_start, i_stop)], dtype=np.float64)
        yy = np.array([yy_file.readline().split() for i in range(i_start, i_stop)], dtype=np.float64)
        tt = np.array([tt_file.readline().split() for i in range(i_start, i_stop)], dtype=np.float64)

    fig = plt.figure()
    fig.canvas.set_window_title('arpra_mpfi_2d')

    # Plot trajectory
    ax1 = fig.add_subplot(121)
    ax1.set_xlabel(x)
    ax1.set_ylabel(y)
    ax1.plot(np.mean(xx, 1), np.mean(yy, 1))#, marker='.')

    # Plot through time
    ax2 = fig.add_subplot(222)
    ax2.set_xlabel('time')
    ax2.set_ylabel(x)
    ax2.plot(np.mean(tt, 1), xx, color='b', label=x)
    ax3 = fig.add_subplot(224, sharex=ax2)
    ax3.set_xlabel('time')
    ax3.set_ylabel(y)
    ax3.plot(np.mean(tt, 1), yy, color='b', label=y)

    return
