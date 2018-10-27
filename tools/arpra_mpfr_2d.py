
import numpy as np
import matplotlib.pyplot as plt

# SETUP
# %load_ext autoreload
# %autoreload 2
# from tools.arpra_mpfr_2d import arpra_mpfr_2d
# #####

def arpra_mpfr_2d (x, y, t, i_start, i_stop, ax_traj, ax_x, ax_y, path='./'):

    with open(path + x, 'r') as xx_file, \
         open(path + y, 'r') as yy_file, \
         open(path + t, 'r') as tt_file:

        for i in range(i_start):
            xx_file.readline();
            yy_file.readline();
            tt_file.readline();

        xx = np.array([xx_file.readline() for i in range(i_start, i_stop)], dtype=np.float)
        yy = np.array([yy_file.readline() for i in range(i_start, i_stop)], dtype=np.float)
        tt = np.array([tt_file.readline() for i in range(i_start, i_stop)], dtype=np.float)

    # Plot (x, y) trajectory
    ax_traj.set_xlabel(x)
    ax_traj.set_ylabel(y)
    ax_traj.plot(xx, yy, color='b')#, marker='.')

    # Plot x through time
    ax_x.set_xlabel('time')
    ax_x.set_ylabel(x)
    ax_x.plot(tt, xx, color='b', label=x)

    # Plot y through time
    ax_y.set_xlabel('time')
    ax_y.set_ylabel(y)
    ax_y.plot(tt, yy, color='b', label=y)

    return


def main ():

    fig = plt.figure()
    fig.canvas.set_window_title('arpra_mpfr_2d')

    ax_traj = fig.add_subplot(121)
    ax_x = fig.add_subplot(222)
    ax_y = fig.add_subplot(224, sharex=ax_x)

    x = 'nrn1_V_000.dat'
    y = 'nrn1_N_000.dat'
    t = 'time_000.dat'

    for i in range(0, 100):
        path = 'experiment_2_out/i_' + str(i) + '/'
        arpra_mpfr_2d (x, y, t, 0, 1000, ax_traj, ax_x, ax_y, path)

    plt.show()
    return


if __name__ == '__main__':
    main()
