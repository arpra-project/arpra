
import numpy as np
import matplotlib.pyplot as plt
from contextlib import ExitStack

# SETUP
# %load_ext autoreload
# %autoreload 2
# from experiments.experiment_2 import experiment_2
# #####


def experiment_2 ():

    x = 'nrn1_V_000.dat'
    y = 'nrn1_N_000.dat'
    t = 'time_000.dat'

    fig = plt.figure()
    ax_traj = fig.add_subplot(121)
    ax_x = fig.add_subplot(222)
    ax_y = fig.add_subplot(224, sharex=ax_x)
    fig.canvas.set_window_title('Experiment 2')
    for i in range(0, 100):
        path = 'experiment_2_out/i_' + str(i) + '/'
        arpra_mpfr_2d (x, y, t, 0, 1000, path, ax_traj, ax_x, ax_y)
    path = 'experiment_2_out/ascending/'
    arpra_mpfr_2d (x, y, t, 0, 1000, path, ax_traj, ax_x, ax_y, 'r')
    path = 'experiment_2_out/descending/'
    arpra_mpfr_2d (x, y, t, 0, 1000, path, ax_traj, ax_x, ax_y, 'r')
    plt.show()

    fig = plt.figure()
    ax_x = fig.add_subplot(211)
    ax_y = fig.add_subplot(212)    
    fig.canvas.set_window_title('Experiment 2')
    arpra_mpfr_mean_std (x, y, t, 100, 0, 1000, ax_x, ax_y)
    plt.show()

    return


def arpra_mpfr_2d (x, y, t, i_start, i_stop, path, ax_traj, ax_x, ax_y, col='b'):

    with ExitStack() as stack:
        x_file = stack.enter_context(open(path + x, 'r'))
        y_file = stack.enter_context(open(path + y, 'r'))
        t_file = stack.enter_context(open(path + t, 'r'))

        for i in range(i_start):
            x_file.readline();
            y_file.readline();
            t_file.readline();

        xx = np.array([x_file.readline() for i in range(i_start, i_stop)], dtype=np.float64)
        yy = np.array([y_file.readline() for i in range(i_start, i_stop)], dtype=np.float64)
        tt = np.array([t_file.readline() for i in range(i_start, i_stop)], dtype=np.float64)

    if ax_traj:
        # Plot (x, y) trajectory
        ax_traj.set_xlabel(x)
        ax_traj.set_ylabel(y)
        ax_traj.plot(xx, yy, color=col)#, marker='.')

    if ax_x:
        # Plot x through time
        ax_x.set_xlabel('time')
        ax_x.set_ylabel(x)
        ax_x.plot(tt, xx, color=col, label=x)

    if ax_y:
        # Plot y through time
        ax_y.set_xlabel('time')
        ax_y.set_ylabel(y)
        ax_y.plot(tt, yy, color=col, label=y)

    return


def arpra_mpfr_mean_std (x, y, t, n, i_start, i_stop, ax_x, ax_y):

    with ExitStack() as stack:

        x_files = [stack.enter_context(
            open('experiment_2_out/i_' + str(i) + '/' + x, 'r')) for i in range(n)]
        y_files = [stack.enter_context(
            open('experiment_2_out/i_' + str(i) + '/' + y, 'r')) for i in range(n)]
        t_files = [stack.enter_context(
            open('experiment_2_out/i_' + str(i) + '/' + t, 'r')) for i in range(n)]

        for i in range(i_start):
            for j in range(n):
                x_files[j].readline();
                y_files[j].readline();
                t_files[j].readline();

        xx = np.empty(n, dtype=np.float64)
        xx_mean = np.empty((i_stop - i_start), dtype=np.float64)
        xx_std = np.empty((i_stop - i_start), dtype=np.float64)
        yy = np.empty(n, dtype=np.float64)
        yy_mean = np.empty((i_stop - i_start), dtype=np.float64)
        yy_std = np.empty((i_stop - i_start), dtype=np.float64)
        tt = np.empty(n, dtype=np.float64)
        tt_mean = np.empty((i_stop - i_start), dtype=np.float64)

        for i in range(i_start, i_stop):
            for j in range(n): xx[j] = x_files[j].readline()
            xx_mean[i] = np.mean(xx)
            xx_std[i] = np.std(xx)
            for j in range(n): yy[j] = y_files[j].readline()
            yy_mean[i] = np.mean(yy)
            yy_std[i] = np.std(yy)
            for j in range(n): tt[j] = t_files[j].readline()
            tt_mean[i] = np.mean(tt)

    if ax_x:
        # Plot x through time
        ax_x.set_xlabel('time')
        ax_x.set_ylabel(x)
        #ax_x.plot(tt_mean, xx_mean, label=x)
        ax_x.plot(tt_mean, xx_std, label=x)

    if ax_y:
        # Plot y through time
        ax_y.set_xlabel('time')
        ax_y.set_ylabel(y)
        #ax_y.plot(tt_mean, yy_mean, label=y)
        ax_y.plot(tt_mean, yy_std, label=y)

    return


if __name__ == '__main__':
    experiment_2()
