
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

    # fig = plt.figure()
    # ax_traj = fig.add_subplot(121)
    # ax_x = fig.add_subplot(222)
    # ax_y = fig.add_subplot(224, sharex=ax_x)
    # fig.canvas.set_window_title('Experiment 2')
    # for i in range(0, 100):
    #     path = 'experiment_2_out/i_' + str(i) + '/'
    #     arpra_mpfr_plot (x, y, t, 0, 1000, path, ax_traj, ax_x, ax_y)
    # path = 'experiment_2_out/ascending/'
    # arpra_mpfr_plot (x, y, t, 0, 1000, path, ax_traj, ax_x, ax_y, 'r')
    # path = 'experiment_2_out/descending/'
    # arpra_mpfr_plot (x, y, t, 0, 1000, path, ax_traj, ax_x, ax_y, 'r')
    # plt.show()

    fig = plt.figure()
    ax_x1 = fig.add_subplot(411)
    ax_x2 = fig.add_subplot(412)
    ax_y1 = fig.add_subplot(413)
    ax_y2 = fig.add_subplot(414)
    fig.canvas.set_window_title('Experiment 2')
    arpra_mpfr_mean_std (x, y, t, 100, 0, 1000, ax_xm=ax_x1, ax_xs=ax_x2, ax_ym=ax_y1, ax_ys=ax_y2)
    arpra_centre_radius (x, y, t, 0, 1000, ax_xc=ax_x1, ax_xr=ax_x2, ax_yc=ax_y1, ax_yr=ax_y2)
    plt.show()
    return


def arpra_mpfr_plot (x, y, t, i_start, i_stop, path,
                     ax_traj=None, ax_x=None, ax_y=None, fmt='b'):

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
        ax_traj.plot(xx, yy, fmt)

    if ax_x:
        # Plot x through time
        ax_x.set_xlabel('time')
        ax_x.set_ylabel(x)
        ax_x.plot(tt, xx, fmt, label=x)

    if ax_y:
        # Plot y through time
        ax_y.set_xlabel('time')
        ax_y.set_ylabel(y)
        ax_y.plot(tt, yy, fmt, label=y)

    return


def arpra_mpfr_mean_std (x, y, t, n, i_start, i_stop,
                         ax_xm=None, ax_xm_fmt='b', ax_xs=None, ax_xs_fmt='b',
                         ax_ym=None, ax_ym_fmt='b', ax_ys=None, ax_ys_fmt='b'):

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

    if ax_xm:
        # Plot mean(x) through time
        ax_xm.set_xlabel('time')
        ax_xm.set_ylabel(x)
        ax_xm.plot(tt_mean, xx_mean, ax_xm_fmt, label=x)

    if ax_xs:
        # Plot std(x) through time
        ax_xs.set_xlabel('time')
        ax_xs.set_ylabel(x)
        ax_xs.plot(tt_mean, xx_std, ax_xs_fmt, label=x)

    if ax_ym:
        # Plot mean(y) through time
        ax_ym.set_xlabel('time')
        ax_ym.set_ylabel(y)
        ax_ym.plot(tt_mean, yy_mean, ax_ym_fmt, label=y)

    if ax_ys:
        # Plot std(y) through time
        ax_ys.set_xlabel('time')
        ax_ys.set_ylabel(y)
        ax_ys.plot(tt_mean, yy_std, ax_ys_fmt, label=y)

    return


def arpra_centre_radius (x, y, t, i_start, i_stop,
                         ax_xc=None, ax_xc_fmt='r', ax_xr=None, ax_xr_fmt='r',
                         ax_yc=None, ax_yc_fmt='r', ax_yr=None, ax_yr_fmt='r'):

    with ExitStack() as stack:
        xc_file = stack.enter_context(open(x + '_c.dat', 'r'))
        xr_file = stack.enter_context(open(x + '_r.dat', 'r'))
        yc_file = stack.enter_context(open(y + '_c.dat', 'r'))
        yr_file = stack.enter_context(open(y + '_r.dat', 'r'))
        tc_file = stack.enter_context(open(t + '_c.dat', 'r'))

        for i in range(i_start):
            xc_file.readline();
            xr_file.readline();
            yc_file.readline();
            yr_file.readline();
            tc_file.readline();

        xc = np.array([xc_file.readline() for i in range(i_start, i_stop)], dtype=np.float64)
        xr = np.array([xr_file.readline() for i in range(i_start, i_stop)], dtype=np.float64)
        yc = np.array([yc_file.readline() for i in range(i_start, i_stop)], dtype=np.float64)
        yr = np.array([yr_file.readline() for i in range(i_start, i_stop)], dtype=np.float64)
        tc = np.array([tc_file.readline() for i in range(i_start, i_stop)], dtype=np.float64)

    if ax_xc:
        # Plot x centre through time
        ax_xc.set_xlabel('time')
        ax_xc.set_ylabel(x)
        ax_xc.plot(tc, xc, ax_xc_fmt, label=x)

    if ax_xr:
        # Plot x radius through time
        ax_xr.set_xlabel('time')
        ax_xr.set_ylabel(x)
        ax_xr.plot(tc, xr, ax_xr_fmt, label=x)

    if ax_yc:
        # Plot y centre through time
        ax_yc.set_xlabel('time')
        ax_yc.set_ylabel(y)
        ax_yc.plot(tc, yc, ax_yc_fmt, label=y)

    if ax_yr:
        # Plot y radius through time
        ax_yr.set_xlabel('time')
        ax_yr.set_ylabel(y)
        ax_yr.plot(tc, yr, ax_yr_fmt, label=y)

    return


if __name__ == '__main__':
    experiment_2()
