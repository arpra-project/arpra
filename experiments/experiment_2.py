
import numpy as np
import matplotlib.pyplot as plt
from contextlib import ExitStack
from itertools import islice

# SETUP
# %load_ext autoreload
# %autoreload 2
# from experiments.experiment_2 import experiment_2
# #####


def experiment_2 ():

    v_name = 'nrn1_V_000'
    n_name = 'nrn1_N_000'
    t_name = 'time_000'

    samples = 100
    v = [None] * samples
    n = [None] * samples
    t = [None] * samples

    fig = plt.figure()
    fig.canvas.set_window_title('Experiment 2')
    ax_traj = fig.add_subplot(121)
    ax_x = fig.add_subplot(222)
    ax_y = fig.add_subplot(224, sharex=ax_x)

    # for i in range(samples):
    #     [v[i], n[i], t[i]] = arpra_mpfr_read(path='experiment_2_out/i_' + str(i))
    # [v_asc, n_asc, t_asc] = arpra_mpfr_read(path='experiment_2_out/ascending')
    # [v_desc, n_desc, t_desc] = arpra_mpfr_read(path='experiment_2_out/descending')
    # [v_hp, n_hp, t_hp] = arpra_mpfr_read(path='experiment_2_out/high_prec')
    [v_arp_c, v_arp_r, n_arp_c, n_arp_r, t_arp_c, t_arp_r] = arpra_read(path='experiment_2_out/arpra')

    # v_mean = np.mean(v, axis=0)
    # v_std = np.std(v, axis=0)
    # n_mean = np.mean(n, axis=0)
    # n_std = np.std(n, axis=0)

    # for i in range(samples):
    #     ax_x.plot(t[i], v[i], 'b')
    #     ax_y.plot(t[i], n[i], 'b')
    # ax_x.plot(t_asc, v_asc, 'y')
    # ax_y.plot(t_asc, n_asc, 'y')
    # ax_x.plot(t_desc, v_desc, 'y')
    # ax_y.plot(t_desc, n_desc, 'y')
    # ax_x.plot(t_hp, v_hp, 'k')
    # ax_y.plot(t_hp, n_hp, 'k')
    # ax_traj.plot(v_hp, n_hp, 'b')

    v_arp_lo = v_arp_c - v_arp_r
    v_arp_hi = v_arp_c + v_arp_r
    n_arp_lo = n_arp_c - n_arp_r
    n_arp_hi = n_arp_c + n_arp_r
    ax_x.plot(t_arp_c, v_arp_lo, 'r')
    ax_x.plot(t_arp_c, v_arp_hi, 'r')
    ax_y.plot(t_arp_c, n_arp_lo, 'r')
    ax_y.plot(t_arp_c, n_arp_hi, 'r')

    # v_mean_std_lo = v_mean - v_std
    # v_mean_std_hi = v_mean + v_std
    # n_mean_std_lo = n_mean - n_std
    # n_mean_std_hi = n_mean + n_std
    # ax_x.plot(t_arp_c, v_mean_std_lo, 'r')
    # ax_x.plot(t_arp_c, v_mean_std_hi, 'r')
    # ax_y.plot(t_arp_c, n_mean_std_lo, 'r')
    # ax_y.plot(t_arp_c, n_mean_std_hi, 'r')

    plt.show()
    return


    fig = plt.figure()
    fig.canvas.set_window_title('Experiment 2')
    ax_v1 = fig.add_subplot(411)
    ax_v2 = fig.add_subplot(412, sharex=ax_v1)
    ax_n1 = fig.add_subplot(413, sharex=ax_v1)
    ax_n2 = fig.add_subplot(414, sharex=ax_v1)

    plt.show()
    return


def arpra_read (path='.', x='nrn1_V_000', y='nrn1_N_000', t='time_000',
                i_start=None, i_stop=None, i_step=None):

    with ExitStack() as stack:
        xc_file = stack.enter_context(open(path + '/' + x + '_c.dat', 'r'))
        xr_file = stack.enter_context(open(path + '/' + x + '_r.dat', 'r'))
        yc_file = stack.enter_context(open(path + '/' + y + '_c.dat', 'r'))
        yr_file = stack.enter_context(open(path + '/' + y + '_r.dat', 'r'))
        tc_file = stack.enter_context(open(path + '/' + t + '_c.dat', 'r'))
        tr_file = stack.enter_context(open(path + '/' + t + '_r.dat', 'r'))

        xc = np.genfromtxt(islice(xc_file, i_start, i_stop, i_step), dtype=np.float64)
        xr = np.genfromtxt(islice(xr_file, i_start, i_stop, i_step), dtype=np.float64)
        yc = np.genfromtxt(islice(yc_file, i_start, i_stop, i_step), dtype=np.float64)
        yr = np.genfromtxt(islice(yr_file, i_start, i_stop, i_step), dtype=np.float64)
        tc = np.genfromtxt(islice(tc_file, i_start, i_stop, i_step), dtype=np.float64)
        tr = np.genfromtxt(islice(tr_file, i_start, i_stop, i_step), dtype=np.float64)

    return [xc, xr, yc, yr, tc, tr]


def arpra_mpfr_read (path='.', x='nrn1_V_000', y='nrn1_N_000', t='time_000',
                     i_start=None, i_stop=None, i_step=None):

    with ExitStack() as stack:
        x_file = stack.enter_context(open(path + '/' + x + '.dat', 'r'))
        y_file = stack.enter_context(open(path + '/' + y + '.dat', 'r'))
        t_file = stack.enter_context(open(path + '/' + t + '.dat', 'r'))

        xx = np.genfromtxt(islice(x_file, i_start, i_stop, i_step), dtype=np.float64)
        yy = np.genfromtxt(islice(y_file, i_start, i_stop, i_step), dtype=np.float64)
        tt = np.genfromtxt(islice(t_file, i_start, i_stop, i_step), dtype=np.float64)

    return [xx, yy, tt]


def arpra_mpfr_plot (x, y, t, path='.', fmt='b',
                     ax_traj=None, ax_x=None, ax_y=None,
                     i_start=None, i_stop=None, i_step=None):

    with ExitStack() as stack:
        x_file = stack.enter_context(open(path + '/' + x + '.dat', 'r'))
        y_file = stack.enter_context(open(path + '/' + y + '.dat', 'r'))
        t_file = stack.enter_context(open(path + '/' + t + '.dat', 'r'))

        xx = np.genfromtxt(islice(x_file, i_start, i_stop, i_step), dtype=np.float64)
        yy = np.genfromtxt(islice(y_file, i_start, i_stop, i_step), dtype=np.float64)
        tt = np.genfromtxt(islice(t_file, i_start, i_stop, i_step), dtype=np.float64)

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


if __name__ == '__main__':
    experiment_2()
