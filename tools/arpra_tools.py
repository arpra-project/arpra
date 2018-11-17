
from contextlib import ExitStack
from itertools import islice
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

# SETUP
# %load_ext autoreload
# %autoreload 2
# from tools.arpra_tools import *
# #####


def arpra_plot (x, y, t, i_start=None, i_stop=None, i_step=None, path='.'):

    with ExitStack() as stack:
        xc_file = stack.enter_context(open(path + '/' + x + '_c.dat', 'r'))
        xr_file = stack.enter_context(open(path + '/' + x + '_r.dat', 'r'))
        yc_file = stack.enter_context(open(path + '/' + y + '_c.dat', 'r'))
        yr_file = stack.enter_context(open(path + '/' + y + '_r.dat', 'r'))
        tc_file = stack.enter_context(open(path + '/' + t + '_c.dat', 'r'))
        tr_file = stack.enter_context(open(path + '/' + t + '_r.dat', 'r'))

        xc = np.genfromtxt(islice(xc_file, i_start, i_stop, i_step), dtype=np.float64)
        xr = np.genfromtxt(islice(xr_file, i_start, i_stop, i_step), dtype=np.float64)
        xlo = xc - xr
        xhi = xc + xr;

        yc = np.genfromtxt(islice(yc_file, i_start, i_stop, i_step), dtype=np.float64)
        yr = np.genfromtxt(islice(yr_file, i_start, i_stop, i_step), dtype=np.float64)
        ylo = yc - yr
        yhi = yc + yr;

        tc = np.genfromtxt(islice(tc_file, i_start, i_stop, i_step), dtype=np.float64)
        tr = np.genfromtxt(islice(tr_file, i_start, i_stop, i_step), dtype=np.float64)

    fig = plt.figure()
    fig.canvas.set_window_title('arpra_plot')

    # Plot trajectory
    ax1 = fig.add_subplot(121)
    ax1.set_xlabel(x)
    ax1.set_ylabel(y)
    ax1.plot(xc, yc)

    # Plot x through time
    ax2 = fig.add_subplot(222)
    ax2.set_xlabel('time')
    ax2.set_ylabel(x)
    ax2.plot(tc, xlo, 'b')
    ax2.plot(tc, xhi, 'b')
    #ax2.plot(tc, xc, 'r', label=x)

    # Plot y through time
    ax3 = fig.add_subplot(224, sharex=ax2)
    ax3.set_xlabel('time')
    ax3.set_ylabel(y)
    ax3.plot(tc, ylo, 'b')
    ax3.plot(tc, yhi, 'b')
    #ax3.plot(tc, yc, 'r', label=y)

    return


def arpra_joint_range (x, y, t, i_start=None, i_stop=None, i_step=None, path='.'):

    with ExitStack() as stack:
        xc_file = stack.enter_context(open(path + '/' + x + '_c.dat', 'r'))
        xr_file = stack.enter_context(open(path + '/' + x + '_r.dat', 'r'))
        xs_file = stack.enter_context(open(path + '/' + x + '_s.dat', 'r'))
        xd_file = stack.enter_context(open(path + '/' + x + '_d.dat', 'r'))
        yc_file = stack.enter_context(open(path + '/' + y + '_c.dat', 'r'))
        yr_file = stack.enter_context(open(path + '/' + y + '_r.dat', 'r'))
        ys_file = stack.enter_context(open(path + '/' + y + '_s.dat', 'r'))
        yd_file = stack.enter_context(open(path + '/' + y + '_d.dat', 'r'))
        tc_file = stack.enter_context(open(path + '/' + t + '_c.dat', 'r'))
        tr_file = stack.enter_context(open(path + '/' + t + '_r.dat', 'r'))

        xc = np.genfromtxt(islice(xc_file, i_start, i_stop, i_step), dtype=np.float64)
        xr = np.genfromtxt(islice(xr_file, i_start, i_stop, i_step), dtype=np.float64)
        xlo = xc - xr
        xhi = xc + xr;

        yc = np.genfromtxt(islice(yc_file, i_start, i_stop, i_step), dtype=np.float64)
        yr = np.genfromtxt(islice(yr_file, i_start, i_stop, i_step), dtype=np.float64)
        ylo = yc - yr
        yhi = yc + yr;

        tc = np.genfromtxt(islice(tc_file, i_start, i_stop, i_step), dtype=np.float64)
        tr = np.genfromtxt(islice(tr_file, i_start, i_stop, i_step), dtype=np.float64)

    fig = plt.figure()
    fig.canvas.set_window_title('arpra_plot')

    ax1 = fig.add_subplot(111)
    ax1.set_xlabel(x)
    ax1.set_ylabel(y)

    # Plot interval regions
    cmap = plt.get_cmap('hsv')
    for i in range(0, (i_stop - i_start)):
        print('Iteration ' + str(i_start + i))

        #pos = ((xc[i] - xr[i]), (yc[i] - yr[i]))
        pos = (-xr[i], -yr[i])
        x_width = 2 * xr[i]
        y_width = 2 * yr[i]
        color = cmap((5 * i) % cmap.N)
        box = patches.Rectangle(pos, x_width, y_width, facecolor='none', edgecolor=color)
        ax1.add_patch(box)
        ax1.autoscale()
        #plt.pause(0.01)

    return


    ## Plot affine regions
    # maxterms = 15;
    # e = zeros(maxterms, 2^maxterms);
    # for i = 1:2^maxterms
    #     e(:, i) = double(bitget(i - 1, 1:maxterms))';
    # end
    # e(e == 0) = -1;
    #
    #for i = i_start:i_stop
    #    disp(num2str(i));
    #
    #    #[xc, ~, err] = sscanf(fgetl(xc_file), '%f');
    #    if ~isempty(err); return; end;
    #    #[xr, ~, err] = sscanf(fgetl(xr_file), '%f');
    #    if ~isempty(err); return; end;
    #    [xs, ~, err] = sscanf(fgetl(xs_file), '%u');
    #    if ~isempty(err); return; end;
    #    [xd, ~, err] = sscanf(fgetl(xd_file), '%f');
    #    if ~isempty(err); return; end;
    #    #[yc, ~, err] = sscanf(fgetl(yc_file), '%f');
    #    if ~isempty(err); return; end;
    #    #[yr, ~, err] = sscanf(fgetl(yr_file), '%f');
    #    if ~isempty(err); return; end;
    #    [ys, ~, err] = sscanf(fgetl(ys_file), '%u');
    #    if ~isempty(err); return; end;
    #    [yd, ~, err] = sscanf(fgetl(yd_file), '%f');
    #    if ~isempty(err); return; end;
    #
    #    us = union(xs, ys);
    #    if isrow(us)
    #        us = us';
    #    end
    #    terms = size(us, 1);
    #
    #    ix = ismember(us, xs);
    #    xxd = zeros(1, terms);
    #    xxd(ix) = xd;
    #    iy = ismember(us, ys);
    #    yyd = zeros(1, terms);
    #    yyd(iy) = yd;
    #
    #    # Get all permutations of noise symbol extremities
    #    xx = zeros(2^terms, 1);
    #    yy = zeros(2^terms, 1);
    #    for j = 1:2^terms
    #        e = double(bitget(j - 1, 1:terms))';
    #        e(e == 0) = -1;
    #        #xx(j) = xc + xxd * e;
    #        #yy(j) = yc + yyd * e;
    #        xx(j) = xxd * e;
    #        yy(j) = yyd * e;
    #
    #    # Get all permutations of noise symbol extremities
    #    #xx = xc + xxd * e(1:terms, 1:2^terms);
    #    #yy = yc + yyd * e(1:terms, 1:2^terms);
    #    xx = xxd * e(1:terms, 1:2^terms);
    #    yy = yyd * e(1:terms, 1:2^terms);
    # 
    #    k = convhull(xx, yy);
    #    plot(xx(k), yy(k));
    #    drawnow;
