
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

# DEBUG
# %load_ext autoreload
# %autoreload 2
# from tools.arpra_joint_range import arpra_joint_range
# x = 'nrn1_N_000'
# y = 'nrn1_V_000'
# x = 'syn_exc_R_000'
# y = 'syn_exc_S_000'
# #####

debug = 0
debug = 1

def arpra_joint_range (x, y, i_start, i_stop):

    fig = plt.figure()

    with open(x + '_c.dat', 'r') as xc_file, \
         open(x + '_r.dat', 'r') as xr_file, \
         open(x + '_s.dat', 'r') as xs_file, \
         open(x + '_d.dat', 'r') as xd_file, \
         open(y + '_c.dat', 'r') as yc_file, \
         open(y + '_r.dat', 'r') as yr_file, \
         open(y + '_s.dat', 'r') as ys_file, \
         open(y + '_d.dat', 'r') as yd_file, \
         open('time_000_c.dat', 'r') as tc_file, \
         open('time_000_r.dat', 'r') as tr_file:

        for i in range(i_start):
            xc_file.readline();
            xr_file.readline();
            xs_file.readline();
            xd_file.readline();
            yc_file.readline();
            yr_file.readline();
            ys_file.readline();
            yd_file.readline();
            tc_file.readline();
            tr_file.readline();

        xc = np.array([float(xc_file.readline()) for i in range(i_start, i_stop)])
        xr = np.array([float(xr_file.readline()) for i in range(i_start, i_stop)])
        xlo = xc - xr; xhi = xc + xr;

        yc = np.array([float(yc_file.readline()) for i in range(i_start, i_stop)])
        yr = np.array([float(yr_file.readline()) for i in range(i_start, i_stop)])
        ylo = yc - yr; yhi = yc + yr;

        tc = np.array([float(tc_file.readline()) for i in range(i_start, i_stop)])
        tr = np.array([float(tr_file.readline()) for i in range(i_start, i_stop)])

        if debug:
            # Plot trajectory
            ax1 = fig.add_subplot(121)
            ax1.plot(xc, yc)#, marker='.')
            ax1.set_xlabel(x)
            ax1.set_ylabel(y)

            # Plot through time
            ax2 = fig.add_subplot(222)
            ax2.plot(tc, xlo, color='b')
            ax2.plot(tc, xhi, color='b')
            ax2.plot(tc, xc, color='r', label=x)
            ax2.set_xlabel('time (msec)')
            ax2.set_ylabel(x)
            ax3 = fig.add_subplot(224, sharex=ax2)
            ax3.plot(tc, ylo, color='b')
            ax3.plot(tc, yhi, color='b')
            ax3.plot(tc, yc, color='r', label=y)
            ax3.set_xlabel('time (msec)')
            ax3.set_ylabel(y)

            return


        fig, ax = plt.subplots(1)
        ax.set_xlabel(x)
        ax.set_ylabel(y)

        # Plot interval regions
        cmap = plt.get_cmap('hsv')
        for i in range(0, (i_stop - i_start)):
            print('Iteration ' + str(i_start + i))

            # CHANGE xc, xr, yc, yr TO NUMPY SCALARS

            #pos = ((xc[i] - xr[i]), (yc[i] - yr[i]))
            pos = (-xr[i], -yr[i])
            x_width = 2 * xr[i]
            y_width = 2 * yr[i]
            color = cmap((5 * i) % cmap.N)
            box = patches.Rectangle(pos, x_width, y_width, facecolor='none', edgecolor=color)
            ax.add_patch(box)
            ax.autoscale()
            plt.pause(0.01)

        return



        # DEBUG
        #p = yc_file.readline()
        #print(p)
        #print(float(p))
        #print([float(pp) for pp in p.split()])
        # #####


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
