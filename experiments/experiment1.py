
import numpy as np
import matplotlib.pyplot as plt

# SETUP
# %load_ext autoreload
# %autoreload 2
# from experiments.experiment1 import experiment1
# #####

def rad_sum (path, x, i_start, i_stop):

    with open(path + x, 'r') as xr_file:

        try:
            for i in range(i_start):
                xr_file.readline();

            xr = np.array([xr_file.readline() for i in range(i_start, i_stop)], dtype=np.float)
            return np.sum(xr)

        except ValueError:
            return np.inf


def experiment1 ():

    xr = 'nrn1_V_000_r.dat'
    yr = 'nrn1_N_000_r.dat'
    xr_sum_mat = np.zeros((25, 25), dtype=np.float)
    yr_sum_mat = np.zeros((25, 25), dtype=np.float)

    for i in range(1, 26):
        inps = i * 2

        for j in range(1, 26):
            freq = j * 2

            path = 'experiment_1_out/in_' + str(inps) + '_freq_' + str(freq) + '/'
            print(path)

            xr_sum_mat[i-1, j-1] = rad_sum (path, xr, 0, 1000)
            yr_sum_mat[i-1, j-1] = rad_sum (path, yr, 0, 1000)

    fig, ax = plt.subplots(1, 2)

    ax[0].set_title('V')
    ax[0].set_xlabel('Input frequency (Hz)')
    ax[0].set_ylabel('Input count')
    ax[0].matshow(np.log10(xr_sum_mat))



    #ax[0].set_xticks(np.arange(0, 25))
    #ax[0].set_yticks(np.arange(0, 25))
    #ax[0].set_xticklabels(np.arange(0, 50, 10))
    #ax[0].set_yticklabels(np.arange(0, 50, 10))



    ax[1].set_title('N')
    ax[1].set_xlabel('Input frequency (Hz)')
    ax[1].set_ylabel('Input count')
    ax[1].matshow(np.log10(yr_sum_mat))

    return


if __name__ == '__main__':
    experiment1()
