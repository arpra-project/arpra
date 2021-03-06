
import numpy as np
import matplotlib.pyplot as plt

# SETUP
# %load_ext autoreload
# %autoreload 2
# from experiments.experiment_1 import experiment_1
# #####

def rad_sum (path, x, i_start, i_stop):

    with open(path + x, 'r') as xr_file:

        try:
            for i in range(i_start):
                xr_file.readline();

            xr = np.array([xr_file.readline() for i in range(i_start, i_stop)], dtype=np.float64)
            return np.sum(xr)

        except ValueError:
            return np.inf


def experiment_1 ():

    xr = 'nrn1_V_000_r.dat'
    yr = 'nrn1_N_000_r.dat'
    xr_sum_mat = np.zeros((26, 26), dtype=np.float64)
    yr_sum_mat = np.zeros((26, 26), dtype=np.float64)

    for i in range(0, 26):
        inps = i * 2

        for j in range(0, 26):
            freq = j * 2

            path = 'experiment_1_out/in_' + str(inps) + '_freq_' + str(freq) + '/'
            print(path)

            xr_sum_mat[i, j] = rad_sum (path, xr, 0, 1000)
            yr_sum_mat[i, j] = rad_sum (path, yr, 0, 1000)

    fig, ax = plt.subplots(1, 2)

    ax[0].matshow(np.log10(xr_sum_mat))
    ax[0].set_xticks(np.arange(0, 26, 5))
    ax[0].set_yticks(np.arange(0, 26, 5))
    ax[0].set_xticklabels(np.arange(0, 51, 10))
    ax[0].set_yticklabels(np.arange(0, 51, 10))
    ax[0].set_xlabel('Input frequency (Hz)')
    ax[0].set_ylabel('Input count')
    ax[0].set_title('V')

    ax[1].matshow(np.log10(yr_sum_mat))
    ax[1].set_xticks(np.arange(0, 26, 5))
    ax[1].set_yticks(np.arange(0, 26, 5))
    ax[1].set_xticklabels(np.arange(0, 51, 10))
    ax[1].set_yticklabels(np.arange(0, 51, 10))
    ax[1].set_xlabel('Input frequency (Hz)')
    ax[1].set_ylabel('Input count')
    ax[1].set_title('N')

    plt.show()
    return


if __name__ == '__main__':
    experiment_1()
