"""
Analyse simulation with ring topology.
Plot syncrhony variables against interconnection strength.

"""


import h5manager as hm
import tables
import matplotlib.pyplot as plt
import numpy as np


def main(dbfile):
    # Get the simulations
    db = tables.openFile(dbfile)
    simus = hm.get_first_level_groups(db.root)

    # Define some function to get specific result values
    def get_strength(simu):
        return hm.get_group_attr(simu, ('paramset', '_v_attrs', 'Common', 'inter_conn_strength', 0, 1))

    def get_mps(simu):
        return hm.get_group_attr(simu, ('results', '_v_attrs', 'MPS', 'whole'))

    def get_sts(simu):
        return hm.get_group_attr(simu, ('results', '_v_attrs', 'STS', 'whole'))

    def get_fftmax(simu):
        return hm.get_group_attr(simu, ('results', '_v_attrs', 'FFTMAX', 'mean'))

    get_indexes = (get_strength, get_mps, get_sts, get_fftmax)

    # Get simulation indexes for each simulation
    res_indexes = np.ndarray((len(simus), len(get_indexes)))
    for i_simu, simu in enumerate(simus):
        for i_index, get_index in enumerate(get_indexes):
            res_indexes[i_simu, i_index] = get_index(simu)
    # Sort index array on `strength` index
    indsort = res_indexes[:, 0].argsort()
    res_indexes = res_indexes[indsort, :]

    # Plot the res_indexes against interconnection strength
    plt.figure()
    plt.plot(res_indexes[:, 0], res_indexes[:, 1], '.', label="MPS (whole)")
    plt.plot(res_indexes[:, 0], res_indexes[:, 2], '.', label="STS (whole)")
    plt.plot(res_indexes[:, 0], res_indexes[:, 3], '.', label="FFTMAX (mean)")
    plt.legend()
    plt.show()

    # Get simulations to plot
    plot_simus = raw_input("Plot simulations #").split(' ')
    if plot_simus != ['']:
         for psimu in plot_simus:
             simu = simus[int(psimu)]
             hm.plot_simulation(simu)

    db.close()


if __name__ == '__main__':
    from sys import argv
    db = argv[1]
    main(db)
