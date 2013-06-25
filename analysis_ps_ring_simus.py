"""
Analyse simulation with ring topology.
Plot syncrhony variables against interconnection strength.

"""


import utils
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

    def get_peakdist(simu):
        peakdist = hm.get_group_attr(simu, ('results', '_v_attrs', 'peak_distances'))
        peakdist_means = utils.get_dict_values(peakdist, 2, "mean")
        peakdist_disps = utils.get_dict_values(peakdist, 2, "disp")
        return peakdist_means, peakdist_disps

    def get_peakdist_mean(simu):
        means, _ = get_peakdist(simu)
        return np.mean(means), np.std(means)

    def get_peakdist_mean_mean(simu):
        mean, _ = get_peakdist_mean(simu)
        return mean

    def get_peakdist_mean_disp(simu):
        _, disp = get_peakdist_mean(simu)
        return disp

    def get_peakdist_disp(simu):
        _, disps = get_peakdist(simu)
        return np.mean(disps), np.std(disps)
    
    def get_peakdist_disp_mean(simu):
        mean, _ = get_peakdist_disp(simu)
        return mean

    def get_peakdist_disp_disp(simu):
        _, disp = get_peakdist_disp(simu)
        return disp

    def get_spiking_rate(simu):
        spikes_it = hm.get_group_attr(simu, ('results', 'spikes_it')).read()
        nspikes = spikes_it.shape[1]
        simu_length = hm.get_group_attr(simu, ('paramset', '_v_attrs', 'Common', 'simu_length'))
        return nspikes/float(simu_length)

    get_indexes = (get_strength, get_mps, get_sts, get_fftmax,
                   get_peakdist_mean_mean, get_peakdist_mean_disp,
                   get_peakdist_disp_mean, get_peakdist_disp_disp,
                   get_spiking_rate)
    index_names = ("strength", "MPS", "STS", "FFTMAX",
                   "Peak Dist mean (mean)", "Peak Dist mean (disp)",
                   "Peak Dist disp (mean)", "Peak Dist disp (disp)",
                   "Spinking rate")

    # Get simulation indexes for each simulation
    res_indexes = np.ndarray((len(simus), len(get_indexes)))
    for i_simu, simu in enumerate(simus):
        for i_index, get_index in enumerate(get_indexes):
            res_indexes[i_simu, i_index] = get_index(simu)
    # Sort index array on `strength` index
    indsort = res_indexes[:, 0].argsort()
    res_indexes = res_indexes[indsort, :]
    # Sort simulations accordingly
    sorted_simulations = []
    for sorted_index in indsort:
        sorted_simulations.append(simus[sorted_index])

    # Plot the res_indexes against interconnection strength
    plt.figure()
    # Plot standard indexes
    data = {1: "MPS (whole)",
            2: "STS (whole)",
            3: "FFTMAX (mean) (Hz)",
            8: "Spiking rate (spikes/sec)"}
    for iplot in data:
        plt.plot(res_indexes[:, 0], res_indexes[:, iplot], '.', label=data[iplot])
    # Plot peak dist index
    plt.errorbar(res_indexes[:, 0], res_indexes[:, 4], yerr=res_indexes[:, 5],
                 fmt=".", label="Peak Dist (mean)")
    plt.errorbar(res_indexes[:, 0], res_indexes[:, 6], yerr=res_indexes[:, 7],
                 fmt=".", label="Peak Dist (disp)")
    # Add some plot information
    plt.xlabel("Interconnection strength")
    plt.legend()
    plt.show()

    # Get simulations to plot
    plot_simus = raw_input("Plot simulations #").split(' ')
    if plot_simus != ['']:
        for psimu in plot_simus:
            simu = sorted_simulations[int(psimu)]
            for index_name, index_fun in zip(index_names, get_indexes):
                print index_name, index_fun(simu)
            print "\n"
            hm.plot_simulation(simu)

    db.close()


if __name__ == '__main__':
    from sys import argv
    db = argv[1]
    main(db)
