"""
Analyse the distribution of peak distances to make appear the bistability.

To be used for two glomeruli simulations.

"""


import tables
import matplotlib.pyplot as plt
import analysis
import utils
import h5manager as hm


def histogram_peaks(simulation):
    dists = get_peak_dists(simulation)
    plt.figure()
    plt.hist(dists)


def get_peak_dists(simulation):
    # First, check that this is two glomeruli simulation
    n_glom = simulation.paramset._v_attrs['Common']['N_subpop']
    assert n_glom == 2, "Number of glomeruli must be 2."
    # Get the interesting signals : s_syn_self
    s_syn_self = simulation.results.s_syn_self.read()
    # Return the peak distances, with distances as angles
    return analysis.get_directional_distances(s_syn_self[0], s_syn_self[1])


def get_ind_rate_strength(simulations):
    """Return a 2 D list of rate and strength"""
    res = []
    for simu in simulations:
        values = []
        rate = simu.paramset._v_attrs['Common']['inter_conn_rate'][0][1]
        strength = simu.paramset._v_attrs['Common']['inter_conn_strength'][0][1]
        res.append([rate, strength])
    res.sort()
    return res


def get_1d_ind_simu(ind, width):
    """Return 1d index given `ind` like 5,3"""
    x, y = ind.split(",")
    x, y = int(x), int(y)
    return utils.to1d(x, y, width)


def main(dbfile):
    # Get the db
    db = tables.openFile(dbfile)
    simulations = hm.get_first_level_groups(db.root)
    # Sort simulations according to `rate` and `strength`
    ind_rate_strength = get_ind_rate_strength(simulations)
    n_rate = len(set([i[0] for i in ind_rate_strength]))

    # Selected simus are input like "1,2 1,5 2,10"
    selected_simus = raw_input("Simulations #").split(' ')
    selected_simus = [get_1d_ind_simu(s, n_rate) for s in selected_simus]

    # Plot the histogram for each selected simulation
    for i_sim in selected_simus:
        histogram_peaks(simulations[i_sim])

    plt.show()
    db.close()


if __name__ == '__main__':
    from sys import argv
    main(argv[1])  # argument must be a HDF5 file with simulations
