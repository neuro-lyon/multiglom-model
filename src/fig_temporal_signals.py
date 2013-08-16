# -*- coding:utf-8 -*-
import tables as tb
import matplotlib.pyplot as plt
import numpy as np

def main():
    max_time = 2000
    time = np.linspace(0, 3, 6000)

    # Select Data
    data = tb.openFile("data/db50_single_glom.h5")
    sim_results = data.root.results

    # Create figure
    plt.figure()

    # Get signals
    # INPUT
    input_sig = sim_results.input.read()[0]
    plt.subplot(4, 1, 1)
    plt.plot(time[:max_time], input_sig[:max_time])
    plt.ylabel('g_input ($S m^{-2}$)')

    # MEMB POT
    memb_pot = sim_results.single_memb_pot.read()
    plt.subplot(4, 1, 2)
    plt.plot(time[:max_time], memb_pot[:max_time])
    plt.ylabel('V (1 mitrale)')

    # S
    s_sig = sim_results.s_granule.read()[0]
    plt.subplot(4, 1, 3)
    plt.plot(time[:max_time], s_sig[:max_time])
    plt.ylabel('s (granule)')

    # S_SYN_SELF
    sss_sig = sim_results.s_syn_self.read()[0]
    plt.subplot(4, 1, 4)
    plt.plot(time[:max_time], sss_sig[:max_time])
    plt.ylabel('s_syn_self')
    plt.xlabel("Temps (s)")

    # Plot the figure
    plt.show()

    # Quit
    data.close()

if __name__ == '__main__':
    main()
