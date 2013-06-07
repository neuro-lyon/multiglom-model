-   `db30x30_new_gamma_non_homeo.h5` (2013-06-07)

    Gamma simulation without homeostatic connection.

    Parameters changed were:
    - interco. rate: 0 to 1, 30 steps
    - interco. strength: 0 to 0.1, 30 steps

-   `db10x10_new_gamma_homeo.h5`

    Same as `db10x10_new_gamma_non_homeo.h5` but with homeostatic connections.
    The new gamma parameters (git rev: `81a48e`) were used.


-   `db10x10_new_gamma_non_homeo.h5`

    Gamma simulation explore the parameter space.
    Take into account the *new gamma parameters*
    (as of git rev `81a48e4479e337d30126cb2838800a3a3425a52d`).

    Parameters were: interco. strength from 0 to 1 (10 steps),
    interco. rate from 0 to 1 (10 steps), in the gamma regime.


-   `db10x10_gamma_homeo.h5`: (2013-05-29)

    Same as `db10x10_gamma_non_homeo.h5` but with homeostatic connections.


-   `db10x10_gamma_non_homeo.h5`: (2013-05-29)

    Gamma simulation to explore the parameter space.

    Parameters were: interco. strength from 0 to 1 (10 steps),
    interco. rate from 0 to 1 (10 steps), in the gamma regime.


-   `db30x30_beta_non_homeo.h5`: (2013-05-28)

    Same as `db30x30_beta_homeostasis.h5` but without homeostasis.


-   `db30x30_beta_homeostasis.h5`: (2013-05-27)

    Homeostasy simulations to be analysed with the new peak distance index.
    Simulations were done for *3 seconds*.

    Parameters of interest were : interconnection strength from 0 to 0.1 (30 steps),
    interconnection rate from 0 to 1 (30 steps).


-   `lite_sim10x10beta.h5`: (2013-05-15)

    Testing the storage space improvement (it worked!). Went from 4.7 GB to 220 MB
    for 10x10 simulations.
 
    Simulation were done for 2 s in the beta regime.
    Interconnection rate from 0 to 1 (10x), interconnection strength from 0
    to *0.01* (10x).


-   `2figs-poster-homeostasy.h5`: (2013-05-13)

    Try to reproduce the results of Nicolas for the homeostasy.
    Interconnection rate fixed at 90%, strength from 0.02 to 0.1.


-   `db10x10beta.h5`:

    First big run done, for a total of 100 simulations.
    Two parameters were changed: interconnection strength (from 0 to 1, 10 values)
    and the interconnection rate (from 0 to 1, 10 values).


-   `db30x30beta.h5`:

    Second big run done, this time for a total of 900 simulations.
    Same thing as db10x10beta.h5 but this time the parameters were from 0 to 1,
    with 30 values each.


-   `db30x30gamma.h5`:


-   `db30x30_two_glom_beta_new_ps_interco_strength0_1_interco_rate0_1.h5`

    Another big run in the Beta regime and varying the interconnection rate and
    strength from 0 to 1 with 30 values each.

    This big run uses the new parameter definition for the Beta regime
    (git rev: `f51ff4b361decf32589e7193c5bc5b6b8b0078eb`)


-   `db5x5beta.h5`:


-   `redo_simu_fft.h5`:


-   `redo_simu_fft_with_seed.h5`:
