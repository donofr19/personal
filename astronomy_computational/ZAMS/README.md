Main Sequence
=============

Project directory for a model of fully convective, low-mass main-sequence star.

Contents
--------

0. `README.md`: this file
1. `astro_const.py`: module containing physical constants, uses `astropy`
2. `eos.py`: starter code for the adiabatic equation of state and mean molecular weight
3. `reactions.py`: starter code for the proton-proton heating rate per mass
4. `structure.py`: starter code to integrate stellar structure equations
5. `structure_for_main.py`: updating central_thermal and teff for convenient handling in main files
6. `zams.py`: starter code for computing the zero-age main sequence, including tabulated values of the surface effective temperature
7. `Main_pp1.py`and `Main_pp10_5.py`: generates plots for the stellar relationships of a zero-age main sequence star
8. `test_convergence.py`: test convergence of delta, eta, xi parameters
9. `rootfind.py`: contains rootfind routine
10. `testing.py`: routine for calling the test suite from an IDE
11. `tests`: battery of acceptance tests
    1. `numcheck.py`: holds routine for computing relative errors
    2. `test_const.py`: tests selected constants against IAU, CODATA
    3. `test_eos.py`: evaluates adiabat and mean molecular weight
    4. `test_reactions.py`: evaluates the reaction rate module
    5. `test_structure.py`: evaluates computation of central P, rho, T in `structure`
    6. `test_zams.py`: evaluates interpolation of effective temperature and surface luminosity in `zams`
12. `ZAMS_Project_Report.pdf`: final explanatory report of project

Testing
-------

To perform the suite of unit tests, you may either run `testing.py` using an IDE such as Spyder or by executing `python -m pytest` from the command line.

At present all tests are skipped except for `test_const`. The following table lists the routines that are to be tested.  As you complete each routine, go into the appropriate test file and delete the line starting with `@pytest.mark.skip` that precedes the test routine.

| Routine | Test File | Test routine |
| :------- | :--------- | :------------ |
| `mean_molecular_weight` | `test_eos.py` | `test_chem` |
| `get_rho_and_T` | `test_eos.py` | `test_adiabat` |
| `pp_rate` | `test_reactions.py` | `test_pp` |
| `central_thermal` | `test_structure.py` |  `test_central_thermal` |
| `Teff` | `test_zams.py` | `test_Teff` |
| `surface_luminosity` | `test_zams.py` | `test_surface_luminosity` |

Instructions
------------

1. Run `astro_const.py`
2. Run `eos.py`
3. Run `ode.py`
4. Run `reactions.py`
5. Run `zams.py`
6. Run `structure.py`
7. Run `structure_for_main.py`
8. Run `rootfind.py`
9. Run `test_convergence.py`
10. Run `Main_pp1.py`
11. Run `Main_pp10_5.py`
