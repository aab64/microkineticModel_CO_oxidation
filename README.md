# Microkinetic model for CO oxidation on Cu

## Contents

This folder contains the following scripts for simulating kinetics of carbon monoxide oxidation on copper. 
1. solve_CO_oxidation_surface_only.m
2. solve_CO_oxidation_tanks_in_series.m
3. get_CO_oxidation_surface_odes.m
4. get_CO_oxidation_odes.m
5. get_CO_oxidation_rates.m
6. get_CO_oxidation_rate_constants.m
7. k_ads.m
8. k_des.m
9. k_arr.m
10. get_CO_oxidation_jac.m
11. solve_parameter_optimisation.m
12. parameter_optimisation.m 

## Running simulations

### 1. Surface only simulations

Files 1 and 3 are for surface kinetics only and assume constant bulk concentrations of CO, O2 and CO2.  

Run file 1 to solve the ODEs provided in file 3. 

### 2. Simulations with mass transport

Files 2 and 4 are for surface kinetics with a tanks-in-series model.  

Run file 2 to solve the ODEs provided in file 4.

### 3. Changing the parameters

Simulation parameters can be modified in the Settings section of files 1 and 2.  

Files 11 and 12 can be used to solve for parameters that meet certain conditions, for example to match experimentally observed phenomena. 

## Computing reaction rates 

Files 5-9 are used to set up the ODEs:
- File 5 computes the surface reaction rates, using file 6 to obtain values for the rate constants.
- File 6 computes the reaction rate constants, using files 7-9 to obtain expressions for adsorption, desorption and Arrhenius reaction.
- The reaction barriers in file 6 are taken from Farsig et al. [1].
- The expressions in files 7-9 are taken from Filot (2018) [2].

File 10 supplies the analytic Jacobian for the tanks-in-series ODEs. This is suggested to improve efficiency of *ode23s* in the Matlab documentation. However, sufficient stability has been experienced using *ode15s* which is significantly faster than *ode23s* for the current cases. 

[1] Falsig, H., Hvolbaek, B., Kristensen, I.S., Jiang, T., Bligaard, T., Christensen, C.H. and Norskov, J.K., 2008. Trends in the catalytic CO oxidation activity of nanoparticles. Angewandte Chemie International Edition, 47(26), pp.4835-4839. 

[2] Filot, I. A. W., 2018, Introduction to Microkinetic Modeling, Technische Universiteit Eindhoven, Eindhoven. 

---

Astrid Boje, 17 April 2020.
