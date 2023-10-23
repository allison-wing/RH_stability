# RH_stability
Analysis scripts and post-processed model data for Wing and Singh (Control of Stability and Relative Humidity in the Radiative-Convective Equilibrium Model Intercomparison Project)

Includes...

## Plotting Codes
- pub_plot_CM1_theory_hur_cape.m: Figures 3-4
- pub_plot_cape.m: Figure 2
- pub_plot_dcape_decompose.m: Figure 9 and S6
- pub_plot_prof_ta_moistadiabat_RH.m: Figure 1
- pub_plot_profiles_ZBP.m: Figures S1-S4
- pub_plot_theory_hur_cape_Romps.m Figures 5-8 and S5

## Functions
- matlab_packages: A subset of needed codes from Martin Singh's matlab packages (https://gitlab.com/martinsingh/matlab_packages)
- calc_PEproxy.m
- calc_precipwaterpath.m
- calculate_CAPE.m
- calculate_CAPE_derivatives.m
- calculate_CAPE_theory.m
- calculate_ZBP.m
- calculate_ZBP_lapse_rate.m
- find_PE_epsilon.m
- fn_plot_CAPE_and_RCE_RCE.m

## Derived Data Files
- precippath.mat: Precipitating condensate path in RCEMIP RCE_small simulations
- PEproxy_cwipwi.mat: Precipitation efficiency proxy in RCEMIP RCE_small simulations
- RH_stab_2-5.mat: Lower-tropospheric relative humidity, CAPE, theory-implied entrainment, theory-implied precipitation efficiency, and associated variables for the RCEMIP RCE_small simulations.

## Model Output Files
- CM1_expts: Parameter perturbation experiments with CM1
- RCEMIP_Table A1: Time and domain-mean RCEMIP 0D data from the RCE_small simulations. Table A1 from Wing et al. 2020. Can also be found at http://hdl.handle.net/21.14101/d4beee8e-6996-453e-bbd1-ff53b6874c0e.
- RCEMIP_1Dprofiles_RCE_small.zip: Time- and domain-mean RCEMIP 1D profile data from the RCE_small simulations. Can also be found at http://hdl.handle.net/21.14101/d4beee8e-6996-453e-bbd1-ff53b6874c0e.
