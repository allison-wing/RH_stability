# matlab_packages

A set of useful packages for MATLAB

Each directory that begins with "+" contains a package. If this package is placed within the
MATLAB path, the functions within it can be accessed using the notation:

  >> package_name.function_name

Packages can contain subpackages which are referenced as might be expected:

  >> package_name.subpackage_name.function_name





Package list:

+atm
+ODE


#########################################################################################################




+atm		################################################################

   Some functions to calculate thermodynamic quantities useful in the atmospheric sciences.
   For instance, calculating saturation humidities, calculating a moist adiabat, and calculating entropy.

   calculate_MSE                  - Function to calculate the Moist static energy
   calculate_adiabat              - Function to calculate an adiabatic parcel ascent with various thermodynamic assumptions
   calculate_dTdp_adiabatic       - Function to calculate the derivative of Temperature with respect to pressure along moist adiabat 
   calculate_dhdp_adiabatic       - Function to calculate the derivative of enthalpy (h) with respect to pressure along a reversible moist adiabat 
   calculate_enthalpy             - Function to calculate the enthalpy
   calculate_entropy              - Function to calculate the entropy
   calculate_frac_ice             - Calculate the fraction of liquid and ice for a saturation adjustment
   calculate_theta_ep             - Calculate pseudo-equivelant potential temperature
   desatdT                        - Function to calculate the derivative of the saturation vapor pressure
   e_sat                          - Function to calculate saturation vapor pressure
   invert_enthalpy                - Invert enthalpy to give temperature given total water content and pressure
   invert_theta_ep                - Invert the pseudo-equivelant potetial temperature
   load_constants                 - Load a table of constants for thermodynamic calculations
   q_sat                          - Function to calculate saturation specific humidity
   r_sat                          - Function to calculate saturation mixing ratio
   saturation_adjustment          - Function to calculate the breakdown of water into vapor, liquid and solid 
   test_script                    - Script to perform some tests to ensure these scripts are working as intended

+ODE	################################################################

    rk_step - Function to perform a single step of a Runge Kutta integration to solve an initial value problem.




