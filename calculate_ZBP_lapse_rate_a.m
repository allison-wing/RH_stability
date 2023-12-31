function [Gamma,Gamma_m,RH,gamma] = calculate_ZBP_lapse_rate_a(T,p,a,PE)
% Function to calculate lapse rate based on the simple zero-buoyancy plume
% assumption.
%
% This function uses the implementation of Romps (2016)
% and assumes epsilon = delta
%
% The function is vectorised, so inputs can be of arbritrary size,
% but all variables must be able to be multiplied together element-wise
%
% Inputs:
%    T = temperature (K)
%    p = pressure (Pa)
%    RH = relative humidity (0-1)
%    a = epsilon*PE/gamma = non-dimensional parameter of Romps (2016)
%
% Outputs:
%    Gamma = -dT/dz (K/m)
%    Gamma_m = moist adiabatic lapse rate (K/m)
%    gamma = fractional gradient of specific humidity
%    RH = relative humidity (0-1)
%
% Requires the +atm package of thermodynamic functions
%


%% Thermodynamics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Thermodynamic constants
c = atm.load_constants;

% Calculate saturation thermodynamics
es = atm.e_sat(T);
qs = c.eps.*es./(p-es);

% Latent heat constant

Lv = c.Lv0.*ones(size(T));
Lv(T<c.T0) = c.Ls0;


%% Calculate moist adiabatic lapse rate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Lapse rate from Romps (2016) equation B19
Gamma = ( (1+a).*c.g + qs.*Lv.*c.g./(c.Rd.*T) )./ ( (1+a).*c.cp + qs.*Lv.^2./(c.Rv.*T.^2) );

% Denominator of the moist adiabatic lapse rate
denom = c.cp + qs.*Lv.^2./(c.Rv.*T.^2);              % J/kg

% Moist adiabatic lapse rate
Gamma_m = c.g .* ( 1 + Lv.*qs./(c.Rd.*T) )./( denom );

% Moist adiabatic frational saturation specific humidity gradient
gamma_m = Lv./(c.Rv.*T.^2) .* Gamma_m - c.g./(c.Rd.*T);



%% Calculate relative humidity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fractional gradient of saturation specific humidity
gamma = Lv.*Gamma./(c.Rv.*T.^2) - c.g./(c.Rd.*T);

% Relative humidity
RH = (1-PE+a)./(1+a);


