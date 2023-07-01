function [theta_ep,Tstar] = calculate_theta_ep(T,r,p)

% Calculate pseudo-equivelant potential temperature
% assuming no ice
% Use Bolton formula, taken from Emanuel (1994), p 132
% This assumes no ice!
%
% [theta_ep,Tstar] = calculate_theta_ep(T,r,p)
%
% T = temperature (K)
% r = mixing ratio (kg/kg)
% p = pressure (Pa)

c = atm.load_constants('bolton');

evap = r.*p./(c.eps + r); 

Tstar = 2840./(3.5.*log(T) - log(evap./100) - 4.805 ) + 55; 

theta_ep = T.*(100000./p).^(0.2854.*(1-0.28.*r)).*exp( r .* (1+0.81.*r) .* (3376./Tstar -2.54) );
