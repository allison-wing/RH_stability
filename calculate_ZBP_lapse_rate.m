function [Gamma,Gamma_m,RH,gamma] = calculate_ZBP_lapse_rate(T,p,epsilon,PE)
% Function to calculate lapse rate based on the simple zero-buoyancy plume
% assumption.
%
% This function uses the implementation of Romps (2014)
% and assumes epsilon = delta
%
% The function is vectorised, so inputs can be of arbritrary size,
% but all variables must be able to be multiplied together element-wise
%
% Inputs:
%    T = temperature (K)
%    p = pressure (Pa)
%    RH = relative humidity (0-1)
%    epsilon = entrianment rate (m^-1)
%
% Outputs:
%    Gamma = -dT/dz (K/m)
%    Gamma_m = moist adiabatic lapse rate (K/m)
%    RH = relative humidity (0-1)
%    gamma = inverse water vapor scale height (1/m)
%
% Requires the +atm package of thermodynamic functions
%

%% Constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Assumption of non-divergent mass flux
delta = epsilon;

% Set alpha
alpha = 1-PE;

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

% Denominator of the moist adiabatic lapse rate
denom = c.cp + qs.*Lv.^2./(c.Rv.*T.^2);              % J/kg

% Moist adiabatic lapse rate
Gamma_m = c.g .* ( 1 + Lv.*qs./(c.Rd.*T) )./( denom );

% Moist adiabatic frational saturation specific humidity gradient
gamma_m = Lv./(c.Rv.*T.^2) .* Gamma_m - c.g./(c.Rd.*T);

% Nondimensional parameter
Q = Lv.^2.*qs./(c.cp.*c.Rv.*T.^2 + Lv.^2.*qs);


%% Calculate ZBP lapse rate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Following Romps (2014) implementation

%Equations from Appendix B2 of Romps (2014)
b1 = c.Rv.*c.cp.*T.^2./Lv + qs.*Lv;
b2 = c.Rv.*c.cp.*T.^2./Lv .* ( delta - alpha.*epsilon + c.g./(c.Rd.*T) ) + qs.*Lv.*(delta-epsilon) - c.g;
b3 = ( c.Rv.*c.cp.*T./(c.Rd.*Lv) - 1).*c.g.*(delta - alpha.*epsilon);

Gamma = c.Rv.*T.^2./Lv .* ( ( -b2 + sqrt(b2.^2 - 4.*b1.*b3) )./ (2.*b1) +c.g./(c.Rd.*T) );


if 0

    % Equations as written in Wing & Singh (2023): These should be exactly the
    % same as above. Tests indicate this is indeed the case.
    a1 = 1;
    a2 = delta + gamma_m - epsilon.*(alpha + (1-alpha).*Q);
    a3 = -epsilon.*gamma_m.*(1-alpha).*Q;

    Gamma_test = Gamma_m + c.Rv.*T.^2./(2.*Lv) .* ( (a2.^2 - 4.*a1.*a3).^0.5 - a2 );


    disp(['Gamma (R14) = ' Gamma])
    disp(['Gamma (W23) = ' Gamma_test])
    
end


%% Calculate relative humidity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fractional gradient of saturation specific humidity
gamma = Lv.*Gamma./(c.Rv.*T.^2) - c.g./(c.Rd.*T);

% Relative humidity
RH = (delta + alpha.*gamma - alpha.*epsilon) ./ (delta + gamma - alpha.*epsilon);


