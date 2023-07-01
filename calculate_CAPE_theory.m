    function [CAPE,RH,CAPE_simple] = calculate_CAPE_theory(Tb,Tt,pb,epsilon,PE)
%
% Calculate the theoretical CAPE value according to the thoery of Romps (2016)
%
% Inputs are:
%
%   Tb = cloud-base temperature
%   pb = cloud-base pressure
%   Tt = temperature of level of neutral buoyancy
%
%   epsilon = entrainment rate (m^-1)
%   PE = preciptiation efficiency


%% A couple of parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% For Romps (2016)
T0 = (Tb+Tt)/2;              % temperature scale  (K)

gamma = 1./4000;             % Inverse of water vapor scale height (m^-1)

%% Thermodynamics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Load thermodynamic constants
c = atm.load_constants;


% saturation specific humidity at cloud base
qs = atm.q_sat(Tb,pb);







%% Romps (2016) theory %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% a is a nondimensional parameter (Eq 3 of R16)
a = epsilon.*PE./gamma;

% Relative humidity (Eq 4 of R16)
RH = ( 1-PE+a )./(1+a);


%% Simplified CAPE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Temperature scale
Tc = c.Rd.*T0./( c.Rd.*c.Lv0./(c.Rv.*T0) -c.cp );				

% Simplified CAPE (Eq 17 of R16)
CAPE_simple = a./(1+a) .* c.Lv0.*qs/T0.*(Tb - Tt - Tc);


%% Full theoretical CAPE calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Eq 10 of R16
f = c.Lv0./(c.Rv.*T0.^2) - c.cp./(c.Rd.*T0);

% Eq 9 of R16
y = c.Lv0.*qs./( (1+a).*c.Rd.*T0 ) .* exp( c.Lv0.*qs./( (1+a).*c.Rd.*T0 ) );
y0 = c.Lv0.*qs./( c.Rd.*T0 ) .* exp( c.Lv0.*qs./( c.Rd.*T0 ) );

% Lambert w function
Wy = lambertw(y);
Wy0 = lambertw(y0);
Wey = lambertw(exp(-f.*(Tb-Tt)).*y);
Wey0 = lambertw(exp(-f.*(Tb-Tt)).*y0);

% Eq 12 of R16
CAPE = c.Rd./(2.*f).*( ...
                       Wy  .* ( 2 - 2.*f.*(Tb-Tt) + Wy ) ...
                     - Wey .* ( 2                 + Wey ) ...
                     - Wy0 .* ( 2 - 2.*f.*(Tb-Tt) + Wy0 ) ...
                     + Wey0.* ( 2                 + Wey0 ) ...
                     );


