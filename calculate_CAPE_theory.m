function [CAPE,RH,CAPE_simple] = calculate_CAPE_theory(Tb,Tt,pb,epsilon,PE,gammaLCL,varargin)
%
% Calculate the theoretical CAPE value according to the thoery of Romps (2016)
%
% Inputs are:
%
%   Tb = cloud-base temperature (K)
%   pb = cloud-base pressure (Pa)
%   Tt = temperature of level of neutral buoyancy (K)
%
%   epsilon = entrainment rate (m^-1)
%   PE = preciptiation efficiency
%   gammaLCL = gamma at the LCL (m^-1)
%   epsilon_type = 'constant' for constant epsilon with height or 'gamma' for epsilon/gamma constant with height (optional)


%% Optional inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Select either constant entrainment or assume epsilon/gamma = constant
epsilon_type = 'constant';
if nargin > 6; epsilon_type = varargin{1}; end


%% Thermodynamic parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% For Romps (2016)
T0 = (Tb+Tt)/2;              % temperature scale  (K)

% Load thermodynamic constants
c = atm.load_constants;

% saturation specific humidity at cloud base
qs = atm.q_sat(Tb,pb);


%% Calculate "a" parameter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(epsilon_type,'constant')

    % This is simple, but slightly inconsistent with the original Romsp
    % (2016) theory
    gammab = 1./4000;             % Inverse of water vapor scale height (m^-1)

else

    % Following Romps (2016), we take a = epsilon*PE/gamma as constant with
    % height. Also taking PE as constant, this implies epsilon varies with
    % height following gamma. We assume input epsilon is the entrainment 
    % rate at cloud base. We take as input the simulated gamma at cloud base to
    % give "a".
    
    gammab = gammaLCL;

end

%disp(['water vapour scale height at cloud base is ', num2str(1./gammab) ' m'])

% a is a nondimensional parameter (Eq 3 of R16)
a = epsilon.*PE./gammab;



%% Relative humidity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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


