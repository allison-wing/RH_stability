function [T_LCL,p_LCL,varargout] = calculate_LCL(T,rt,p,varargin)
%
% Function to calculate the lifted condensation level (LCL) for a parcel.
%
%
% [T_LCL,p_LCL,[T_LCLl,p_LCLl,T_LCLi,p_LCLi]] = calculate_LCL(T,r,p,[type,ice,deltaT])
%
%
% Finds the temperature (T_LCL) and pressure (p_LCL) at which a parcel
% with a given temperature (T), mixing ratio (r) and pressure (p) becomes
% saturated.
%
% Phase changes are assumed to have a mixed phase range, with consitent
% alterations to the saturation vapour pressure. These options are set by
% the optional input arguments (see e_sat.m for details).
%
% Optional output arguments give the lifted condensation level (LCL; liquid
% water) and the lifted deposition level (LDL; for ice).
%
% The parcel is assumed to conserve its potential temperature and its
% mixing ratio as it is lifted. These assumptions allow the derivation of
% an analytic expression for the pressure of the LCL in terms of Lambert W
% functions. See Romps (2017) for details.
%
% The analytic expression is strictly only valid for a thermodynamic
% formulation in which the specific heat capacities are assumed constant.
% If a formulation of the saturation vapour pressure that is not consistent
% with this assumption is used, then the LCL calculated here will not
% exactly match an accurate numerical computation, However, the error
% induced by this apprixmation is likely to be very small.
%
%
% The MATLAB library for the Lambert W function appears to take some time
% to load on the first call to the W function within a session. This means
% that calling calculate_LCL will be rather slow on the first call, but
% will speed up for subsequent calls. This should be taken into account
% when testing functions that use calculate_LCL (for example
% calculate_adiabat) for speed.
%
% If the LCL is within the mixed-phase range, the solution for the LCL is
% no longer analytic, but an iteration must be performed. This has the
% potenial to slow down the calculation further, but testing indicates that
% only a few iterations are normally required.
%

% Load constants
c = atm.load_constants(varargin{1:end});

% Calculate water species
[r,~,~] = atm.saturation_adjustment(p,T,rt,varargin{1:end});

% Vapor dependent variables
cpm = c.cp + r.*c.cpv;
Rm = c.Rd + r.*c.Rv;

% Vapor pressures
e = p.*r./(r+c.eps);
[es,esl,esi] = atm.e_sat(T,varargin{1:end});

% Relative humidities
RHl = e./esl;
RHi = e./esi;

% Calculate lifted condensation temperature
ac = cpm./Rm + (c.cpl-c.cpv)./c.Rv;
bc = - ( c.Lv0 - (c.cpv-c.cpl).*c.T0 )./(c.Rv.*T);
cc = bc./ac;

% Calculate LCL (see Romps 2017)
T_LCLl = cc.*lambertw(-1,RHl.^(1./ac).*cc.*exp(cc)).^(-1).*T;


% Calculate lifted deposition temperature
ac = cpm./Rm + (c.cpi-c.cpv)./c.Rv;
bc = - ( c.Ls0 - (c.cpv-c.cpi).*c.T0 )./(c.Rv.*T);
cc = bc./ac;

% Calculate LDL (see Romps 2017)
T_LCLi = cc.*lambertw(-1,RHi.^(1./ac).*cc.*exp(cc)).^(-1).*T;


%% Set the temperature of lifted saturation level

% First set LCL to condensation level
T_LCL = T_LCLl;

if c.ice > 0

  % Below homogenous freezing temperature set T_LCL to deposition level
  T_LCL(T_LCLl < c.T0-c.deltaT)   = T_LCLi(T_LCLl < c.T0-c.deltaT);
   
  % For Temperatre of LCL is in mixed-phase range
  I = T_LCLl > c.T0-c.deltaT & T_LCLl < c.T0;

  % Need to iterate to find correct value
 
  iter = 0;
  dT = ones(size(T_LCL));

  while max(abs(dT(:)))>c.epsT

      [es,~,esi] = atm.e_sat(T_LCL(I));
      RHi_LCL = es./esi;
      T_LCL_new = cc(I).*lambertw(-1,(RHi(I)./RHi_LCL).^(1./ac(I)).*cc(I).*exp(cc(I))).^(-1).*T(I);

      dT = T_LCL_new-T_LCL(I);
      T_LCL(I) = T_LCL_new;
      iter = iter+1;

  end

else
 
  T_LCLi = T_LCLl;

end

p_LCL = p.*(T_LCL./T).^(cpm./Rm);



if nargout > 2; varargout{1} = T_LCLl; end
if nargout > 3; p_LCLl = p.*(T_LCLl./T).^(cpm./Rm); varargout{2} = p_LCLl; end
if nargout > 3; varargout{3} = T_LCLi; end
if nargout > 4; p_LCLi = p.*(T_LCLi./T).^(cpm./Rm); varargout{4} = p_LCLi; end
