function [dTdp,varargout] = calculate_dTdp_adiabatic(T,rt,p,varargin)
%
% Function to calculate the derivative of Temperature with respect to pressure along moist adiabat 
% dTdp = calculate_dTdp(T,rt,p,[gamma,type,ice,deltaT]) 
%
% We assume that the liquid and solid water mixing ratios (rl and ri respectively)
% are divided based on the temperature such that:
%
% rl = f*rc
% ri = (1-f)*rc
% 
% with f = f(T) is a funciton of temperature only, and rc is the condensed 
% water amount.
%
% Then, an adiabatically lifted parcel obeys the following equations 
% for the enthalpy, h as a function of pressure, p:
%
% dh/dp = 1/rho_d
% dh/dp  = cpm*dT/dp + Lv*dr_sat/dp + (1-f)*Lf*dr_sat/dp + Lf*rc*df/dT*dT/dp
%
% where cpm = cpd + rv*cpv + rl*cl+ri*ci, Lv and Lf are the latent heats of
% vaporization and melting, respectively and rho_d  is the dry air density.
%
% Combining the above equations, an expression for dT/dp may be derived.
%
% If gamma = 1, it is assumed that rc = 0 at all times.
%
% If rt is replaced by the string 'sat', we assume that rt = rs, the
% saturation mixing ratio, and rc = 0;




  c = atm.load_constants(varargin{2:end});

 gamma = c.gamma;
  if nargin >= 4; gamma = varargin{1}; end
  
  if ischar(rt)
      rs = atm.r_sat(T,p,varargin{2:end});
      rv = rs;
      rt = rv;
      rl = zeros(size(rv));
      ri = zeros(size(rv));
  else
     [rv,rl,ri,rs] = atm.saturation_adjustment(p,T,rt,varargin{2:end});

     if gamma == 1
        rl = zeros(size(rv));
        ri = zeros(size(rv));
     end
  end
  
  rc = rl+ri;
  
  cpm = c.cp + rv.*c.cpv + rl.*c.cpl + ri.*c.cpi;
  Rm  = c.Rd + rv.*c.Rv;
  
  
  % Unsaturated
  dTdp = Rm.*T./(p.*cpm);
 

  % Need to check the sensitivity to this number 
  sat =  rt./rs > 0.999;
  
  
  Lv = c.Lv0 + (c.cpv-c.cpl).*(T(sat)-c.T0);
  Lf = c.Ls0 - c.Lv0 + (c.cpl-c.cpi).*(T(sat)-c.T0);

   
  
  desdT = atm.desatdT(T(sat),varargin{2:end});
  [~,fice,dfliqdT] = atm.calculate_frac_ice(T(sat),varargin{2:end});

  Lm = Lv + fice.*Lf;
  
  numer =  1 + Lm.*rv(sat)./(c.Rd.*T(sat))  ;
  denom =  1 + Lm./cpm(sat) .* ( (1+rv(sat)./c.eps) ) .*(c.eps+rv(sat)).* desdT./p(sat) + Lf.*dfliqdT.*rc(sat)./cpm(sat);
  
  dTdp(sat) = dTdp(sat).*numer./denom;
  
  if nargout > 1; varargout{1} = rv; end
  if nargout > 2; varargout{2} = sat; end
      


end
