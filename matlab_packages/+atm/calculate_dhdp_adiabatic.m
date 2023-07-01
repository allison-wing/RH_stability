function dhdp = calculate_dhdp(h,rt,p,Tguess,varargin)
%
% Function to calculate the derivative of enthalpy (h) with respect to pressure along a reversible moist adiabat 
% given the enthalpy, total water mixing ratio (rt) and the pressure (p) and an iniital guess for T (Tguess).
% dhdp = calculate_dhdp(h,rt,p,Tguess,[type,ice,deltaT]) 


  c = atm.load_constants(varargin{1:end});

  [T,rv,~,~] = atm.invert_enthalpy(h,rt,p,Tguess,varargin{1:end});
  Rm = c.Rd + rv.*c.Rv;
  rho = p./(Rm.*T);
  dhdp = 1./rho;
  


end
