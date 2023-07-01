function [r,rl,ri,varargout] = saturation_adjustment(p,T,r_t,varargin)
%
% Function to calculate the breakdown of water into vapor, liquid and solid 
% based on different microphysical assumptions
%
%  [r,rl,ri,[rs]] = saturation_adjustment(p,T,r_t [type,ice,deltaT])
%
% Outputs: 
% r = water vapor mixing ratio (kg/kg)
% rl = liquid water mixing ratio (kg/kg)
% ri = solid water mixing ratio (kg/kg)
% rs = saturation mixing ratio (kg/kg)
%
% Inputs:
% p = pressure (Pa)
% T = temperature (K)
% r_t = total water mixing ratio (kg/kg)
% deltaT = mixed phase range (K) (default = 40)


c = atm.load_constants(varargin{1:end});

    es = atm.e_sat(T,varargin{1:end});

    rs = c.eps.*(es./(p-es));

    r = min(r_t,rs);

    r(r<0) = 0;
   
    [fliq,fice] = atm.calculate_frac_ice(T,varargin{1:end});
    

    rl = fliq.*(r_t-r);
    ri = fice.*(r_t-r);

    rl(abs(r_t-r)<1e-8) = 0;
    ri(abs(r_t-r)<1e-8) = 0;
    
    if nargout > 3; varargout{1} = rs; end
