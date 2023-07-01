function rs = r_sat(T,p,varargin)
%
% Function to calculate saturation mixing ratio
%  
% Calling format:
%
% qs = q_sat(T,p,[type,ice,deltaT)
%
% qs = saturation mixing ratio (kg/kg)
% T = temperature (K)
% p = pressure (Pa)
%
% See e_sat(T) function for info about optional arguments

c = atm.load_constants(varargin{1:end});

es = atm.e_sat(T,varargin{1:end});

rs = c.eps .* es./(p-es);

