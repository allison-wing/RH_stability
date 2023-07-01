function [desdT,desldT,desidT] = desatdT(T,varargin)
%
% Function to calculate the derivative of the saturation vapor pressure
% with respect to temperature.
%
% Assumes a mixed-phase range of deltaT K.
%
% [dessp[,deslsp,desisp]] = desatsp(T[,type,ice,deltaT])
%
% desdp = derivative of saturation vapor pressure w.r.t temperature (Pa/K)
% desldp = derivative for saturation over liquid (Pa/K)
% desidp = derivative for saturation over ice (Pa/K)
%
% T = temperature (K)
%
% type = {'default','bolton','teten','sam'}
% ice =  {[0] , [1] }
% deltaT = mixed-phase temperature ranges
%
%
% This code uses the Clausius-Clapeyron relation. If the default thermodynamics
% is used, this is consistent with the definition of saturation vapor pressure.
% For other thermo types, the saturation vapor pressure may not satisfy the 
% Clausius-Clapeyron equation exactly. But this is a numerical issue inherent
% to using approximate formulas for the saturation vapor curve.
%

c = atm.load_constants(varargin{1:end});



[~,esl,esi] = atm.e_sat(T,varargin{1:end});
[fliq,fice,dfliqdT] = atm.calculate_frac_ice(T,varargin{1:end});

Lv = c.Lv0 + (c.cpv-c.cpl).*(T-c.T0);
Ls = c.Ls0 + (c.cpl-c.cpi).*(T-c.T0);

desldT = Lv./(c.Rv.*T.^2).*esl;
desidT = Ls./(c.Rv.*T.^2).*esi;


% Remember to add the component associated with changes in ice fraction
desdT = desldT.*fliq + desidT.*fice + dfliqdT.*(esl-esi);




end
