function [T,rv,rl,ri] = invert_enthalpy(h,p,rt,varargin)
%
% Invert enthalpy to give temperature given total water content and pressure
% Uses an iteration procedure. Should converge in a few iterations
%
% [T,rv,rl,ri] = invert_enthalpy(h,p,rt,[T,type,ice,deltaT])
%
% h = enthalpy (J/kg/K; as given by calculate_enthalpy.m)
% p = pressure (Pa)
% rt = total water mixing ratio
%      if rt = 'sat', we invert the enthalpy assuming the atmosphere is at saturation
% T = initial guess of temperature (otherwise assumed to be triple point)

if nargin >= 4 && isnumeric(varargin{1})
   T = varargin{1};   
   opts = {varargin{2:end}}; 
else 
   opts = {varargin{1:end}};
   T = 273.16;
end

c = atm.load_constants(opts{:});

T_err_max = 1;
i_iter = 0;


% Calculate the enthalpy given our guess of the temperature
h_guess = atm.calculate_enthalpy(T,p,rt,opts{1:end});

while T_err_max > 0.0001
 i_iter = i_iter+1;


 % Improve the guess using a linear approximation
 % ISSUE: the tangent approximation has a disconinuity at saturation
 dhdT = c.cp.*ones(size(T));
 if ischar(rt)
   rv = atm.r_sat(T,p,opts{1:end});
   sat = true(size(T));
 else
    [rv,rl,ri] = atm.saturation_adjustment(p,T,rt,opts{1:end});
    sat = rt>rv;
 end

 dhdT(sat) = c.cp + (c.Lv0.^2./(c.Rv.*T(sat).^2)).*rv(sat);
 
 
 % Calculate a new temperature guess
 T_new = T + (h - h_guess)./dhdT;

 h_guess = atm.calculate_enthalpy(T_new,p,rt,opts{1:end});
  
 
 % Get the error for this guess (convert to temperature)
 T_err= abs((h-h_guess)./c.cp);
 [T_err_max,Ierr] = max(T_err(:));
 T = T_new;

 if i_iter == 10
    % If we get to 10 iterations, we are probably in a loop between two values of T far from the correct one.
    % Try to salvage a solution with a last ditch perturbation to the temperature field.
    T(T_err>0.1) = T(T_err>0.1)+0.1;
 elseif i_iter > 20 & i_iter <= 30
    % Now we are really likely to be not converged. Print some output to warn the user
    disp(['iter: ' num2str(i_iter) ', error = ' num2str(T_err_max) ' K'])
    disp(['iter: ' num2str(i_iter) ', (I,h,T,r,p) = (' num2str(Ierr) ',' num2str(h(Ierr)) ',' num2str(T(Ierr)) ',' num2str(rv(Ierr)) ',' num2str(p(Ierr)) ')'])
 elseif i_iter > 30
     % Pull the plug
    error('enthalpy inversion did not converge. This may be because the input data is corrupt.')
 end
 
end

if isstr(rt)
   rv = atm.r_sat(T,p,opts{1:end});
   rl = 0;
   ri = 0;
else
   [rv,rl,ri] = atm.saturation_adjustment(p,T,rt,opts{1:end});
end
