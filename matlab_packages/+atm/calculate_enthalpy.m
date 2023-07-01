function h = calculate_enthalpy(T,p,rt,varargin)
%
% Function to calculate the enthalpy
%
% s = calculate_enthalpy(T,p,rt,{rl,ri,type,ice,deltaT})
%
% Calculates entropy based on pressure (p) Temperature (T) and 

% either:
%         vapor mixing ratio (rt) and given values of liquid and and solid water (rl,ri) 
% or:
%         using total water content (rt) with mixed-phase to find
%         the amount of liquid and solid.
%



%% Read arguments
if nargin == 3                          % No optional inputs
   thermo_opts = {};
   condensed_input = 0;

elseif ischar(varargin{1})              % no condensed water input
   thermo_opts = {varargin{1:end}};
   condensed_input = 0;

else                                    % Both liquid and solid water given
   thermo_opts = {varargin{3:end}};
   condensed_input = 1;
end

c = atm.load_constants(thermo_opts{:});


if ischar(rt)
    rv = atm.r_sat(T,p,thermo_opts{:});
    rl = zeros(size(rv));
    ri = zeros(size(rv));
    rt = rv;

elseif condensed_input
    rl = varargin{1};
    ri = varargin{2};
    rv = rt - rl-ri;

else
  % calculate proportions of vapor, liquid and solid
  [rv,rl,ri] = atm.saturation_adjustment(p,T,rt,thermo_opts{:});

end



    % calculate enthalpy per unit dry air mass:
    
    % Calculate Enthalpy
    hd = c.cp.*(T - c.T0);
    hv = c.cpv.*(T - c.T0) + c.Lv0;
    hl = c.cpl.*(T - c.T0);
    hi = c.cpi.*(T - c.T0) + c.Lv0 - c.Ls0;
    
    %h = (1-qi+ql+qv).*hd + qv.*hv + ql.*hl + qi.*hi;
    h = hd + rv.*hv + rl.*hl + ri.*hi;

    
    
    
