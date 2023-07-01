function s = calculate_MSE(T,p,z,rt,varargin)
%
% Function to calculate the Moist static energy
%
% s = calculate_MSE(T,p,z,rt [rl,ri,type,ice,deltaT])
% s = moist static energy       (J/kg)
%
% T = temperature               (K)
% rt = total water mixing ratio (kg/kg)
% p = pressure                  (Pa)
% z = height                    (m)
%
% Optional arguments
%
% rl = liquid mixing ratio (kg/kg)
% ri = ice mixing ratio (kg/kg)
%
% If these are not given, the liquid and solid are divided according 
% to the saturation adjustment scheme. See saturation_adjustment.m.
%
% If rt is replaced by 'sat', the saturation MSE is given
%

%% Read arguments
if nargin == 4 				% No optional inputs
   thermo_opts = {};
   condensed_input = 0;

elseif ischar(varargin{1}) 		% no condensed water input
   thermo_opts = {varargin{1:end}};
   condensed_input = 0;

else					% Both liquid and solid water given
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

% Convert to specific humidities
qv = rv./(1+rt);
ql = rl./(1+rt);
qi = ri./(1+rt);
qt = rt./(1+rt);

% calculate moist static energy
sd = c.cp.*(T-c.T0)  + c.g.*z;
sv = c.cpv.*(T-c.T0) + c.g.*z + c.Lv0;
sl = c.cpl.*(T-c.T0) + c.g.*z;
si = c.cpi.*(T-c.T0) + c.g.*z - (c.Ls0-c.Lv0);

s = sd.*(1-qt) + qv.*sv + ql.*sl + qi.*si;

