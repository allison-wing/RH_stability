function [s,varargout] = calculate_entropy(T,p,rt,varargin)
%
% Function to calculate the entropy
%
% s = calculate_entropy(T,p,rt [,rl,ri,type,ice,deltaT])
%
% Calculates entropy based on pressure (p; Pa) Temperature (T; K) and total water mixing ratio (rt; kg/kg))
% Output is entropy per unit mass of dry air 
%
% Optional arguments
%
% rl = liquid mixing ratio (kg/kg)
% ri = ice mixing ratio (kg/kg)
%
% If these are not given, the liquid and solid are divided according
% to the saturation adjustment scheme. See saturation_adjustment.m.
%
% If rt is replaced by 'sat', the saturation entropy is given
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

    e = p.*rv./(c.eps+rv);
    pd = p-e;

    % calculate specific entropy:
    
    s_d = c.cp   .* log(T./c.T0) - c.Rd.* log(pd./c.p00);
    s_v = c.cpv  .* log(T./c.T0) - c.Rv.* log(e./c.e0) + c.Lv0./c.T0;
    s_l = c.cpl  .* log(T./c.T0) - c.Lv0./c.T0         + c.Lv0./c.T0;
    s_i = c.cpi  .* log(T./c.T0) - c.Ls0./c.T0         + c.Lv0./c.T0;
    
   
    % Kludge to deal with zero vapour pressure
    s_v(rv<1e-10) = 0; 
    
     s = s_d + rv.*s_v + rl.*s_l + ri.*s_i;


    if nargout > 1
      % Alternative s to check
      [es,esl,esi] = atm.e_sat(T,thermo_opts{:});
      [Lv,Lf,Ls] = atm.Lv(T,thermo_opts{:});
      s_alt = (c.cp + rt.*c.cpl).*log(T./c.T0) - c.Rd.*log(pd./c.p00) - c.Rv.*rv.*log(e./esl) + c.Rv.*ri.*log(esl./esi) + rv.*Lv./T - ri.*Lf./T;

      varargout{1} = s_alt; 
    end
    
