function [Lcond,Lfrz,Ldep] = Lv(T,varargin)
% Calculate the latent heats of vaporization and freezing at given
% temperature
%
%

% Get constants
c = atm.load_constants(varargin{:});

Lcond = c.Lv0 + (c.cpv - c.cpl ).*(T-c.T0);
Ldep = c.Ls0 + (c.cpv - c.cpi ).*(T-c.T0);
Lfrz = Ldep-Lcond;


end
