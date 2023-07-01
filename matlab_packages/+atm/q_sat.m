function qs = q_sat(T,p,varargin)
%
% Function to calculate saturation specific humidity
%  
% qs = q_sat(T,p [,qt,type,ice,deltaT])
%
% qs = saturation specific humidity (kg/kg)
% T = temperature (K)
% p = pressure (Pa)
% qt = mass fraction of total water 
%
% Other optional arguments are described in e_sat

qt_arg   = 0;
type_arg = 1;

if nargin >=3 & isnumeric(varargin{1}); qt_arg = 1; type_arg = 2; end

c = atm.load_constants(varargin{type_arg:end});

    % Calculate q_sat assuming no liquid/solid water
    es = atm.e_sat(T,varargin{type_arg:end});
    qs = c.eps.*(es./(p-es.*(1-c.eps)));

    % If include qt as an argument, adjust for the existence of
    % liquid/solid water
    if qt_arg
       qt = varargin{1};
       qs(qt>qs) = (1-qt(qt>qs)).*c.eps.*(es(qt>qs)./(p(qt>qs)-es(qt>qs)));
    end

end
