function [CAPE,p_LNB,T_rho_LNB,T_LNB] = calculate_CAPE(Tinit,rinit,pinit,T_rho_env,p_env,T_env,z,varargin)
%
% Calculate the CAPE based on input parcel properties and an environmental
% density temperature profile
%
% [CAPE,p_LNB,T_rho_LNB,T_LNB] = calculate_CAPE(Tinit,rinit,pinit,T_rho_env,p_env,T_env,z,varargin)
%T_env added as input, T_LNB added as output by AAW
%
%

%% Check arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make sure the inputs are vectors
if ~isvector(p_env) || ~isvector(T_rho_env)
    error('Input environmental properties must be Nx1 vector')
else
    p_env = p_env(:);
    T_rho_env = T_rho_env(:);
    T_env = T_env(:); %AAW
end

% Defaults
gamma = 0; %  0 = reversible, 1 = pseudoadiabatic %Marty had as 1, AAW set to 0
deltaT = 40;
ice = 1;
use_Tv = 1;
z_top = 0; % lift to LNB. if z_top is something other than 0, then only lift to that point.

if nargin >= 8; gamma = varargin{1}; end
if nargin >= 9; deltaT = varargin{2}; end
if nargin >= 10; ice = varargin{3}; end
if nargin >= 11; use_Tv = varargin{4}; end
if nargin >= 12; z_top = varargin{5}; end
%AAW added 1 to all nargin because added another required input



%% Load constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c = atm.load_constants;


%% Calculate parcel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make a pressure vector that only includes levels above pinit
T_rho_env = T_rho_env(p_env<=pinit);
z = z(p_env<=pinit);
p_env = p_env(p_env<=pinit);
T_env = T_env(p_env<=pinit); %AAW



if gamma == 2
    
    if use_Tv == 1; error('cannot use virtual temperature with ZBP model'); end
    if ~ismember(pinit,p_env); error('Need to have parcel pressure included in pressure vector for ZBP model'); end
    
    
    
    % Use a zero-entrainment plume as the "adiabat"
    [T_par,~,z,RH,Tm] = calculate_ZBP(z,p_env,Tinit,0,1);
    
    % If using ZBP model, don't use virtual temp.
    T_rho_par = T_par;
    
else
    % Calculate parcel ascent based on input parcel properties
    [T_par,r_par,rl_par,ri_par,T_rho_par] = atm.calculate_adiabat(Tinit,rinit,pinit,p_env,gamma,'default',ice,deltaT);
end

if ~use_Tv; T_rho_par = T_par; end



%% Estimate level of neutral buoyancy %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if z_top == 0
    
    % Find last positive buoyancy point
    I = find(T_rho_par>T_rho_env,1,'last');
    
    % Interpolate to find the LNB
    p_LNB = interp1(T_rho_par(I:I+1)-T_rho_env(I:I+1),p_env(I:I+1),0);
    T_rho_LNB = interp1(p_env(I:I+1),T_rho_env(I:I+1),p_LNB);
    T_rho_env_LNB = T_rho_LNB;
    T_LNB = interp1(p_env(I:I+1),T_env(I:I+1),p_LNB); %AAW added

else
    % Find last height below z-top
    I = find(z<=z_top,1,'last');
    
    % Interpolate to z_top
    p_LNB = interp1(z(I:I+1),p_env(I:I+1),z_top);
    T_rho_LNB = interp1(z(I:I+1),T_rho_par(I:I+1),z_top);
    T_rho_env_LNB = interp1(z(I:I+1),T_rho_env(I:I+1),z_top);
    T_LNB = interp1(z(I:I+1),T_env(I:I+1),z_top);

    
    % remove everything above z_top
    T_rho_env = T_rho_env(1:I);
    T_rho_par = T_rho_par(1:I);
    p_env = p_env(1:I);
end

%% Include the LNB in the vectors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Sort pressure vector
[p_env,p_order] = sort([p_env; p_LNB],'descend');

% Sort density temperature
T_rho_par = [T_rho_par; T_rho_LNB];
T_rho_par = T_rho_par(p_order);

T_rho_env = [T_rho_env; T_rho_env_LNB];
T_rho_env = T_rho_env(p_order);


%% Calculate CAPE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Integrate positive buoyancy with trapezoidal method.
CAPE = -c.Rd.*trapz(log(p_env),(T_rho_par - T_rho_env).*((T_rho_par-T_rho_env)>0));


