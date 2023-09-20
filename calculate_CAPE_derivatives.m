function [CAPE,dC] = calculate_CAPE_derivatives(Tb,pb,Tt,epsilon,PE,varargin)
%
% Calculate the CAPE according to Romps (2016) analytic formula
% and calculate its derivatives with respect to its input arguments
% 
% Input paramaters:
%
%   Tb       - cloud base temperature (K) (LCL, typically)
%   pb       - cloud base pressure (pa) (LCL, typically)
%   Tt       - temperature of level of neutral buoyancy (K)
%   epsilon  - entrainment rate (m^{-1})
%   PE       - precipitation efficiency
%
% Outputs:
% 
%   CAPE - the CAPE according to Romps (2016)
%   dC   - of size [5, size(Tb)], with the five values being:
%
%          dC(1,:) = dCAPE/dTb          [J/kg/K]
%          dC(2,:) = dCAPE/dpb          [J/kg/Pa]
%          dC(3,:) = dCAPE/dTt          [J/kg/K]
%          dC(4,:) = dCAPE/depsilon     [J/kg m]
%          dC(5,:) = dCAPE/dPE          [J/kg]
%       
%   where the derivatives are taken while holding other variables constant.

%% Optional inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Select either constant entrainment or assume epsilon/gamma = constant
epsilon_type = 'constant';
if nargin > 5; epsilon_type = varargin{1}; end

%% Constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           

% Load thermodynamic constants
type = 'default';
c = atm.load_constants(type);


%% Temperature difference %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate its value
dT = Tb - Tt;

% Set the derivative of dT with respect to the inputs
ddT = zeros([5,size(Tb)]);
ddT(1,:) = 1;
ddT(3,:) = -1;


%% temperature scale %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tropospheric temperature (K)
% T0 = T0(Tb,Tt)

% Calculate its value
T0 = (Tb + Tt)./2; 

% Calculate the derivatives
dT0 = zeros([5,size(Tb)]);

dT0(1,:) = 0.5;
dT0(3,:) = 0.5;



%% saturation specific humidity at cloud base %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% qs = qs(Tb,pb)

% Calculate its value
qs = atm.q_sat(Tb,pb,type);

% calculate the derivatives
dqs = zeros([5,size(Tb)]);

% With respect to temperature
dqs(1,:) = ( 1 + (1/c.eps-1).*qs ).*atm.Lv(Tb)./(c.Rv.*Tb.^2).*qs;

% With respect to pressure
dqs(2,:) = - ( 1 + (1/c.eps-1).*qs ).*qs./pb;


%% Non-dimensional parameter a
% a = a(epsilon,PE)
% a is a nondimensional parameter (Eq 3 of R16)

dgamma = zeros([5,size(Tb)]);
if strcmp(epsilon_type,'constant')

    % This is simple, but slightly inconsistent with the original Romsp
    % (2016) theory
    gamma = 1./4000;             % Inverse of water vapor scale height (m^-1)

else

    % Following Romps (2016), we take a = epsilon*PE/gamma as constant with
    % height. Also taking PE as constant, this implies epsilon varies with
    % height following gamma. We assume input epsilon is the entrainment 
    % rate at cloud base. We then need to calculate gamma at cloud base to
    % give "a".


    % Calculate gamma based on cloud-base temperature 
    % (see Appendix B of Romps, 2014).

    % Non-dimensional constants
    A = c.Lv0./(c.Rd.*Tb);
    B = c.cp.*c.Rv.*Tb.^2/c.Lv0.^2;

    % Derivatives
    dA = zeros([5,size(Tb)]);
    dB = zeros([5,size(Tb)]);

    dA(1,:) = -c.Lv0./(c.Rd.*Tb.^2) ;
    dB(1,:) = 2.*c.cp.*c.Rv.*Tb/c.Lv0.^2 ;

    % Quadratic co-efficiencts
    a1 = c.Lv0 .* ( B + qs );
    a2 = B .* ( c.Lv0.*epsilon.*PE + c.g.*A ) - c.g;
    a3 = c.g.* ( B.*A - 1 ).*epsilon.*PE;

    % Derivatives
    da2 = zeros([5,size(Tb)]);
    da3 = zeros([5,size(Tb)]);

    da1 = c.Lv0 .* ( dB + dqs );
    for i = 1:5
       da2(i,:) = dB(i,:) .* ( c.Lv0.*epsilon.*PE + c.g.*A ) + B.*c.g.*dA(i,:);
       da3(i,:) = c.g.* ( dB(i,:).*A ).*epsilon.*PE + c.g.* ( B.*dA(i,:) ).*epsilon.*PE;
    end
    da2(4,:) = da2(4,:) + B.*c.Lv0.*PE;
    da2(5,:) = da2(5,:) + B.*c.Lv0.*epsilon;

    
    da3(4,:) = da3(4,:) + c.g.* ( B.*A - 1 ).*PE;
    da3(5,:) = da3(5,:) + c.g.* ( B.*A - 1 ).*epsilon;

    % Quadratic formula
    gamma = ( -a2 + ( a2.^2 - 4.*a1.*a3).^0.5 )./(2.*a1);


    % Derivatives
    dgamma_da1 = -0.5.*(a2.^2-4.*a1.*a3).^-0.5.*4.*a3./(2.*a1) - ( -a2 + ( a2.^2 - 4.*a1.*a3).^0.5 )./(2.*a1.^2);
    dgamma_da2 = ( -1 +(a2.^2-4.*a1.*a3).^-0.5.*a2 )./(2.*a1);
    dgamma_da3 = -( (a2.^2-4.*a1.*a3).^-0.5.*a1 )./(a1);

    for i = 1:5 
       dgamma(i,:) = dgamma_da1.*da1(i,:) + dgamma_da2.*da2(i,:) + dgamma_da3.*da3(i,:);
    end
end


% Calculate the value
a = epsilon.*PE./gamma;

% Derivatives of a
da = zeros([5,size(Tb)]);
for i = 1:5 
   da(i,:) = - epsilon.*PE./gamma.^2 .* dgamma(i,:);
end
da(4,:) = da(4,:) + PE./gamma;
da(5,:) = da(5,:) + epsilon./gamma;


%% Parameter f 
% f = f(T0)
% Eq 10 of R16

% Calculate the value
f = c.Lv0./(c.Rv.*T0.^2) - c.cp./(c.Rd.*T0);

% Derivatives of f
df_dT0 = -2.*c.Lv0./(c.Rv.*T0.^3) + c.cp./(c.Rd.*T0.^2);

df = zeros([5,size(Tb)]);
for i = 1:5 
    df(i,:) = df_dT0.*dT0(i,:); 
end


%% z values
% z = z(qs,T0,a)
% z0 = z0(qs,T0)

% Calculate the values of z
z = c.Lv0.*qs./( (1+a).*c.Rd.*T0 );
z0 = c.Lv0.*qs./( c.Rd.*T0 );

% Derivatives
dz_dqs = c.Lv0./( (1+a).*c.Rd.*T0 );
dz_da = - c.Lv0.*qs./( (1+a).^2.*c.Rd.*T0 );
dz_dT0 = - c.Lv0.*qs./( (1+a).*c.Rd.*T0.^2 );

dz0_dqs = c.Lv0./( c.Rd.*T0 );
dz0_dT0 = - c.Lv0.*qs./( c.Rd.*T0.^2 );

dz = zeros([5,size(Tb)]);
dz0  = zeros([5,size(Tb)]);
for i = 1:5
    dz(i,:) = dz_dqs.*dqs(i,:) + dz_dT0.*dT0(i,:) + dz_da.*da(i,:);
    dz0(i,:) = dz0_dqs.*dqs(i,:) + dz0_dT0.*dT0(i,:);
end


%% Arguments to the Lambert W functions

% y = y(z)

% Calculate the values of y
y  = z .* exp( z  );
y0 = z0 .* exp( z0 );

% Derivatives
dy_dz = exp(z) + z.*exp(z);
dy0_dz0 = exp(z0) + z0.*exp(z0);

dy = zeros([5,size(Tb)]);
dy0 = zeros([5,size(Tb)]);
for i = 1:5
    dy(i,:) = dy_dz.*dz(i,:);
    dy0(i,:) = dy0_dz0.*dz0(i,:);   
end


% ye = ye(y,f,dT)

% Calculate the values of ye
ye = exp(-f.*(dT)).*y;
ye0 = exp(-f.*(dT)).*y0;

% Derivatives
dye_dy = exp(-f.*(dT));
dye_df = -dT.*ye;
dye_ddT = - f.*ye;

dye0_dy0 = exp(-f.*(dT));
dye0_df = -dT.*ye0;
dye0_ddT = -f.*ye0;

dye = zeros([5,size(Tb)]);
dye0 = zeros([5,size(Tb)]);
for i = 1:4
    dye(i,:) = dye_dy.*dy(i,:) + dye_df.*df(i,:) + dye_ddT.*ddT(i,:);
    dye0(i,:) = dye0_dy0.*dy0(i,:)+ dye0_df.*df(i,:) + dye0_ddT.*ddT(i,:);   
end


%% Lambert W functions

% Calculate the values
Wy = lambertw(y);
Wy0 = lambertw(y0);
Wye = lambertw(ye);
Wye0 = lambertw(ye0);

% Calculate the derivatives
dWy_dy = Wy./(y.*(1+Wy));
dWy0_dy0 = Wy0./(y0.*(1+Wy0));
dWye_dye = Wye./(ye.*(1+Wye));
dWye0_dye0 = Wye0./(ye0.*(1+Wye0));

dWy = zeros([5,size(Tb)]);
dWy0 = zeros([5,size(Tb)]);
dWye = zeros([5,size(Tb)]);
dWye0 = zeros([5,size(Tb)]);
for i = 1:5
    dWy(i,:)   = dWy_dy.*dy(i,:);
    dWy0(i,:)  = dWy0_dy0.*dy0(i,:);
    dWye(i,:)  = dWye_dye.*dye(i,:);
    dWye0(i,:) = dWye0_dye0.*dye0(i,:);
end

%% Finally the CAPE

% Calculate the CAPE
% Eq 12 of R16
CAPE = c.Rd./(2.*f).*( ...
                       Wy  .* ( 2 - 2.*f.*(dT) + Wy ) ...
                     - Wye .* ( 2              + Wye ) ...
                     - Wy0 .* ( 2 - 2.*f.*(dT) + Wy0 ) ...
                     + Wye0.* ( 2              + Wye0 ) ...
                     );


% Now calculate the derivatives
dC_dWy   =  c.Rd./(2.*f).*( 2 - 2.*f.*(dT) + 2.*Wy );
dC_dWye  = -c.Rd./(2.*f).*( 2              + 2.*Wye );
dC_dWy0  = -c.Rd./(2.*f).*( 2 - 2.*f.*(dT) + 2.*Wy0 );
dC_dWye0 =  c.Rd./(2.*f).*( 2              + 2.*Wye0 );
dC_ddT  =   c.Rd./(2.*f).*(-Wy.*2.*f + Wy0.*2.*f );

dC_df   = -c.Rd./(2.*f.^2) .* ( ...
                       Wy  .* ( 2  + Wy ) ...
                     - Wye .* ( 2  + Wye ) ...
                     - Wy0 .* ( 2  + Wy0 ) ...
                     + Wye0.* ( 2  + Wye0 ) ...
                     );



dC = zeros([5,size(Tb)]);
for i = 1:5
    dC(i,:)   = dC_dWy.*dWy(i,:) + dC_dWy0.*dWy0(i,:) + dC_dWye.*dWye(i,:) + dC_dWye0.*dWye0(i,:) + ...
                dC_df.*df(i,:) + dC_ddT.*ddT(i,:);
end



