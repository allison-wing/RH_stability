function [CAPE,dC] = calculate_CAPE_derivatives_with_gamma(Tb,pb,Tt,epsilon,PE,gamma)
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
%   gamma    - fractional vertical gradient of saturation specific humidity
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
%          dC(6,:) = dCAPE/dgamma       [J/kg m]
%       
%   where the derivatives are taken while holding other variables constant.




%% Constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           

% Load thermodynamic constants
type = 'default';
c = atm.load_constants(type);


%% Temperature difference %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate its value
dT = Tb - Tt;

% Set the derivative of dT with respect to the inputs
ddT = zeros([6,size(Tb)]);
ddT(1,:) = 1;
ddT(3,:) = -1;


%% temperature scale %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tropospheric temperature (K)
% T0 = T0(Tb,Tt)

% Calculate its value
T0 = (Tb + Tt)./2; 

% Calculate the derivatives
dT0 = zeros([6,size(Tb)]);

dT0(1,:) = 0.5;
dT0(3,:) = 0.5;



%% saturation specific humidity at cloud base %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% qs = qs(Tb,pb)

% Calculate its value
qs = atm.q_sat(Tb,pb,type);

% calculate the derivatives
dqs = zeros([6,size(Tb)]);

% With respect to temperature
dqs(1,:) = ( 1 + (1/c.eps-1).*qs ).*atm.Lv(Tb)./(c.Rv.*Tb.^2).*qs;

% With respect to pressure
dqs(2,:) = - ( 1 + (1/c.eps-1).*qs ).*qs./pb;


%% Non-dimensional parameter a
% a = a(epsilon,PE)
% a is a nondimensional parameter (Eq 3 of R16)


% Calculate the value of a
a = epsilon.*PE./gamma;

% Derivatives of a
da = zeros([6,size(Tb)]);
da(4,:) = da(4,:) + PE./gamma;
da(5,:) = da(5,:) + epsilon./gamma;
da(6,:) = da(6,:) - epsilon.*PE./gamma.^2;


%% Parameter f 
% f = f(T0)
% Eq 10 of R16

% Calculate the value
f = c.Lv0./(c.Rv.*T0.^2) - c.cp./(c.Rd.*T0);

% Derivatives of f
df_dT0 = -2.*c.Lv0./(c.Rv.*T0.^3) + c.cp./(c.Rd.*T0.^2);

df = zeros([5,size(Tb)]);
for i = 1:6 
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

dz = zeros([6,size(Tb)]);
dz0  = zeros([6,size(Tb)]);
for i = 1:6
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

dy = zeros([6,size(Tb)]);
dy0 = zeros([6,size(Tb)]);
for i = 1:6
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

dye = zeros([6,size(Tb)]);
dye0 = zeros([6,size(Tb)]);
for i = 1:6
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
for i = 1:6
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



dC = zeros([6,size(Tb)]);
for i = 1:6
    dC(i,:)   = dC_dWy.*dWy(i,:) + dC_dWy0.*dWy0(i,:) + dC_dWye.*dWye(i,:) + dC_dWye0.*dWye0(i,:) + ...
                dC_df.*df(i,:) + dC_ddT.*ddT(i,:);
end



