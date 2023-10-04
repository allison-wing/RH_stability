function [CAPE,RH] = calculate_CAPE_theory_const_epsilon(Tb,Tt,pb,epsilon,PE)
%
% Calculate the theoretical CAPE value according to the thoery of Romps (2016)
%
% Inputs are:
%
%   Tb = cloud-base temperature
%   pb = cloud-base pressure
%   Tt = temperature of level of neutral buoyancy
%
%   epsilon = entrainment rate (m^-1)
%   PE = preciptiation efficiency

c = atm.load_constants;

%% Calculate the zero-buoyancy plume solution for constant entrainment

[T,p,z,RH,Tm,gamma] = calculate_ZBP([500 20000],pb,Tb,epsilon,PE,'const');

% Assume dz is constant for now 
dz = z(2)-z(1);

Istrat = false(size(T));

% This is a super annoying loop
for i = 1:size(T,1)
   Istrat(i,:) = squeeze(T(i,:))'<Tt(:);
end
Tm(Istrat) = T(Istrat);

%% Calculate CAPE using a simple sum
CAPE = squeeze(sum( c.g.*(Tm-T)./T .*dz ));

%% Calculate relative humidity
rho = p./(c.Rd.*T);
rho(z<2000 | z>5000) = 0;
RH = squeeze(sum(rho.*RH.*dz)./sum(rho.*dz));



