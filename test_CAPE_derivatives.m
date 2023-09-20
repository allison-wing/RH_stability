% This script tests the derivatives produced by the function
% "calculate_CAPE_derivatives" by comparing them to a numerical
% approximation.
%

% Type
%ent_type = 'constant';
ent_type = 'gamma';

% Reference values
Tb0 = 295;
pb0 = 95000;
Tt0 = 210;
epsilon0 = 0.1e-3;
PE0 = 0.5;

% Testing range
Tb = 273:0.01:310;
Tt = 180:0.01:250;
epsilon = (0.1:0.01:2).*1e-3;
PE = 0.1:0.01:0.99;

%% Figure 

fig.bfig(15,15)


%% Test the derivative with respect to temperature

[CAPE,dC] = calculate_CAPE_derivatives(Tb,pb0.*ones(size(Tb)),Tt0.*ones(size(Tb)),epsilon0.*ones(size(Tb)),PE0.*ones(size(Tb)),ent_type);

% Numerically estimate dCAPE/dTb
dCAPE_numerical = gradient(CAPE,Tb);

subplot(2,2,1)
plot(Tb,dC(1,:)./CAPE.*100)
hold on
plot(Tb,dCAPE_numerical./CAPE.*100,'--')

xlabel('cloud base temp (K)')
ylabel('% CAPE increase with Tb (%/K)')
l = legend('analytic','numerical');
set(l,'box','off')
box off


%% Test derivative with respect to convection depth

[CAPE,dC] = calculate_CAPE_derivatives(Tb0.*ones(size(Tt)),pb0.*ones(size(Tt)),Tt,epsilon0.*ones(size(Tt)),PE0.*ones(size(Tt)),ent_type);

% Numerically estimate dCAPE/dTb
dCAPE_numerical = gradient(CAPE,Tt);

subplot(2,2,2)
plot(Tt,dC(3,:)./CAPE.*100)
hold on
plot(Tt,dCAPE_numerical./CAPE.*100,'--')

xlabel('LNB temperature (K)')
ylabel('% CAPE increase with T_t (%/K)')
box off


%% Test derivative with respect to entrainment

[CAPE,dC] = calculate_CAPE_derivatives(Tb0.*ones(size(epsilon)),pb0.*ones(size(epsilon)),Tt0.*ones(size(epsilon)),epsilon,PE0.*ones(size(epsilon)),ent_type);

% Numerically estimate dCAPE/dTb
dCAPE_numerical = gradient(CAPE,epsilon);

subplot(2,2,3)
plot(epsilon.*1000,dC(4,:)./CAPE.*100./1000)
hold on
plot(epsilon.*1000,dCAPE_numerical./CAPE.*100./1000,'--')

xlabel('entrainment rate (km^{-1})')
ylabel('% CAPE increase with \epsilon (% km)')
box off


%% Test derivative with respect to entrainment

[CAPE,dC] = calculate_CAPE_derivatives(Tb0.*ones(size(PE)),pb0.*ones(size(PE)),Tt0.*ones(size(PE)),epsilon0.*ones(size(PE)),PE,ent_type);

% Numerically estimate dCAPE/dTb
dCAPE_numerical = gradient(CAPE,PE);

subplot(2,2,4)
plot(PE.*100,dC(5,:)./CAPE)
hold on
plot(PE.*100,dCAPE_numerical./CAPE,'--')

xlabel('precipitation efficiency (%)')
ylabel('% CAPE increase with PE (%/%)')
box off




