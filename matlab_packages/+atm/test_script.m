%
% Script to perform some tests to ensure these scripts are working as intended
%

c = atm.load_constants;

test_sat_vap = 1;
test_adiabat = 1;

test_multi_adiabat = 0;


%% Saturation vapor pressure testing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if test_sat_vap

T = 200:300;

deltaT = 40;
ice  = 1;

[e_sat_cc,esl,esi] = atm.e_sat(T,c.type,1,deltaT);
e_sat_b = atm.e_sat(T,'bolton',1,deltaT);
e_sat_t = atm.e_sat(T,'teten',1,deltaT);

fig.bfig('a5p')

ax1 = axes('position',[0.15 0.62 0.7 0.35])
plot(T,e_sat_cc,'k')
hold on
plot(T,esl,'b--')
plot(T,esi,'r--')


box off
ylabel('vapor pressure (Pa)')
set(gca,'yscale','log')

l = legend('e^*','e^*_l','e^*_i');
set(l,'box','off','location','northwest')

ax1 = axes('position',[0.15 0.2 0.7 0.35])
plot(T,e_sat_b./e_sat_cc,'r')
hold on
plot(T,e_sat_t./e_sat_cc,'g')



l = legend('Bolton','Teten');
set(l,'box','off','location','northwest')

box off

xlabel('temperature (K)')
ylabel('vapor pressure ratio')

set(gca,'xlim',[230 300])


caption = 'Fig. 1 (Top) saturation vapor pressure uisng constant c_p method (black) and saturation vapor pressure over liquid (blue) and ice (red). (Bottom) Ratio of saturation vapor pressure calculated using Bolton (red) and Teten (green) methods compared to the default method (const. c_p).';
aa = annotation('textbox',[0.1 0.003 0.8 0.1],'string',caption,'linestyle','none','horizontalalignment','center','fontsize',8)

end

%% Adiabat calculations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if test_adiabat

Ts = 298.15;
%rs = 0.021428240835937;
rs = 0.02;
ps = 95000;

p = 95000:-2000:10000;
p = p(:);

Ts1000 = 300;
rs1000 = 0.02;
ps1000 = 100000;
p1000 = 100000:-1000:10000;



gamma = 0;
type = 'default';
ice = 0;
deltaT = 40;
[T_noice_r,rv,rl_noice_r,ri_noice_r,T_rho_noice_r] = atm.calculate_adiabat(Ts,rs,ps,p,gamma,type,ice,deltaT);

gamma = 1;
type = 'default';
ice = 0;
deltaT = 40;
[T_noice_p,rv,rl_noice_p,ri_noice_p,T_rho_noice_p] = atm.calculate_adiabat(Ts,rs,ps,p,gamma,type,ice,deltaT);

gamma = 0;
type = 'default';
ice = 1;
deltaT = 40;
[T_r,rv,rl_r,ri_r,T_rho_r] = atm.calculate_adiabat(Ts,rs,ps,p,gamma,type,ice,deltaT);

gamma = 1;
type = 'default';
ice = 1;
deltaT = 40;
[T_p,rv,rl_p,ri_p,T_rho_p] = atm.calculate_adiabat(Ts,rs,ps,p,gamma,type,ice,deltaT);

gamma = 0.5;
type = 'default';
ice = 1;
deltaT = 40;
[T_h,rv,rl_h,ri_h,T_rho_h] = atm.calculate_adiabat(Ts,rs,ps,p,gamma,type,ice,deltaT);



% Calculate adiabat using no-ice entropy to compare

type = 'default';
ice = 0;
deltaT = 40;

ss = atm.calculate_entropy(Ts,ps,rs,type,ice,deltaT);
T_entropy = zeros(size(p));
T_entropy(1) = Ts;
for k = 2:length(p)
   T_entropy(k) = fzero(@(x) ss-atm.calculate_entropy(x,p(k),rs,type,ice,deltaT),T_entropy(k-1));
end

[rv,rl,ri] = atm.saturation_adjustment(p,T_entropy,rs.*ones(size(p)),type,ice,deltaT);
T_rho_entropy = T_entropy.*(1 + rv./c.eps)./(1+rv+rl+ri);

[theta_ep,Tstar] = atm.calculate_theta_ep(Ts,rs,ps);


T_psu_entropy = zeros(size(p));
rv_psu_entropy = zeros(size(p));
T_psu_entropy(1) = Ts;
rv_psu_entropy(1) = rs;

for k = 2:length(p)
   [T_psu_entropy(k),rv_psu_entropy(k)] = atm.invert_theta_ep(theta_ep,rs,Tstar,p(k)); 
end

T_rho_psu_entropy = T_psu_entropy.*(1 + rv_psu_entropy./c.eps)./(1+rv_psu_entropy);


% Calculate some adiabats with a dry-adiabatic layer as well
type = 'default';
ice = 1;
deltaT = 40;

gamma = 0;
[T_r1000,rv1000,rl_r1000,ri_r1000,T_rho_r1000] = atm.calculate_adiabat(Ts1000,rs1000,ps1000,p1000,gamma,type,ice,deltaT);

gamma = 1;
[T_p1000,rv1000,rl_r1000,ri_r1000,T_rho_p1000] = atm.calculate_adiabat(Ts1000,rs1000,ps1000,p1000,gamma,type,ice,deltaT);

gamma = 0.5;
[T_h1000,rv1000,rl_r1000,ri_r1000,T_rho_h1000] = atm.calculate_adiabat(Ts1000,rs1000,ps1000,p1000,gamma,type,ice,deltaT);





%%%%%%%%%%%%%%%%%


fig.bfig('a5l')

ax1 = axes('position',[0.1 0.25 0.8 0.7])

plot(T_r1000,p1000./100,'k')
hold on
plot(T_p1000,p1000./100,'r')
plot(T_h1000,p1000./100,'g')

set(gca,'ydir','reverse')
l = legend('reversible','pseudoadiabatic','\gamma = 0.5')
set(l,'box','off','location','northwest')
xlabel('temperature (K)')
ylabel('pressure (hPa)')
box off

caption = ['Fig. 1: Temperature for adiabatic parcel ascents assuming reversible (black) and pseudoadiabatic (red) thermodynamics with a mixed-phase range between 233.15 and 273.15 K.. Adiabat initialized with (T,r,p) = (' num2str(Ts1000) ' K, ' num2str(rs1000) ' kg/kg, ' num2str(p1000(1)./100) ' hPa)'];
aa = annotation('textbox',[0.1  0.06 0.8 0.1],'string',caption,'linestyle','none','horizontalalignment','center','fontsize',8)

%%%%%%%%%%%%%%%%

fig.bfig('a5l')

ax1 = axes('position',[0.1 0.25 0.35 0.7])
plot(T_noice_r - T_entropy,p./100,'k')
hold on
plot(T_noice_p - T_entropy,p./100,'r')
plot(T_psu_entropy - T_entropy,p./100,'b')

set(gca,'ydir','reverse')
xlabel('temperarture (K)')
ylabel('pressure (hPa)')
box off

l = legend('reversible','pseudoadiabatic','Bolton-pseudo')
set(l,'box','off','location','southwest')

ax2 = axes('position',[0.55 0.25 0.35 0.7])

plot(T_rho_noice_r - T_rho_entropy,p./100,'k')
hold on
plot(T_rho_noice_p- T_rho_entropy,p./100,'r')
plot(T_rho_psu_entropy - T_rho_entropy,p./100,'b')
    
box off    
set(gca,'ydir','reverse')
xlabel('density temperarture (K)')
ylabel('pressure (hPa)')

caption = ['Fig. 3: Temperature (left) and density temperature (right) for adiabatic parcel ascents assuming reversible (black) and pseudoadiabatic (red) thermodynamics and assuming no ice formation. The profiles are plotted as an anomaly from a control ascent that is calculated based on conservation of entropy, assuming no ice formation (an exact solution for the reversible case). Approximate pseudoadiabatic ascent calculated by assuming conservation of pseudo-equivalent potential temperature as defined by Bolton (1980) is shown in blue. Adiabat initialized with (T,r,p) = (' num2str(Ts) ' K, ' num2str(rs) ' kg/kg, ' num2str(p(1)./100) ' hPa)'];
aa = annotation('textbox',[0.1  0.06 0.8 0.1],'string',caption,'linestyle','none','horizontalalignment','center','fontsize',8)



%%%%%%%%%%%%%%%

fig.bfig('a5l')

ax1 = axes('position',[0.1 0.25 0.25 0.7])

plot(T_r - T_entropy,p./100,'k')
hold on
plot(T_p - T_entropy,p./100,'r')
plot(T_h - T_entropy,p./100,'g')

set(gca,'ydir','reverse')
l = legend('reversible','pseudoadiabatic','\gamma = 0.5')
set(l,'box','off','location','northwest')
xlabel('temperature (K)')
ylabel('pressure (hPa)')
box off

ax2 = axes('position',[0.4 0.25 0.25 0.7])

plot(T_rho_r - T_rho_entropy,p./100,'k')
hold on
plot(T_rho_p - T_rho_entropy,p./100,'r')
plot(T_rho_h - T_rho_entropy,p./100,'g')

set(gca,'ydir','reverse')
l = legend('reversible','pseudoadiabatic','\gamma = 0.5')
set(l,'box','off','location','northwest')
xlabel('density temperature (K)')
set(gca,'yticklabel','')
box off

ax3 = axes('position',[0.7 0.25 0.25 0.7])



plot(rl_r,p./100,'k:')
hold on
plot(ri_r,p./100,'k--')
plot(ri_r+rl_r,p./100,'k-')

plot(rl_p,p./100,'r:')
plot(ri_p ,p./100,'r--')
plot(ri_p+rl_p,p./100,'r-')

plot(rl_h,p./100,'g:')
plot(ri_h ,p./100,'g--')
plot(ri_h+rl_h,p./100,'g-')

set(gca,'ydir','reverse')
xlabel('mixing ratio (kg/kg)')
set(gca,'yticklabel','')

l = legend('liquid','ice','total')
set(l,'box','off','location','southeast')
box off

caption = ['Fig. 4: Temperature (left), density temperature (middle), and condensate mixing ratio (right) for parcel ascents including ice and using reversible (blue) pseudo-adiabatic (red) and intermediate thermodynamics. The intermediate case is calculated by assuming half of the condensate created at a given level is precipitated out. The temperature profiles are plotted as ananomaly from a control parcel ascent calculated assuming conservation of entropy with no ice. Adiabat initialized with (T,r,p) = (' num2str(Ts) ' K, ' num2str(rs) ' kg/kg, ' num2str(p(1)./100) ' hPa)'];
aa = annotation('textbox',[0.1 0.05 0.8 0.1],'string',caption,'linestyle','none','horizontalalignment','center','fontsize',8)



end


if test_multi_adiabat



Ts = [280 290 300 310; 280 290 300 310];
rs = [0.002 0.01 0.02 0.03; 0.02 0.02 0.02 0.02];
ps = [98000 97000 94000 100000; 98000 97000 94000 100000];

sigma =1:-0.02:0.05;
sigma = sigma(:);
sigma_mat = repmat(sigma,[1 size(ps)]);
ps_mat = repmat(ps,[1 1 length(sigma)]);
ps_mat = permute(ps_mat,[3 1 2]);
p = sigma_mat.*ps_mat;


gamma = 0;
type = 'default';
ice = 1;
deltaT = 40;
[T,rv,rl,ri,T_rho] = atm.calculate_adiabat(Ts,rs,ps,p,gamma,type,ice,deltaT);


end


%% Test the wetbulb calculations
c = atm.load_constants;

T = 260:305;
p = 100000:-1000:50000;
[T,p] = meshgrid(T,p);

e = 0.6.*atm.e_sat(T);
r = c.eps.*e./(p-e);


Td = atm.calculate_dewpoint(T,p,r,'SAM');



[Tw,Tw_orig] = atm.calculate_isobaric_wetbulb(T,p,r,'SAM');

return

fig.bfig('a5p')
subplot(211)
plot(T,T-Td,'k')
hold on
plot(T,T-Tw,'r')
plot(T,T-Tw_orig,'r--')
l = legend('T-T_d','T-T_w','T-T_w(simple)','location','northwest');
set(l,'box','off')
box off
ylabel('depression (K)')

subplot(212)
plot(T,Tw-Tw_orig,'r--')
title('Difference between simple & comprehensive method')
xlabel('temperature (K)')
ylabel('\delta T (K)')

box off





