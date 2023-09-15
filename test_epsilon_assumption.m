%
% Test how the assumption of constant epsilon with height affects the phase
% space
%

%% Calculate the phase space

% Set the lower and uppre level
Tb = 295;
Tt = 220;
pb = 95000;

% Set the range of entrainment and precip. efficiency that we test
epsilon = linspace(0,2e-3,201);
PE = linspace(0,1,101);

% Set the values of entrianment and precip. efficiency that we plot
epsilon_plot = [0.05 0.1 0.2 0.4 0.8].*1e-3;
PE_plot = [0.2 0.4 0.6 0.8 1];


clear RH CAPE CAPEv RHv 
% Calculate the CAPE and RH following Romps (2016)
for i = 1:length(epsilon)
    for j = 1:length(PE)

        % Simple calculation where we pick a value for gamma
        [CAPE(i,j),RH(i,j),CAPE_simple] = calculate_CAPE_theory(Tb,Tt,pb,epsilon(i),PE(j),'constant');

        % Consitent calculation where gamma (and epsilon) vary with height
        % The input value of epsilon is that for cloud base
        [CAPEv(i,j),RHv(i,j),CAPE_simple] = calculate_CAPE_theory(Tb,Tt,pb,epsilon(i),PE(j),'gamma');
    end
end

%% Make the phase space figure
fig.bfig(18,18)
axes
set(gca,'xlim',[0.5 1],'ylim',[0 6000])
hold on

lab_locx = [0.55 0.65 0.75 0.85 0.9];
for i = 1:length(epsilon_plot)

    I = find(epsilon == epsilon_plot(i));
    pp = plot(RH(I,:),CAPE(I,:),'color',[0.5 0.5 0.5]);
    fig.inline_label(pp,['\epsilon = ' num2str(epsilon_plot(i).*1000)],lab_locx(i),[],'color',[0.5 0.5 0.5]);

    pp = plot(RHv(I,:),CAPEv(I,:),'--','color',[0.5 0.5 0.5]);
    fig.inline_label(pp,['\epsilon = ' num2str(epsilon_plot(i).*1000)],lab_locx(i),[],'color',[0.5 0.5 0.5]);
   

    
end


lab_locy = [1000 2000 3000 4000 5000];
for j = 1:length(PE_plot)

    I = find(PE == PE_plot(j));
    pp = plot(RH(:,I),CAPE(:,I),'color',[0 0 0]);
    fig.inline_label(pp,['PE = ' num2str(PE_plot(j))],[],lab_locy(j));

    pp = plot(RHv(:,I),CAPEv(:,I),'--','color',[0 0 0]);

end

xlabel('relative humidity')
ylabel('CAPE (J/kg)')

print -dpdf ../Figures/phase_space_constant_a.pdf

%% Now calculate the ZBP model with height for both varying and constant entrainment

z = [500 15000];
epsilon = 0.5e-3;
PE = 0.3;
[T,p,z,RH,Tm,gamma] = calculate_ZBP(z,pb,Tb,epsilon,PE,'const');
[Tv,pv,z,RHv,Tmv,gammav] = calculate_ZBP(z,pb,Tb,epsilon,PE,'gamma');

[Tm,pv,z,RHv,~,gammam] = calculate_ZBP(z,pb,Tb,0,PE,'const');
[Tmv,pv,z,RHv,~,gammamv] = calculate_ZBP(z,pb,Tb,0,PE,'gamma');



fig.bfig(25,15)

subplot(131)
plot(T,z./1000)
hold on
plot(Tv,zv./1000)
plot(Tm,z./1000)
l = legend('\epsilon = const','\epsilon \propto gamma','undilute');
set(l,'box','off')
xlabel('temperature (K)')
ylabel('height (km)')

subplot(132)
plot(1./gamma,z./1000)
hold on
plot(1./gammav,z./1000)
plot(1./gammamv,z./1000)
xlabel('q^* scale height')


subplot(133)
plot(epsilon.*ones(size(z)).*1000,z./1000)
hold on
plot(epsilon.*gammav./gammav(1).*1000,z./1000)
xlabel('entrainment rate (km^{-1})')


print -dpdf ../Figures/profiles_constant_a.pdf





