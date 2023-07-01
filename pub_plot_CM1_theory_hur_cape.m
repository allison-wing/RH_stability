clear all

%Define lower troposphere
z1 = 2*1000; %m
z2 = 5*1000; %m

c = atm.load_constants;

%% Read in data for CM1 simulations with precipitation efficiency output and altered microphysics
CM1_expts = {'RCEMIPN96_dx1000.0_SST300_mean.mat' ...
    'RCEMIPN96_dx1000.0_SST300_vt-1_asnow1.0_arain1.0_thic0.01_thcl0.01_mean.mat' ...
    'RCEMIPN96_dx1000.0_SST300_vt-1_asnow10.0_arain10.0_thic0.01_thcl0.01_mean.mat' ...
    'RCEMIPN96_dx1000.0_SST300_vt-1_asnow3.0_arain3.0_thic0.01_thcl0.01_mean.mat' ...
    'RCEMIPN96_dx1000.0_SST300_vt-1_asnow5.0_arain5.0_thic0.01_thcl0.01_mean.mat' ...
    'RCEMIPN96_dx1000.0_SST300_vt15.0_asnow1.0_arain1.0_thic0.01_thcl0.01_mean.mat' ...
    'RCEMIPN96_dx1000.0_SST300_vt15.0_asnow1.0_arain1.0_thic0.1_thcl0.1_mean.mat' ...
    'RCEMIPN96_dx1000.0_SST300_vt15.0_asnow1.0_arain1.0_thic1.0_thcl1.0_mean.mat' ...
    'RCEMIPN96_dx1000.0_SST300_vt15.0_asnow10.0_arain10.0_thic0.01_thcl0.01_mean.mat' ...
    'RCEMIPN96_dx1000.0_SST300_vt15.0_asnow10.0_arain10.0_thic1.0_thcl1.0_mean.mat' ...
    'RCEMIPN96_dx1000.0_SST300_vt15.0_asnow5.0_arain5.0_thic0.01_thcl0.01_mean.mat' ...
    'RCEMIPN96_dx1000.0_SST300_vt15.0_asnow5.0_arain5.0_thic1.0_thcl1.0_mean.mat' ...
    'RCEMIPN96_dx1000.0_SST300_vt2.0_asnow1.0_arain1.0_thic1.0_thcl1.0_mean.mat' ...
    'RCEMIPN96_dx1000.0_SST300_vt8.0_asnow1.0_arain1.0_thic1.0_thcl1.0_mean.mat' ...
    'RCEMIPN96_dx1000.0_SST300_vt8.0_asnow5.0_arain5.0_thic1.0_thcl1.0_mean.mat'};

for i = 1:length(CM1_expts)
    
    %load data
    load(['/Users/awing/Dropbox/RCEMIP_local/RH_stability/CM1_altered_microphysics/' CM1_expts{i}])
    profCM1(i).ta = T_mean;
    profCM1(i).tv = Trho_mean;
    profCM1(i).z = z;
    profCM1(i).dz = dz;
    profCM1(i).p = p_mean;
    profCM1(i).hur = RH_mean;
    profCM1(i).hus = qv_mean;
    profCM1(i).rv = rv_mean;
    profCM1(i).rho = rho_mean;
    profCM1(i).pe = P_mean/P_gross; %precipitation efficiency
    peCM1(i) = profCM1(i).pe;
    profCM1(i).peproxy =  P_mean./(qc_path_mean+qi_path_mean+qs_path_mean+qr_path_mean+qg_path_mean); %precipitation efficiency proxy
    %     profCM1(i).peproxy =  P_mean./(qc_path_mean+qi_path_mean); %precipitation efficiency proxy
    peproxyCM1(i) = profCM1(i).peproxy;
    
    %Calculate the LCL, based on lowest level T,rv,p
    [T_LCL(i),p_LCL(i)] = atm.calculate_LCL(profCM1(i).ta(1),profCM1(i).rv(1),profCM1(i).p(1));
    
    %Calculate CAPE & LNB, based on lowest level T,rv,p
    [CAPE(i),p_LNB(i),T_rho_LNB(i),T_LNB(i)] = calculate_CAPE(profCM1(i).ta(1),profCM1(i).rv(1),profCM1(i).p(1),profCM1(i).tv,profCM1(i).p,profCM1(i).ta,profCM1(i).z);
    
    % Calculate lower-troposphere indices
    [~,I2] = min(abs(profCM1(i).z-z2));
    [~,I1] = min(abs(profCM1(i).z-z1));
    I = false(length(profCM1(i).z),1); I(I1:I2) = true;
    
    % Calculate lower-tropospheric humidity
    hurCM1_avg(i) = sum(profCM1(i).hur.*profCM1(i).rho.* I.*profCM1(i).dz)./sum(profCM1(i).rho.* I.*profCM1(i).dz);
    
    %calculate theory for CAPE based on LCL and LNB, don't plot
    [PE_R16,epsilon_R16,RH_R16,CAPE_R16] = fn_plot_CAPE_and_RH_RCE(0,T_LCL(i),p_LCL(i),T_LNB(i));
    
    %Diagnose PE & epsilon implied by theory
    [PE_impCM1(i),eps_impCM1(i)] = find_PE_epsilon(PE_R16,epsilon_R16,double(CAPE_R16),double(RH_R16),double(CAPE(i)),double(hurCM1_avg(i)));  
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PLOTS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  CM1 Microphysics Simulations, CAPE vs. RH on top of theory, colored by real PE
[rCM1_mic,pCM1_mic] = corrcoef(CAPE,hurCM1_avg);
pCM1 =polyfit(hurCM1_avg,CAPE,1);

figure;
[PE_R16,epsilon_R16,RH_R16,CAPE_R16] = fn_plot_CAPE_and_RH_RCE(1,nanmean(T_LCL),nanmean(p_LCL),nanmean(T_LNB));
hold on
scatter(hurCM1_avg,CAPE,100,peCM1,'filled','HandleVisibility','off')
plot(hurCM1_avg,pCM1(1)*hurCM1_avg+pCM1(2),'Color',[0.5 0 0.5],'LineWidth',3,'DisplayName',sprintf('r = %3.2f; p = %3.2f',rCM1_mic(1,2),pCM1_mic(1,2)))

set(gca,'FontSize',16)
legend('location','northwest')
ylabel('CAPE (J kg^{-1})')
xlabel('Lower-Tropospheric Relative Humidity')
title('Perturbed Physics Ensemble')
cc = colorbar;
ylabel(cc,'Precip. Efficiency')
xlim([0.5 1])
% ylim([0 300])

gcfsavepdf(['Fig03-theory-hur-cape' num2str(z1/1000) '-' num2str(z2/1000) '-scatter-CM1.pdf'])

%% PE vs. proxy and theory
figure('position',[100 100 400 800]);
subplot(2,1,1)
pp = polyfit(peCM1,PE_impCM1,1);
[rPE,pPE] = corrcoef(peCM1,PE_impCM1);
plot(peCM1,PE_impCM1,'o','MarkerSize',10,'MarkerFaceColor',[0 0 0], 'MarkerEdgeColor',[0 0 0],'HandleVisibility','off')
hold on
plot(0:0.1:1,0:0.1:1,'k-','HandleVisibility','off') %1:1 line
x = 0:0.1:1;
plot(x,pp(1).*x+pp(2),'k--','DisplayName',sprintf('r = %3.2f; p = %3.2f',rPE(1,2),pPE(1,2))) %best fit line

set(gca,'FontSize',16)
xlabel('Actual Precip. Efficiency')
ylabel('Theory-Implied Precip. Efficiency')
title('(a)')
xlim([0 1])
ylim([0 1])
legend('location','northwest')

subplot(2,1,2)
pp = polyfit(peCM1,peproxyCM1,1);
[rPE,pPE] = corrcoef(peCM1,peproxyCM1);
plot(peCM1,peproxyCM1,'o','MarkerSize',10,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0],'HandleVisibility','off')
hold on
x = 0:0.1:1;
plot(x,pp(1).*x+pp(2),'k--','DisplayName',sprintf('r = %3.2f; p = %3.2f',rPE(1,2),pPE(1,2))) %best fit line

set(gca,'FontSize',16)
xlabel('Actual Precip. Efficiency')
ylabel('Proxy Precip. Efficiency')
title('(b)')
xlim([0 1])
ylim([0 4e-3])
legend('location','northwest')

gcfsavepdf(['Fig04-PE-CM1-' num2str(z1/1000) '-' num2str(z2/1000) '.pdf'])

%%  CM1 Microphysics Simulations, CAPE vs. RH on top of theory, colored by implied entrianment
[rCM1_mic,pCM1_mic] = corrcoef(CAPE,hurCM1_avg);
pCM1 =polyfit(hurCM1_avg,CAPE,1);

figure;
[PE_R16,epsilon_R16,RH_R16,CAPE_R16] = fn_plot_CAPE_and_RH_RCE(1,nanmean(T_LCL),nanmean(p_LCL),nanmean(T_LNB));
hold on
scatter(hurCM1_avg,CAPE,100,eps_impCM1*10^3,'filled','HandleVisibility','off')
plot(hurCM1_avg,pCM1(1)*hurCM1_avg+pCM1(2),'Color',[0.5 0 0.5],'LineWidth',3,'DisplayName',sprintf('r = %3.2f; p = %3.2f',rCM1_mic(1,2),pCM1_mic(1,2)))

set(gca,'FontSize',16)
legend('location','northwest')
ylabel('CAPE (J kg^{-1})')
xlabel('Lower-Tropospheric Relative Humidity')
title('Perturbed Physics Ensemble')
cc = colorbar;
ylabel(cc,'Entrainment')
xlim([0.5 1])
% ylim([0 300])
caxis([0.1 0.4])

gcfsavepdf(['FigXX-theory-hur-cape' num2str(z1/1000) '-' num2str(z2/1000) '-scatter-CM1-eps.pdf'])





