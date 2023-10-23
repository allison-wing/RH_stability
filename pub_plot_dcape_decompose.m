
clear all

%heights to average RH over
z1 = 2; %km
z2 = 5; %km

%load data
load(['./heightsensitivity/RCEMIP/RH_stab_' num2str(z1) '-' num2str(z2) '.mat'])

%rewrite symbols & sizes
for i = 1:length(small_model_type)
    if strcmp(small_model_type{i},'CRM')
        symbs{i} = '.';
        sizes{i} = 30;
    elseif strcmp(small_model_type{i},'VER')
        symbs{i} = 'o';
        sizes{i} = 10;
    elseif strcmp(small_model_type{i},'LES')
        symbs{i} = 'd';
        sizes{i} = 10;
    end
end

%Scale by dSST or dTLCL?
doLCL = 0; %if 1, LCL. otherwise, SST..

% Load thermodynamic constants
c = atm.load_constants;


%% Calculate decomposition of CAPE changes with warming

%Scale by dSST or dTLCL?
dTLCL = T_LCL305-T_LCL295;
if doLCL
    dT = dTLCL;
else
    dT = 10; %SST difference
end

%Calculate theoretical CAPE
%%295%%
[CAPE_theory295,RH_theory295] = calculate_CAPE_theory(T_LCL295,T_LNB295,p_LCL295,eps_imp295small,PE_imp295small,gammaLCL295,'gamma');
%%300%%
[CAPE_theory300,RH_theory300] = calculate_CAPE_theory(T_LCL300,T_LNB300,p_LCL300,eps_imp300small,PE_imp300small,gammaLCL300,'gamma');
%%305%%
[CAPE_theory305,RH_theory305] = calculate_CAPE_theory(T_LCL305,T_LNB305,p_LCL305,eps_imp305small,PE_imp305small,gammaLCL305,'gamma');

%Calculate change in simulated and theoretical CAPE with warming
dCAPEdT = (100./dT).*(cape305small-cape295small)./cape300small; %as %/K with respect to delta T
dCAPEtheorydT = (100./dT).*(CAPE_theory305-CAPE_theory295)./CAPE_theory300; % %/K

%Calculate CAPE derivatives, using 300K values of the parameters
[CAPE_theory300,dC] = calculate_CAPE_derivatives_with_gamma(T_LCL300,p_LCL300,T_LNB300,eps_imp300small,PE_imp300small,gammaLCL300);
%     dCAPEdTLCL = dC(1,:);
%     dCAPEdpLCL = dC(2,:);
%     dCAPEdTLNB = dC(3,:);
%     dCAPEdeps = dC(4,:);
%     dCAPEdPE = dC(5,:);
%     dCAPEdgamma = dC(6,:);

%Calculate change in parameters with warming
dTLCLdT = (T_LCL305-T_LCL295)./dT;
dpLCLdT = (p_LCL305-p_LCL295)./dT;
dTLNBdT = (T_LNB305-T_LNB295)./dT;
depsdT = (eps_imp305small - eps_imp295small)./dT;
dPEdT = (PE_imp305small - PE_imp295small)./dT;
dgammadT = (gammaLCL305 - gammaLCL295)./dT;

%Calculate decomposition, as %/K
dCAPE_TLCL = (100./cape300small).*dTLCLdT.*dC(1,:);
dCAPE_pLCL = (100./cape300small).*dpLCLdT.*dC(2,:);
dCAPE_TLNB = (100./cape300small).*dTLNBdT.*dC(3,:);
dCAPE_eps = (100./cape300small).*depsdT.*dC(4,:);
dCAPE_PE = (100./cape300small).*dPEdT.*dC(5,:);
dCAPE_gamma = (100./cape300small).*dgammadT.*dC(6,:);

%% Explanations for CAPE change with warming
%residual
dCAPE_resid = dCAPEdT - dCAPE_TLCL - dCAPE_pLCL - dCAPE_TLNB - dCAPE_eps - dCAPE_PE - dCAPE_gamma;

%compute CC scaling
if doLCL
    %compute CC scaling with respect to average TLCL, rather than SST
    %100*1/es des/dT = 100*L/(RvT^2) ==> %/K where T = TLCL
    CCscale = 100*c.Lv0/(c.Rv*nanmean(T_LCL300)^2);
else
    %compute CC scaling with respect to SST
    %100*1/es des/dT = 100*L/(RvT^2) ==> %/K where T = SST
    CCscale = 100*c.Lv0/(c.Rv*300^2);
end

% if doLCL
%     %compute CC scaling with respect to average TLCL, rather than SST
%     qs295 = atm.q_sat(nanmean(T_LCL295),nanmean(p_LCL295));
%     qs300 = atm.q_sat(nanmean(T_LCL300),nanmean(p_LCL300));
%     qs305 = atm.q_sat(nanmean(T_LCL305),nanmean(p_LCL305));
%     CCscale = (100/nanmean(dTLCL))*(qs305-qs295)/qs300;
% else
%     %compute CC scaling as dqsat(SST)
%     qs295 = atm.q_sat(295,1014.8*100);
%     qs300 = atm.q_sat(300,1014.8*100);
%     qs305 = atm.q_sat(305,1014.8*100);
%     %CCscale = (nanmean(cape300small)./qs300).*(qs305-qs295)/10; % J/kg/K over 10K SST range
%     CCscale = (100/10).*(qs305-qs295)./qs300; % %/K over 10K SST range
% end

%% plot

%Force things at 305 to be NaN for DALES and DALES-damping to ignore them for
%changes with warming, because they are aggregated at 305
dCAPEdT(strcmp(small_model_list,'DALES')==1)=NaN;
dCAPEdT(strcmp(small_model_list,'DALES-damping')==1)=NaN;
dCAPE_TLCL(strcmp(small_model_list,'DALES')==1)=NaN;
dCAPE_TLCL(strcmp(small_model_list,'DALES-damping')==1)=NaN;
dCAPE_pLCL(strcmp(small_model_list,'DALES')==1)=NaN;
dCAPE_pLCL(strcmp(small_model_list,'DALES-damping')==1)=NaN;
dCAPE_TLNB(strcmp(small_model_list,'DALES')==1)=NaN;
dCAPE_TLNB(strcmp(small_model_list,'DALES-damping')==1)=NaN;
dCAPE_eps(strcmp(small_model_list,'DALES')==1)=NaN;
dCAPE_eps(strcmp(small_model_list,'DALES-damping')==1)=NaN;
dCAPE_PE(strcmp(small_model_list,'DALES')==1)=NaN;
dCAPE_PE(strcmp(small_model_list,'DALES-damping')==1)=NaN;
dCAPE_gamma(strcmp(small_model_list,'DALES')==1)=NaN;
dCAPE_gamma(strcmp(small_model_list,'DALES-damping')==1)=NaN;
dCAPE_resid(strcmp(small_model_list,'DALES')==1)=NaN;
dCAPE_resid(strcmp(small_model_list,'DALES-damping')==1)=NaN;


% figure('position',[100 100 700 1000]); %for 2x1
figure('position',[100 100 1500 500]); %for 1 x2
%%%
%Panel 1: Parameterized
%%%
subplot(1,2,1)
hold on
%plot box plot
h=boxplot([dCAPEdT(iPAR)' dCAPE_TLCL(iPAR)' dCAPE_pLCL(iPAR)' dCAPE_TLNB(iPAR)' dCAPE_eps(iPAR)' dCAPE_PE(iPAR)' dCAPE_gamma(iPAR)' dCAPE_resid(iPAR)'],'Labels',{'Full','T_LCL','p_LCL','T_LNB','Entrain.','PE','gamma','Resid.'},'Symbol','');
set(h,{'linew'},{2})
%plot individual models
for i = 1:length(iPAR)
    plot([1],dCAPEdT(iPAR(i)),'.','MarkerSize',30,'Color',prof300small(iPAR(i)).color)
    %     plot([2],dCAPEtheorydT(iPAR(i)),'.','MarkerSize',30,'Color',prof300small(iPAR(i)).color)
    plot([2],dCAPE_TLCL(iPAR(i)),'.','MarkerSize',30,'Color',prof300small(iPAR(i)).color)
    plot([3],dCAPE_pLCL(iPAR(i)),'.','MarkerSize',30,'Color',prof300small(iPAR(i)).color)
    plot([4],dCAPE_TLNB(iPAR(i)),'.','MarkerSize',30,'Color',prof300small(iPAR(i)).color)
    plot([5],dCAPE_eps(iPAR(i)),'.','MarkerSize',30,'Color',prof300small(iPAR(i)).color)
    plot([6],dCAPE_PE(iPAR(i)),'.','MarkerSize',30,'Color',prof300small(iPAR(i)).color)
    plot([7],dCAPE_gamma(iPAR(i)),'.','MarkerSize',30,'Color',prof300small(iPAR(i)).color)
    plot([8],dCAPE_resid(iPAR(i)),'.','MarkerSize',30,'Color',prof300small(iPAR(i)).color)
end
%plot mean
plot([1.125],nanmean(dCAPEdT(iPAR)),'k*','MarkerSize',12,'LineWidth',2)
plot([2.125],nanmean(dCAPE_TLCL(iPAR)),'k*','MarkerSize',12,'LineWidth',2)
plot([3.125],nanmean(dCAPE_pLCL(iPAR)),'k*','MarkerSize',12,'LineWidth',2)
plot([4.125],nanmean(dCAPE_TLNB(iPAR)),'k*','MarkerSize',12,'LineWidth',2)
plot([5.125],nanmean(dCAPE_eps(iPAR)),'k*','MarkerSize',12,'LineWidth',2)
plot([6.125],nanmean(dCAPE_PE(iPAR)),'k*','MarkerSize',12,'LineWidth',2)
plot([7.125],nanmean(dCAPE_gamma(iPAR)),'k*','MarkerSize',12,'LineWidth',2)
plot([8.125],nanmean(dCAPE_resid(iPAR)),'k*','MarkerSize',12,'LineWidth',2)

title('(a) Parameterized Convection')
ylim([-5 20])
% xlim([0.5 7.5])
% ylabel('Change in CAPE (J kg^{-1} K^{-1})')
ylabel('Change in CAPE (% K^{-1})')
set(gca,'FontSize',20)
hline(CCscale,'k--')
hline(0,'k')

%%%
%Panel 2: Explicit
%%%
subplot(1,2,2)
hold on
%plot box plot
h=boxplot([dCAPEdT(iEXP)' dCAPE_TLCL(iEXP)' dCAPE_pLCL(iEXP)' dCAPE_TLNB(iEXP)' dCAPE_eps(iEXP)' dCAPE_PE(iEXP)' dCAPE_gamma(iEXP)' dCAPE_resid(iEXP)'],'Labels',{'Full','T_LCL','p_LCL','T_LNB','Entrain.','PE','gamma','Resid.'},'Symbol','');
set(h,{'linew'},{2})
%plot individual models
for i = 1:length(iEXP)
    plot([1],dCAPEdT(iEXP(i)),'.','MarkerSize',30,'Color',prof300small(iEXP(i)).color)
    %     plot([2],dCAPEtheorydT(iEXP(i)),'.','MarkerSize',30,'Color',prof300small(iEXP(i)).color)
    plot([2],dCAPE_TLCL(iEXP(i)),'.','MarkerSize',30,'Color',prof300small(iEXP(i)).color)
    plot([3],dCAPE_pLCL(iEXP(i)),'.','MarkerSize',30,'Color',prof300small(iEXP(i)).color)
    plot([4],dCAPE_TLNB(iEXP(i)),'.','MarkerSize',30,'Color',prof300small(iEXP(i)).color)
    plot([5],dCAPE_eps(iEXP(i)),'.','MarkerSize',30,'Color',prof300small(iEXP(i)).color)
    plot([6],dCAPE_PE(iEXP(i)),'.','MarkerSize',30,'Color',prof300small(iEXP(i)).color)
    plot([7],dCAPE_gamma(iEXP(i)),'.','MarkerSize',30,'Color',prof300small(iEXP(i)).color)
    plot([8],dCAPE_resid(iEXP(i)),'.','MarkerSize',30,'Color',prof300small(iEXP(i)).color)
end
%plot mean
plot([1.125],nanmean(dCAPEdT(iEXP)),'k*','MarkerSize',12,'LineWidth',2)
plot([2.125],nanmean(dCAPE_TLCL(iEXP)),'k*','MarkerSize',12,'LineWidth',2)
plot([3.125],nanmean(dCAPE_pLCL(iEXP)),'k*','MarkerSize',12,'LineWidth',2)
plot([4.125],nanmean(dCAPE_TLNB(iEXP)),'k*','MarkerSize',12,'LineWidth',2)
plot([5.125],nanmean(dCAPE_eps(iEXP)),'k*','MarkerSize',12,'LineWidth',2)
plot([6.125],nanmean(dCAPE_PE(iEXP)),'k*','MarkerSize',12,'LineWidth',2)
plot([7.125],nanmean(dCAPE_gamma(iEXP)),'k*','MarkerSize',12,'LineWidth',2)
plot([8.125],nanmean(dCAPE_resid(iEXP)),'k*','MarkerSize',12,'LineWidth',2)

title('(b) Explicit Convection')
ylim([-5 20])
% xlim([0.5 7.5])
% ylabel('Change in CAPE (J kg^{-1} K^{-1})')
ylabel('Change in CAPE (% K^{-1})')
set(gca,'FontSize',20)
hline(CCscale,'k--')
hline(0,'k')

gcfsavepdf(['Fig09-dCAPE-decompose-' num2str(z1) '-' num2str(z2) '.pdf'])

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot of CRM-VER-LES families
figure('position',[100 100 700 1200]);

%remove index corresponding to DALES-LES because it is NaN but its CRM
    %partner is present - want to ignore the latter. 
    iignore = find(strcmp(small_model_list(iLES),'DALES-LES')==1);
    iLES(iignore) = NaN;
    iLES = iLES(~isnan(iLES));

%%%
%Panel 1: CRM of LES counterparts (iLES-2)
%%%
subplot(3,1,1)
hold on
%plot box plot
h=boxplot([dCAPEdT(iLES-2)' dCAPE_TLCL(iLES-2)' dCAPE_pLCL(iLES-2)' dCAPE_TLNB(iLES-2)' dCAPE_eps(iLES-2)' dCAPE_PE(iLES-2)' dCAPE_gamma(iLES-2)' dCAPE_resid(iLES-2)'],'Labels',{'Full','T_LCL','p_LCL','T_LNB','Entrain.','PE','gamma','Resid.'},'Symbol','');
set(h,{'linew'},{2})
%plot individual models
for i = 1:length(iLES-2)
    plot([1],dCAPEdT(iLES(i)-2),'.','MarkerSize',30,'Color',prof300small(iLES(i)-2).color)
    plot([2],dCAPE_TLCL(iLES(i)-2),'.','MarkerSize',30,'Color',prof300small(iLES(i)-2).color)
    plot([3],dCAPE_pLCL(iLES(i)-2),'.','MarkerSize',30,'Color',prof300small(iLES(i)-2).color)
    plot([4],dCAPE_TLNB(iLES(i)-2),'.','MarkerSize',30,'Color',prof300small(iLES(i)-2).color)
    plot([5],dCAPE_eps(iLES(i)-2),'.','MarkerSize',30,'Color',prof300small(iLES(i)-2).color)
    plot([6],dCAPE_PE(iLES(i)-2),'.','MarkerSize',30,'Color',prof300small(iLES(i)-2).color)
    plot([7],dCAPE_gamma(iLES(i)-2),'.','MarkerSize',30,'Color',prof300small(iLES(i)-2).color)
    plot([8],dCAPE_resid(iLES(i)-2),'.','MarkerSize',30,'Color',prof300small(iLES(i)-2).color)
end
%plot mean
plot([1.125],nanmean(dCAPEdT(iLES-2)),'k*','MarkerSize',12,'LineWidth',2)
plot([2.125],nanmean(dCAPE_TLCL(iLES-2)),'k*','MarkerSize',12,'LineWidth',2)
plot([3.125],nanmean(dCAPE_pLCL(iLES-2)),'k*','MarkerSize',12,'LineWidth',2)
plot([4.125],nanmean(dCAPE_TLNB(iLES-2)),'k*','MarkerSize',12,'LineWidth',2)
plot([5.125],nanmean(dCAPE_eps(iLES-2)),'k*','MarkerSize',12,'LineWidth',2)
plot([6.125],nanmean(dCAPE_PE(iLES-2)),'k*','MarkerSize',12,'LineWidth',2)
plot([7.125],nanmean(dCAPE_gamma(iLES-2)),'k*','MarkerSize',12,'LineWidth',2)
plot([8.125],nanmean(dCAPE_resid(iLES-2)),'k*','MarkerSize',12,'LineWidth',2)

title('(a) CRM')
ylim([-5 20])
% xlim([0.5 7.5])
% ylabel('Change in CAPE (J kg^{-1} K^{-1})')
ylabel('Change in CAPE (% K^{-1})')
set(gca,'FontSize',20)
hline(CCscale,'k--')
hline(0,'k')

%%%
%Panel 2: VER of LES counterpart (iLES-1)
%%%
subplot(3,1,2)
hold on
%plot box plot
h=boxplot([dCAPEdT(iLES-1)' dCAPE_TLCL(iLES-1)' dCAPE_pLCL(iLES-1)' dCAPE_TLNB(iLES-1)' dCAPE_eps(iLES-1)' dCAPE_PE(iLES-1)' dCAPE_gamma(iLES-1)' dCAPE_resid(iLES-1)'],'Labels',{'Full','T_LCL','p_LCL','T_LNB','Entrain.','PE','gamma','Resid.'},'Symbol','');
set(h,{'linew'},{2})
%plot individual models
for i = 1:length(iLES-1)
    plot([1],dCAPEdT(iLES(i)-1),'.','MarkerSize',30,'Color',prof300small(iLES(i)-1).color)
    plot([2],dCAPE_TLCL(iLES(i)-1),'.','MarkerSize',30,'Color',prof300small(iLES(i)-1).color)
    plot([3],dCAPE_pLCL(iLES(i)-1),'.','MarkerSize',30,'Color',prof300small(iLES(i)-1).color)
    plot([4],dCAPE_TLNB(iLES(i)-1),'.','MarkerSize',30,'Color',prof300small(iLES(i)-1).color)
    plot([5],dCAPE_eps(iLES(i)-1),'.','MarkerSize',30,'Color',prof300small(iLES(i)-1).color)
    plot([6],dCAPE_PE(iLES(i)-1),'.','MarkerSize',30,'Color',prof300small(iLES(i)-1).color)
    plot([7],dCAPE_gamma(iLES(i)-1),'.','MarkerSize',30,'Color',prof300small(iLES(i)-1).color)
    plot([8],dCAPE_resid(iLES(i)-1),'.','MarkerSize',30,'Color',prof300small(iLES(i)-1).color)
end
%plot mean
plot([1.125],nanmean(dCAPEdT(iLES-1)),'k*','MarkerSize',12,'LineWidth',2)
plot([2.125],nanmean(dCAPE_TLCL(iLES-1)),'k*','MarkerSize',12,'LineWidth',2)
plot([3.125],nanmean(dCAPE_pLCL(iLES-1)),'k*','MarkerSize',12,'LineWidth',2)
plot([4.125],nanmean(dCAPE_TLNB(iLES-1)),'k*','MarkerSize',12,'LineWidth',2)
plot([5.125],nanmean(dCAPE_eps(iLES-1)),'k*','MarkerSize',12,'LineWidth',2)
plot([6.125],nanmean(dCAPE_PE(iLES-1)),'k*','MarkerSize',12,'LineWidth',2)
plot([7.125],nanmean(dCAPE_gamma(iLES-1)),'k*','MarkerSize',12,'LineWidth',2)
plot([8.125],nanmean(dCAPE_resid(iLES-1)),'k*','MarkerSize',12,'LineWidth',2)

title('(b) VER')
ylim([-5 20])
% xlim([0.5 7.5])
% ylabel('Change in CAPE (J kg^{-1} K^{-1})')
ylabel('Change in CAPE (% K^{-1})')
set(gca,'FontSize',20)
hline(CCscale,'k--')
hline(0,'k')

%%%
%Panel 3: LES
%%%
subplot(3,1,3)
hold on
%plot box plot
h=boxplot([dCAPEdT(iLES)' dCAPE_TLCL(iLES)' dCAPE_pLCL(iLES)' dCAPE_TLNB(iLES)' dCAPE_eps(iLES)' dCAPE_PE(iLES)' dCAPE_gamma(iLES)' dCAPE_resid(iLES)'],'Labels',{'Full','T_LCL','p_LCL','T_LNB','Entrain.','PE','gamma','Resid.'},'Symbol','');
set(h,{'linew'},{2})
%plot individual models
for i = 1:length(iLES)
    plot([1],dCAPEdT(iLES(i)),'.','MarkerSize',30,'Color',prof300small(iLES(i)).color)
    %     plot([2],dCAPEtheorydT(iLES(i)),'.','MarkerSize',30,'Color',prof300small(iLES(i)).color)
    plot([2],dCAPE_TLCL(iLES(i)),'.','MarkerSize',30,'Color',prof300small(iLES(i)).color)
    plot([3],dCAPE_pLCL(iLES(i)),'.','MarkerSize',30,'Color',prof300small(iLES(i)).color)
    plot([4],dCAPE_TLNB(iLES(i)),'.','MarkerSize',30,'Color',prof300small(iLES(i)).color)
    plot([5],dCAPE_eps(iLES(i)),'.','MarkerSize',30,'Color',prof300small(iLES(i)).color)
    plot([6],dCAPE_PE(iLES(i)),'.','MarkerSize',30,'Color',prof300small(iLES(i)).color)
    plot([7],dCAPE_gamma(iLES(i)),'.','MarkerSize',30,'Color',prof300small(iLES(i)).color)
    plot([8],dCAPE_resid(iLES(i)),'.','MarkerSize',30,'Color',prof300small(iLES(i)).color)
end
%plot mean
plot([1.125],nanmean(dCAPEdT(iLES)),'k*','MarkerSize',12,'LineWidth',2)
plot([2.125],nanmean(dCAPE_TLCL(iLES)),'k*','MarkerSize',12,'LineWidth',2)
plot([3.125],nanmean(dCAPE_pLCL(iLES)),'k*','MarkerSize',12,'LineWidth',2)
plot([4.125],nanmean(dCAPE_TLNB(iLES)),'k*','MarkerSize',12,'LineWidth',2)
plot([5.125],nanmean(dCAPE_eps(iLES)),'k*','MarkerSize',12,'LineWidth',2)
plot([6.125],nanmean(dCAPE_PE(iLES)),'k*','MarkerSize',12,'LineWidth',2)
plot([7.125],nanmean(dCAPE_gamma(iLES)),'k*','MarkerSize',12,'LineWidth',2)
plot([8.125],nanmean(dCAPE_resid(iLES)),'k*','MarkerSize',12,'LineWidth',2)

title('(c) LES')
ylim([-5 20])
% xlim([0.5 7.5])
% ylabel('Change in CAPE (J kg^{-1} K^{-1})')
ylabel('Change in CAPE (% K^{-1})')
set(gca,'FontSize',20)
hline(CCscale,'k--')
hline(0,'k')


gcfsavepdf('FigSX-dCAPE-decompose-LES.pdf')
