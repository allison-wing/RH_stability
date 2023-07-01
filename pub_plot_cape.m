clear all

small_model_list={'CAM5-GCM','CAM6-GCM','CM1','CM1-VER','CM1-LES','CNRM-CM6','DALES',...
    'DALES-VER','DALES-LES','DALES-damping','DALES-damping-VER','DAM','GEOS-GCM',...
    'ICON-LEM-CRM','ICON-LEM-VER','ICON-LEM-LES','ICON-NWP-CRM','MESONH',...
    'MESONH-VER','MESONH-LES','MicroHH','MicroHH-VER','MicroHH-LES','MPAS','SAM-CRM',...
    'SAM-CRM-VER','SAM-CRM-LES','SCALE','UCLA-CRM','UKMO-GA7.1','UKMO-CASIM',...
    'UKMO-RA1-T-hrad','UKMO-RA1-T','UKMO-RA1-T-nocloud','WRF-COL-CRM','WRF-CRM',...
    'WRF-GCM-cps0','WRF-GCM-cps1','WRF-GCM-cps2','WRF-GCM-cps3','WRF-GCM-cps4',...
    'WRF-GCM-cps6'};

small_model_type={'GCM','GCM','CRM','VER','LES','GCM','CRM',...
    'VER','LES','CRM','VER','CRM','GCM',...
    'CRM','VER','LES','CRM','CRM',...
    'VER','LES','CRM','VER','LES','GCRM','CRM',...
    'VER','LES','CRM','CRM','GCM','CRM',...
    'CRM','CRM','CRM','CRM','CRM',...
    'GCM','GCM','GCM','GCM','GCM',...
    'GCM'};

iPAR = find(strcmp(small_model_type,'GCM')==1);
iCRM = find(strcmp(small_model_type,'CRM')==1);
iVER = find(strcmp(small_model_type,'VER')==1);
iLES = find(strcmp(small_model_type,'LES')==1);
iEXP = [iCRM iVER iLES];

symbs = {'.','.','o','o','o','.','d',...
    'd','d','.','.','.','.',...
    'hexagram','hexagram','hexagram','.','^',...
    '^','^','v','v','v','.','s',...
    's','s','.','.','.','.',...
    '.','.','.','.','.',...
    '.','.','.','.','.',...
    '.'}; %for CRM-VER-LES familes

sizes = {25,25,200,100,50,25,200,...
    100,50,25,25,25,25,...
    200,100,50,25,200,...
    100,50,200,100,50,25,200,...
    100,50,25,25,25,25,...
    25,25,25,25,25,...
    25,25,25,25,25,...
    25}; %for CRM-VER-LES familes


colors_model_type=[1 0 0;1 0 0;0 0 1;0.5 0.5 1;0.75 0.75 1;1 0 0;0 0 1;...
    0.5 0.5 1;0.75 0.75 1;0 0 1;0.5 0.5 1;0 0 1;1 0 0;...
    0 0 1;0.5 0.5 1;0.75 0.75 1;0 0 1;0 0 1;...
    0.5 0.5 1;0.75 0.75 1;0 0 1;0.5 0.5 1;0.75 0.75 1;0 0 1;0 0 1;...
    0.5 0.5 1;0.75 0.75 1;0 0 1;0 0 1;1 0 0;0 0 1;...
    0 0 1;0 0 1;0 0 1;0 0 1;0 0 1;...
    1 0.5 0.5;1 0.5 0.5;1 0.5 0.5;1 0.5 0.5;1 0.5 0.5;...
    1 0.5 0.5;];

filebase = '_cfv0-profiles.nc';

SSTs = [295 300 305];


%% Read in small profiles
for i = 1:length(small_model_list)
    for j = 1:length(SSTs)
        filename = ['/Users/awing/Dropbox/RCEMIP_local/var_files/small_v6/' small_model_list{i} '_RCE_small' num2str(SSTs(j)) filebase]
        ncid = netcdf.open(filename);
        
        if SSTs(j)==295
            prof295small(i).model = small_model_list{i};
            prof295small(i).z = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'zg_avg')); %km
            prof295small(i).p = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'pa_avg')); %hPa
            prof295small(i).ta = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'ta_avg')); % K
            prof295small(i).hur = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'hur_avg'));% percent out of 100
            prof295small(i).hus = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'hus_avg')); %g/kg
            prof295small(i).Cr = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'Cr'));
            prof295small(i).Cg = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'Cg'));
            prof295small(i).Cb = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'Cb'));
            prof295small(i).color = [prof295small(i).Cr prof295small(i).Cg prof295small(i).Cb];
            colors_model(i,:) = prof295small(i).color;
            netcdf.close(ncid)
        elseif SSTs(j)==300
            prof300small(i).model = small_model_list{i};
            prof300small(i).z = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'zg_avg')); %km
            prof300small(i).p = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'pa_avg')); %hPa
            prof300small(i).ta = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'ta_avg')); % K
            prof300small(i).hur = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'hur_avg'));% percent out of 100
            prof300small(i).hus = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'hus_avg')); %g/kg
            prof300small(i).Cr = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'Cr'));
            prof300small(i).Cg = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'Cg'));
            prof300small(i).Cb = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'Cb'));
            prof300small(i).color = [prof300small(i).Cr prof300small(i).Cg prof300small(i).Cb];
            colors_model(i,:) = prof300small(i).color;
            netcdf.close(ncid)
        elseif SSTs(j) == 305
            prof305small(i).model = small_model_list{i};
            prof305small(i).z = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'zg_avg')); %km
            prof305small(i).p = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'pa_avg')); %hPa
            prof305small(i).ta = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'ta_avg')); % K
            prof305small(i).hur = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'hur_avg'));% percent out of 100
            prof305small(i).hus = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'hus_avg')); %g/kg
            prof305small(i).Cr = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'Cr'));
            prof305small(i).Cg = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'Cg'));
            prof305small(i).Cb = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'Cb'));
            prof305small(i).color = [prof305small(i).Cr prof305small(i).Cg prof305small(i).Cb];
            colors_model(i,:) = prof305small(i).color;
            netcdf.close(ncid)
        end
    end %end loop over SSTs for reading in data
    
    %Fix GEOS level ordering (flip from top to bottom to bottom to top to match everything else)
    if strcmp(small_model_list(i),'GEOS-GCM')==1
        prof295small(i).z = flipud(prof295small(i).z);
        prof300small(i).z = flipud(prof300small(i).z);
        prof305small(i).z = flipud(prof305small(i).z);
        prof295small(i).p = flipud(prof295small(i).p);
        prof300small(i).p = flipud(prof300small(i).p);
        prof305small(i).p = flipud(prof305small(i).p);
        prof295small(i).ta = flipud(prof295small(i).ta);
        prof300small(i).ta = flipud(prof300small(i).ta);
        prof305small(i).ta = flipud(prof305small(i).ta);
        prof295small(i).hur = flipud(prof295small(i).hur);
        prof300small(i).hur = flipud(prof300small(i).hur);
        prof305small(i).hur = flipud(prof305small(i).hur);
        prof295small(i).hus = flipud(prof295small(i).hus);
        prof300small(i).hus = flipud(prof300small(i).hus);
        prof305small(i).hus = flipud(prof305small(i).hus);
    end
    %Fix UK GA7 units for hus from g/g to g/kg
    if strcmp(small_model_list(i),'UKMO-GA7.1')==1
        prof295small(i).hus = prof295small(i).hus*1000;
        prof300small(i).hus = prof300small(i).hus*1000;
        prof305small(i).hus = prof305small(i).hus*1000;
    end
    
    %Compute moist adiabat & CAPE (wants rv in g/g, p in Pa), skipping models who are missing data
    %Ignore UKMA-RA1-T because it is aggregated
    if strcmp(small_model_list(i),'DALES-VER')==1 || strcmp(small_model_list(i),'DALES-LES')==1 || ...
            strcmp(small_model_list(i),'DALES-damping-VER')==1 || ...
            strcmp(small_model_list(i) ,'UKMO-RA1-T')==1
        %if missing data, set to nan
        prof295small(i).tma = nan(length(prof295small(i).ta),1);
        prof300small(i).tma = nan(length(prof300small(i).ta),1);
        prof305small(i).tma = nan(length(prof305small(i).ta),1);
        prof295small(i).tmarho = nan(length(prof295small(i).ta),1);
        prof300small(i).tmarho = nan(length(prof300small(i).ta),1);
        prof305small(i).tmarho = nan(length(prof305small(i).ta),1);
        cape295small(i)=NaN; cape300small(i)=NaN; cape305small(i)=NaN;
        Tlnb295(i) = NaN; Tlnb300(i) = NaN; Tlnb305(i) = NaN;
        T_LCL295(i) = NaN; T_LCL300(i) = NaN; T_LCL305(i) = NaN;
        p_LCL295(i) = NaN; p_LCL300(i) = NaN; p_LCL305(i) = NaN;
    else
        %convert from specific humidity to mixing ratio
        prof295small(i).rv = prof295small(i).hus/1000./(1-prof295small(i).hus/1000);
        prof300small(i).rv = prof300small(i).hus/1000./(1-prof300small(i).hus/1000);
        prof305small(i).rv = prof305small(i).hus/1000./(1-prof305small(i).hus/1000);
        
        %compute virtual temperature
        prof295small(i).tv = prof295small(i).ta.*(1+prof295small(i).rv/.622)./(1+prof295small(i).rv); %virtual temperature
        prof300small(i).tv = prof300small(i).ta.*(1+prof300small(i).rv/.622)./(1+prof300small(i).rv); %virtual temperature
        prof305small(i).tv = prof305small(i).ta.*(1+prof305small(i).rv/.622)./(1+prof305small(i).rv); %virtual temperature
        
        %compute moist adiabat
        
        %%%%%%%%%%%%%%%%%
        %295%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%
        %Lift from lowest model level
        Tinit = prof295small(i).ta(1);
        rvinit = prof295small(i).rv(1);
        pinit = prof295small(i).p(1);
        
        %Calcluate CAPE 
        [cape295small(i),p_LNB295(i),Tv_LNB295(i),T_LNB295(i)] = calculate_CAPE(Tinit,rvinit,pinit*100,prof295small(i).tv,prof295small(i).p*100,prof295small(i).ta,prof295small(i).z*1000);
        
       
        %%%%%%%%%%%%%%%%%
        %300%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%
        %Lift from lowest model level
        Tinit = prof300small(i).ta(1);
        rvinit = prof300small(i).rv(1);
        pinit = prof300small(i).p(1);
        
        %Calcluate CAPE 
        [cape300small(i),p_LNB300(i),Tv_LNB300(i),T_LNB300(i)] = calculate_CAPE(Tinit,rvinit,pinit*100,prof300small(i).tv,prof300small(i).p*100,prof300small(i).ta,prof300small(i).z*1000);
        
        %%%%%%%%%%%%%%%%%
        %305%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%
        %Lift from lowest model level
        Tinit = prof305small(i).ta(1);
        rvinit = prof305small(i).rv(1);
        pinit = prof305small(i).p(1);
        
        %Calcluate CAPE
        [cape305small(i),p_LNB305(i),Tv_LNB305(i),T_LNB305(i)] = calculate_CAPE(Tinit,rvinit,pinit*100,prof305small(i).tv,prof305small(i).p*100,prof305small(i).ta,prof305small(i).z*1000);
        
    end
    
end

%Force cape at 305 to be NaN for DALES and DALES-damping to ignore them for
%changes with warming, because they are aggregated at 305
cape305small(strcmp(small_model_list,'DALES')==1)=NaN;
cape305small(strcmp(small_model_list,'DALES-damping')==1)=NaN;

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PLOT %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('position',[100 100 1200 400])
%%%%%%%%%%%%%%%%%
%PANEL 1 Histogram of CAPE
%parameterized models
[Nparam,Edgeparam] = histcounts(cape300small(iPAR),'BinWidth',400,'BinLimits',[0 4000],'Normalization','pdf');
for i = 1:length(Nparam)
    xparam(i) = mean(Edgeparam(i:i+1));
end

%explicit models
[Nexpl,Edgeexpl] = histcounts(cape300small(iEXP),'BinWidth',400,'BinLimits',[0 4000],'Normalization','pdf');
for i = 1:length(Nexpl)
    xexpl(i) = mean(Edgeexpl(i:i+1));
end

%CRMs
[NCRM,EdgeCRM] = histcounts(cape300small(iCRM),'BinWidth',400,'BinLimits',[0 4000],'Normalization','pdf');
for i = 1:length(NCRM)
    xCRM(i) = mean(EdgeCRM(i:i+1));
end

%VER
[NVER,EdgeVER] = histcounts(cape300small(iVER),'BinWidth',400,'BinLimits',[0 4000],'Normalization','pdf');
for i = 1:length(NVER)
    xVER(i) = mean(EdgeVER(i:i+1));
end

%LES
[NLES,EdgeLES] = histcounts(cape300small(iLES),'BinWidth',400,'BinLimits',[0 4000],'Normalization','pdf');
for i = 1:length(NLES)
    xLES(i) = mean(EdgeLES(i:i+1));
end

%plot
subplot(1,2,1)
plot(xparam,Nparam,'Color',[1 0 0],'LineWidth',3,'DisplayName','Parameterized')
hold on
plot(xCRM,NCRM,'Color',[0 0 1],'LineWidth',3,'DisplayName','Explicit-CRM')
plot(xVER,NVER,'Color',[0 180/255 1],'LineWidth',3,'DisplayName','Explicit-VER')
plot(xLES,NLES,'Color',[0 1 1],'LineWidth',3,'DisplayName','Explicit-LES')

plot(nanmean(cape300small(iPAR)),0,'o','MarkerSize',10,'MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[1 0 0],'DisplayName',sprintf('PAR: %4.0f J kg^{-1}',nanmean(cape300small(iPAR))))
plot(nanmean(cape300small(iCRM)),0,'o','MarkerSize',10,'MarkerEdgeColor',[0 0 1],'MarkerFaceColor',[0 0 1],'DisplayName',sprintf('CRM: %4.0f J kg^{-1}',nanmean(cape300small(iCRM))))
plot(nanmean(cape300small(iVER)),0,'o','MarkerSize',10,'MarkerEdgeColor',[0 180/255 1],'MarkerFaceColor',[0 180/255 1],'DisplayName',sprintf('VER: %4.0f J kg^{-1}',nanmean(cape300small(iVER))))
plot(nanmean(cape300small(iLES)),0,'o','MarkerSize',10,'MarkerEdgeColor',[0 1 1],'MarkerFaceColor',[0 1 1],'DisplayName',sprintf('LES: %4.0f J kg^{-1}',nanmean(cape300small(iLES))))

set(gca,'FontSize',20)
xlabel('J kg^{-1}')
ylabel('Probability Density')
title('(a) CAPE at 300 K')
legend('location','northwest')

%%%%%%%%%%%%%%
%PANEL 2 Histogram of CAPE change with warming
clear xparam Nparam xCRM NCRM xVER NVER xLES NLES
%compute CC scaling as dqsat(SST)
qs295 = atm.q_sat(295,1014.8*100);
qs300 = atm.q_sat(300,1014.8*100);
qs305 = atm.q_sat(305,1014.8*100);
CCscale = (100/10)*(qs305-qs295)/qs300;

%parameterized models
dCAPE_param = (100/10)*(cape305small(iPAR)-cape295small(iPAR))./cape300small(iPAR); % percent/K over 10K range
[Nparam,Edgeparam] = histcounts(dCAPE_param,'BinWidth',2,'BinLimits',[0 20],'Normalization','pdf');
for i = 1:length(Nparam)
    xparam(i) = mean(Edgeparam(i:i+1));
end

%explicit models
dCAPE_expl = (100/10)*(cape305small(iEXP)-cape295small(iEXP))./cape300small(iEXP); % percent/K over 10K range
[Nexpl,Edgeexpl] = histcounts(dCAPE_expl,'BinWidth',2,'BinLimits',[0 20],'Normalization','pdf');
for i = 1:length(Nexpl)
    xexpl(i) = mean(Edgeexpl(i:i+1));
end

%CRMs
dCAPE_CRM = (100/10)*(cape305small(iCRM)-cape295small(iCRM))./cape300small(iCRM); % percent/K over 10K range
[NCRM,EdgeCRM] = histcounts(dCAPE_CRM,'BinWidth',2,'BinLimits',[0 20],'Normalization','pdf');
for i = 1:length(NCRM)
    xCRM(i) = mean(EdgeCRM(i:i+1));
end

%VER
dCAPE_VER = (100/10)*(cape305small(iVER)-cape295small(iVER))./cape300small(iVER); % percent/K over 10K range
[NVER,EdgeVER] = histcounts(dCAPE_VER,'BinWidth',2,'BinLimits',[0 20],'Normalization','pdf');
for i = 1:length(NVER)
    xVER(i) = mean(EdgeVER(i:i+1));
end

%LES
dCAPE_LES = (100/10)*(cape305small(iLES)-cape295small(iLES))./cape300small(iLES); % percent/K over 10K range
[NLES,EdgeLES] = histcounts(dCAPE_LES,'BinWidth',2,'BinLimits',[0 20],'Normalization','pdf');
for i = 1:length(NLES)
    xLES(i) = mean(EdgeLES(i:i+1));
end

%plot
subplot(1,2,2)
plot(xparam,Nparam,'Color',[1 0 0],'LineWidth',3,'DisplayName','Parameterized')
hold on
plot(xCRM,NCRM,'Color',[0 0 1],'LineWidth',3,'DisplayName','Explicit-CRM')
plot(xVER,NVER,'Color',[0 180/255 1],'LineWidth',3,'DisplayName','Explicit-VER')
plot(xLES,NLES,'Color',[0 1 1],'LineWidth',3,'DisplayName','Explicit-LES')

plot(nanmean(dCAPE_param),0,'o','MarkerSize',10,'MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[1 0 0],'DisplayName',sprintf('PAR: %2.1f %%K^{-1}',nanmean(dCAPE_param)))
plot(nanmean(dCAPE_CRM),0,'o','MarkerSize',10,'MarkerEdgeColor',[0 0 1],'MarkerFaceColor',[0 0 1],'DisplayName',sprintf('CRM: %2.1f %%K^{-1}',nanmean(dCAPE_CRM)))
plot(nanmean(dCAPE_VER),0,'o','MarkerSize',10,'MarkerEdgeColor',[0 180/255 1],'MarkerFaceColor',[0 180/255 1],'DisplayName',sprintf('VER: %2.1f %%K^{-1}',nanmean(dCAPE_VER)))
plot(nanmean(dCAPE_LES),0,'o','MarkerSize',10,'MarkerEdgeColor',[0 1 1],'MarkerFaceColor',[0 1 1],'DisplayName',sprintf('LES: %2.1f %%K^{-1}',nanmean(dCAPE_LES)))
plot(CCscale,0,'o','MarkerSize',10,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0],'DisplayName',sprintf('CC: %2.1f %%K^{-1}',CCscale))

set(gca,'FontSize',20)
xlabel('% K^{-1}')
ylabel('Probability Density')
title('(b) CAPE Response to Warming')
legend('location','northwest')

gcfsavepdf('Fig02-CAPE-hist.pdf')







