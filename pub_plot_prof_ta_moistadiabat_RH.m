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
    end
    
    %Missing pressure data: DALES, DALES-VER, DALES-LES, DALES-damping,
    %DALES-damping-VER, UCLA-CRM
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
    
    %Compute moist adiabat (wants qhus in g/g, p in Pa), skipping models who are missing data
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
        prof295small(i).tv = nan(length(prof295small(i).ta),1);
        prof300small(i).tv = nan(length(prof300small(i).ta),1);
        prof305small(i).tv = nan(length(prof305small(i).ta),1);
    else
        %convert from specific humidity to mixing ratio
        prof295small(i).rv = prof295small(i).hus/1000./(1-prof295small(i).hus/1000);
        prof300small(i).rv = prof300small(i).hus/1000./(1-prof300small(i).hus/1000);
        prof305small(i).rv = prof305small(i).hus/1000./(1-prof305small(i).hus/1000);
        %compute virtual temperature
        prof295small(i).tv = prof295small(i).ta.*(1+prof295small(i).rv/.622)./(1+prof295small(i).rv);
        prof300small(i).tv = prof300small(i).ta.*(1+prof300small(i).rv/.622)./(1+prof300small(i).rv);
        prof305small(i).tv = prof305small(i).ta.*(1+prof305small(i).rv/.622)./(1+prof305small(i).rv);
    
        %compute moist adiabat, lifting parcel from lowest model level
        [prof295small(i).tma,rv,rl,ri,prof295small(i).tmarho] = atm.calculate_adiabat(prof295small(i).ta(1),prof295small(i).rv(1),prof295small(i).p(1)*100,prof295small(i).p*100);
        [prof300small(i).tma,rv,rl,ri,prof300small(i).tmarho] = atm.calculate_adiabat(prof300small(i).ta(1),prof300small(i).rv(1),prof300small(i).p(1)*100,prof300small(i).p*100);
        [prof305small(i).tma,rv,rl,ri,prof305small(i).tmarho] = atm.calculate_adiabat(prof305small(i).ta(1),prof305small(i).rv(1),prof305small(i).p(1)*100,prof305small(i).p*100);


    end    
end



%% 295K: Plot deviation from moist adiabat
figure('Position',[100 100 1200 2000])
panel1=subplot(2,3,1);
for i = 1:length(small_model_list)
    
    if strcmp(small_model_list(i),'DALES-VER')==1 || strcmp(small_model_list(i),'DALES-LES')==1 || ...
            strcmp(small_model_list(i),'DALES-damping-VER')==1 || ...
            strcmp(small_model_list(i) ,'UKMO-RA1-T')==1
    else
        %plot profile
        
        hold on
        plot(prof295small(i).tmarho-prof295small(i).tv,prof295small(i).z,'Color',prof295small(i).color,'LineWidth',3,'DisplayName',small_model_list{i})
        
    end
   
    set(gca,'FontSize',26)
    xlabel('T_{m,\rho} - T_\rho (K)')
    ylabel('Height (km)')
    ylim([0 20])
    xlim([-22 22])
    vline(0,'k--')
    title('(a) RCE\_small295')
%     legend('location','eastoutside')
   
end

% % set unit for figure size to inches
% set(panel1, 'unit', 'inches');
% % get the original size of figure before the legends are added
% figure_size =  get(panel1, 'position');
% % reset figure position
% figure_size(1) = figure_size(1)-1;
% % set new figure size
% set(panel1, 'position', figure_size)

%% 300K: Plot deviation from moist adiabat

panel2=subplot(2,3,2);
for i = 1:length(small_model_list)
    
    if strcmp(small_model_list(i),'DALES-VER')==1 || strcmp(small_model_list(i),'DALES-LES')==1 || ...
            strcmp(small_model_list(i),'DALES-damping-VER')==1 || ...
            strcmp(small_model_list(i) ,'UKMO-RA1-T')==1
    else
        %plot profile
        
        hold on
        plot(prof300small(i).tmarho-prof300small(i).tv,prof300small(i).z,'Color',prof300small(i).color,'LineWidth',3,'DisplayName',small_model_list{i})
        
    end
    
    set(gca,'FontSize',26)
    xlabel('T_{m,\rho} - T_\rho (K)')
    ylabel('Height (km)')
    ylim([0 20])
    xlim([-22 22])
    vline(0,'k--')
    title('(b) RCE\_small300')
%     legend('location','eastoutside')
end

% % set unit for figure size to inches
% set(panel2, 'unit', 'inches');
% % get the original size of figure before the legends are added
% figure_size =  get(panel2, 'position');
% % reset figure position
% figure_size(1) = figure_size(1)-1;
% % set new figure size
% set(panel2, 'position', figure_size)

%% 305K: Plot deviation from moist adiabat
panel3=subplot(2,3,3);
for i = 1:length(small_model_list)
    
    if strcmp(small_model_list(i),'DALES-VER')==1 || strcmp(small_model_list(i),'DALES-LES')==1 || ...
            strcmp(small_model_list(i),'DALES-damping-VER')==1 || ...
            strcmp(small_model_list(i) ,'UKMO-RA1-T')==1
    else
        %plot profile
        
        hold on
        plot(prof305small(i).tmarho-prof305small(i).tv,prof305small(i).z,'Color',prof305small(i).color,'LineWidth',3,'DisplayName',small_model_list{i})
        
    end

    set(gca,'FontSize',26)
    xlabel('T_{m,\rho} - T_\rho (K)')
    ylabel('Height (km)')
    ylim([0 20])
    xlim([-22 22])
    vline(0,'k--')
    title('(c) RCE\_small305') 
end

% % set unit for figure size to inches
% set(panel3, 'unit', 'inches');
% % get the original size of figure before the legends are added
% figure_size =  get(panel3, 'position');
% % add legends and get its handle
% h_legend = legend('location','eastoutside','FontSize',20);
% % set unit for legend size to inches
% set(h_legend, 'unit', 'inches')
% % get legend size
% legend_size = get(h_legend, 'position');
% % get new size of figure
% figure_size_new =  get(panel3, 'position');
% % reset figure width
% figure_size_new(3) = figure_size(3);
% % figure_size_new(1) = figure_size_new(1)+0.15;
% % set new figure size
% set(panel3, 'position', figure_size_new)

%% 295: Plot relative humidity
panel4 = subplot(2,3,4);
for i = 1:length(small_model_list)
    
    if strcmp(small_model_list(i),'DALES-VER')==1 || strcmp(small_model_list(i),'DALES-LES')==1 || ...
            strcmp(small_model_list(i),'DALES-damping-VER')==1 || ...
            strcmp(small_model_list(i) ,'UKMO-RA1-T')==1
    else
        %plot profile
        
        hold on
        plot(prof295small(i).hur,prof295small(i).z,'Color',prof295small(i).color,'LineWidth',3,'DisplayName',small_model_list{i})
        
    end
   
    set(gca,'FontSize',26)
    xlabel('RH (%)')
    ylabel('Height (km)')
    ylim([0 20])
    xlim([0 120])
    vline(100,'k--')
    title('(d) RCE\_small295')
%     legend('location','eastoutside')
   
end

% % set unit for figure size to inches
% set(panel4, 'unit', 'inches');
% % get the original size of figure before the legends are added
% figure_size =  get(panel4, 'position');
% % reset figure position
% figure_size(1) = figure_size(1)-1;
% % set new figure size
% set(panel4, 'position', figure_size)

%% 300K: Plot relative humidity

panel5=subplot(2,3,5);
for i = 1:length(small_model_list)
    
    if strcmp(small_model_list(i),'DALES-VER')==1 || strcmp(small_model_list(i),'DALES-LES')==1 || ...
            strcmp(small_model_list(i),'DALES-damping-VER')==1 || ...
            strcmp(small_model_list(i) ,'UKMO-RA1-T')==1
    else
        %plot profile
        
        hold on
        plot(prof300small(i).hur,prof300small(i).z,'Color',prof300small(i).color,'LineWidth',3,'DisplayName',small_model_list{i})
        
    end
    
    set(gca,'FontSize',26)
    xlabel('RH (%)')
    ylabel('Height (km)')
    ylim([0 20])
    xlim([0 120])
    vline(100,'k--')
    title('(e) RCE\_small300')
%     legend('location','eastoutside')
end

% % set unit for figure size to inches
% set(panel5, 'unit', 'inches');
% % get the original size of figure before the legends are added
% figure_size =  get(panel5, 'position');
% % reset figure position
% figure_size(1) = figure_size(1)-1;
% % set new figure size
% set(panel5, 'position', figure_size)

%% 305K: Plot relative humidity
panel6=subplot(2,3,6);
for i = 1:length(small_model_list)
    
    if strcmp(small_model_list(i),'DALES-VER')==1 || strcmp(small_model_list(i),'DALES-LES')==1 || ...
            strcmp(small_model_list(i),'DALES-damping-VER')==1 || ...
            strcmp(small_model_list(i) ,'UKMO-RA1-T')==1
    else
        %plot profile
        
        hold on
        plot(prof305small(i).hur,prof305small(i).z,'Color',prof305small(i).color,'LineWidth',3,'DisplayName',small_model_list{i})
        
    end

    set(gca,'FontSize',26)
    xlabel('RH (%)')
    ylabel('Height (km)')
    ylim([0 20])
    xlim([0 120])
    vline(100,'k--')
    title('(f) RCE\_small305') 
end

% % set unit for figure size to inches
% set(panel6, 'unit', 'inches');
% % get the original size of figure before the legends are added
% figure_size =  get(panel6, 'position');
% % add legends and get its handle
% h_legend = legend('location','eastoutside','FontSize',20);
% % set unit for legend size to inches
% set(h_legend, 'unit', 'inches')
% % get legend size
% legend_size = get(h_legend, 'position');
% % get new size of figure
% figure_size_new =  get(panel6, 'position');
% % reset figure width
% figure_size_new(3) = figure_size(3);
% figure_size_new(1) = figure_size_new(1)-1;
% % set new figure size
% set(panel6, 'position', figure_size_new)

gcfsavepdf('Fig-moistadia-tv-rh_nolegend.pdf')






