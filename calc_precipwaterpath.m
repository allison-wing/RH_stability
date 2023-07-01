%% Compute vertically integrated precipitating water paths
% vert. integral of pli + plw

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

%% Read in profiles
for i = 1:length(small_model_list)
    for j = 1:length(SSTs)
        filename = ['/Users/awing/Dropbox/RCEMIP_local/var_files/small_v6/' small_model_list{i} '_RCE_small' num2str(SSTs(j)) filebase]
        ncid = netcdf.open(filename);
        
        if SSTs(j)==295
            prof295small(i).model = small_model_list{i};
            prof295small(i).z = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'zg_avg')); %km
            prof295small(i).p = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'pa_avg')); %hPa
            prof295small(i).ta = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'ta_avg')); % K
            prof295small(i).plwi = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'total_precip')); %g/kg
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
            prof300small(i).plwi = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'total_precip')); %g/kg
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
            prof305small(i).plwi = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'total_precip')); %g/kg
            prof305small(i).Cr = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'Cr'));
            prof305small(i).Cg = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'Cg'));
            prof305small(i).Cb = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'Cb'));
            prof305small(i).color = [prof305small(i).Cr prof305small(i).Cg prof305small(i).Cb];
            colors_model(i,:) = prof305small(i).color;
            netcdf.close(ncid)
        end
    end
    
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
        prof295small(i).plwi = flipud(prof295small(i).plwi);
        prof300small(i).plwi = flipud(prof300small(i).plwi);
        prof305small(i).plwi = flipud(prof305small(i).plwi);
    end
    
    %     %Skip models with missing pressure data
    %     if strcmp(small_model_list(i),'DALES-VER')==1 || ...
    %             strcmp(small_model_list(i),'DALES-LES')==1  || ...
    %             strcmp(small_model_list(i),'DALES-damping-VER')==1     %             
    %     end
    % skip models with missing plw and/or pli data (GEOS, UKMO-GA7.1)
    %Calculate density; rho = p/RT
    prof295small(i).rho = prof295small(i).p*100./(287*prof295small(i).ta);
    prof300small(i).rho = prof300small(i).p*100./(287*prof300small(i).ta);
    prof305small(i).rho = prof305small(i).p*100./(287*prof305small(i).ta);
    
    %Calculate total precipitating water path
    dz = diff(prof295small(i).z*1000); %m
    plwivi295(i) = sum((prof295small(i).plwi(1:end-1)/1000).*prof295small(i).rho(1:end-1).*dz); %kg/m^2 = mm
    
    dz = diff(prof300small(i).z*1000); %m
    plwivi300(i) = sum((prof300small(i).plwi(1:end-1)/1000).*prof300small(i).rho(1:end-1).*dz); %kg/m^2 = mm
    
    dz = diff(prof305small(i).z*1000); %m
    plwivi305(i) = sum((prof305small(i).plwi(1:end-1)/1000).*prof305small(i).rho(1:end-1).*dz); %kg/m^2 = mm
    
end

save('precippath.mat','plwivi295','plwivi300','plwivi305','small_model_list')



