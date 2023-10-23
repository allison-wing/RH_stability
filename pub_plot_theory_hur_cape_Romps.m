%this version uses the diagnosed PE and epsilon (from full Romps CAPE, not LT Cape from ZBP) to feed back into the Romps
%model to calculate theoretical CAPE, and does so holding different things
%vary or varying across models to explain the intermodel spread
%LIFT FROM LOWEST MODEL LEVEL
%uses simulated gamma at LCL as an input to the Romps model (to diagnose a). Diagnosed epsilon then represents epsilon at LCL.

clear all

%Load data with PE proxy
load PEproxy_cwipwi.mat

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

%heights to average RH over
z1 = 2; %km
z2 = 5; %km

c = atm.load_constants;

%% Read in small profiles
tic;
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
    %Missing pressure data: DALES-VER, DALES-damping-VER
    %missing hus data: DALES-LES
    %Ignore UKMA-RA1-T because it is aggregated
    if strcmp(small_model_list(i),'DALES-VER')==1 || strcmp(small_model_list(i),'DALES-LES')==1 || ...
            strcmp(small_model_list(i),'DALES-damping-VER')==1 || ...
            strcmp(small_model_list(i) ,'UKMO-RA1-T')==1
        %     if strcmp(small_model_list(i),'DALES')==1 || strcmp(small_model_list(i),'DALES-VER')==1 || ...
        %             strcmp(small_model_list(i),'DALES-LES')==1 || strcmp(small_model_list(i),'DALES-damping')==1 || ...
        %             strcmp(small_model_list(i),'DALES-damping-VER')==1 || ...
        %             strcmp(small_model_list(i) ,'UCLA-CRM')==1 || strcmp(small_model_list(i) ,'UKMO-RA1-T')==1
        
        %if missing data, set to nan
        cape295small(i)=NaN; cape300small(i)=NaN; cape305small(i)=NaN;
        Tv_LNB295(i) = NaN; Tv_LNB300(i) = NaN; Tv_LNB305(i) = NaN;
        T_LNB295(i) = NaN; T_LNB300(i) = NaN; T_LNB305(i) = NaN;
        p_LNB295(i) = NaN; p_LNB300(i) = NaN; p_LNB305(i) = NaN;
        T_LCL295(i) = NaN; T_LCL300(i) = NaN; T_LCL305(i) = NaN;
        p_LCL295(i) = NaN; p_LCL300(i) = NaN; p_LCL305(i) = NaN;
        PE_imp295small(i) = NaN; PE_imp300small(i) = NaN; PE_imp305small(i) = NaN;
        eps_imp295small(i) = NaN; eps_imp300small(i) = NaN; eps_imp305small(i) = NaN;
        T1_295(i) = NaN;  T1_300(i) = NaN; T1_305(i) = NaN;
        p1_295(i) = NaN; p1_300(i) = NaN; p1_305(i) = NaN;
        hur295small_avg(i) = NaN; hur300small_avg(i) = NaN; hur305small_avg(i) = NaN;
        prof295small(i).gamma = NaN; prof300small(i).gamma = NaN; prof305small(i).gamma = NaN;
        prof295small(i).qs = NaN; prof300small(i).qs = NaN; prof305small(i).qs = NaN;
        gammaLCL295(i) = NaN;  gammaLCL300(i) = NaN; gammaLCL305(i) = NaN;
        
    else
        %convert from specific humidity to mixing ratio
        prof295small(i).rv = prof295small(i).hus/1000./(1-prof295small(i).hus/1000); %g/g
        prof300small(i).rv = prof300small(i).hus/1000./(1-prof300small(i).hus/1000); %g/g
        prof305small(i).rv = prof305small(i).hus/1000./(1-prof305small(i).hus/1000); %g/g
        
        %compute virtual temperature
        prof295small(i).tv = prof295small(i).ta.*(1+prof295small(i).rv/.622)./(1+prof295small(i).rv); %virtual temperature, K
        prof300small(i).tv = prof300small(i).ta.*(1+prof300small(i).rv/.622)./(1+prof300small(i).rv); %virtual temperature, K
        prof305small(i).tv = prof305small(i).ta.*(1+prof305small(i).rv/.622)./(1+prof305small(i).rv); %virtual temperature, K
        
        %compute density
        prof295small(i).rho = prof295small(i).p*100./(c.Rd*prof295small(i).ta); %kg/m^3
        prof300small(i).rho = prof300small(i).p*100./(c.Rd*prof300small(i).ta); %kg/m^3
        prof305small(i).rho = prof305small(i).p*100./(c.Rd*prof305small(i).ta); %kg/m^3
        
        %calculate saturation specific humidity
        prof295small(i).qs = atm.q_sat(prof295small(i).ta,prof295small(i).p*100); %g/g
        prof300small(i).qs = atm.q_sat(prof300small(i).ta,prof300small(i).p*100); %g/g
        prof305small(i).qs = atm.q_sat(prof305small(i).ta,prof305small(i).p*100); %g/g
        
        %%%%%%%%%%%%%%%%%
        %295%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%
        
        %calculate lower-troposphere indices
        [~,I2] = min(abs(prof295small(i).z-z2));
        [~,I1] = min(abs(prof295small(i).z-z1));
        ILT = false(length(prof295small(i).z),1); ILT(I1:I2) = true;
        
        % Calculate dz, ensuring it has same length as variables
        if prof295small(i).z(1) > 0 %if lowest model level is > 0, add 0 as the first level for dz calcluation
            dz = 1000*diff([0 ;prof295small(i).z]); %m
        else %otherwise add on a level to the top
            topdz = prof295small(i).z(end)-prof295small(i).z(end-1);
            dz = 1000*diff([prof295small(i).z ;prof295small(i).z(end)+topdz]); %m
        end
        
        % Calculate lower-tropospheric humidity
        hur295small_avg(i) = sum(prof295small(i).hur.*prof295small(i).rho.* ILT.*dz)./sum(prof295small(i).rho.* ILT.*dz);
        
        %Lift from lowest model level
        Tinit = prof295small(i).ta(1); %K
        rvinit = prof295small(i).rv(1); %g/g
        pinit = prof295small(i).p(1); %hPa
        
        %Calculate LCL 
        [T_LCL295(i),p_LCL295(i)] = atm.calculate_LCL(Tinit,rvinit,pinit*100); %T_LCL in K, p_LCL in Pa
        
        %Calculate gamma at LCL
        nz = length(prof295small(i).z);
        clear dqsdz
        prof295small(i).qs = smooth(prof295small(i).qs,5); % smooth before taking derivative
        dqsdz(1) = (prof295small(i).qs(2) - prof295small(i).qs(1))/(prof295small(i).z(2) - prof295small(i).z(1)); %forward difference
        dqsdz(2:nz-1) = (prof295small(i).qs(3:nz) - prof295small(i).qs(1:nz-2))./(prof295small(i).z(3:nz) - prof295small(i).z(1:nz-2)); %centered difference
        dqsdz(nz) = (prof295small(i).qs(nz) - prof295small(i).qs(nz-1))/(prof295small(i).z(nz) - prof295small(i).z(nz-1)); %backward difference
        prof295small(i).gamma = -1*dqsdz'./prof295small(i).qs; %km^-1
        prof295small(i).gamma = prof295small(i).gamma/1000; %m^-1
        gammaLCL295(i) = interp1(prof295small(i).p*100,prof295small(i).gamma,p_LCL295(i));
        
        %Calcluate CAPE & LNB 
        [cape295small(i),p_LNB295(i),Tv_LNB295(i),T_LNB295(i)] = calculate_CAPE(Tinit,rvinit,pinit*100,prof295small(i).tv,prof295small(i).p*100,prof295small(i).ta,prof295small(i).z*1000);
        %CAPE in J/kg, p_LNB in Pa, Tv_LNB in K, T_LNB in K
        
        %calculate theory for CAPE based on LCL and LNB, don't plot
        [PE_R16,epsilon_R16,RH_R16_295,CAPE_R16_295] = fn_plot_CAPE_and_RH_RCE(0,T_LCL295(i),p_LCL295(i),T_LNB295(i),gammaLCL295(i));
        
        % Diagnose PE & epsilon implied by theory
        [PE_imp295small(i),eps_imp295small(i)] = find_PE_epsilon(PE_R16,epsilon_R16,double(CAPE_R16_295),double(RH_R16_295),double(cape295small(i)),double(hur295small_avg(i))/100);
        
        
        %%%%%%%%%%%%%%%%%
        %300%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%
        
        %calculate lower-troposphere indices
        [~,I2] = min(abs(prof300small(i).z-z2));
        [~,I1] = min(abs(prof300small(i).z-z1));
        ILT = false(length(prof300small(i).z),1); ILT(I1:I2) = true;
        
        % Calculate dz, ensuring it has same length as variables
        if prof300small(i).z(1) > 0 %if lowest model level is > 0, add 0 as the first level for dz calcluation
            dz = 1000*diff([0 ;prof300small(i).z]); %m
        else %otherwise add on a level to the top
            topdz = prof300small(i).z(end)-prof300small(i).z(end-1);
            dz = 1000*diff([prof300small(i).z ;prof300small(i).z(end)+topdz]); %m
        end
        
        % Calculate lower-tropospheric humidity
        hur300small_avg(i) = sum(prof300small(i).hur.*prof300small(i).rho.* ILT.*dz)./sum(prof300small(i).rho.* ILT.*dz);
        
        %Lift from lowest model level
        Tinit = prof300small(i).ta(1); %K
        rvinit = prof300small(i).rv(1); %g/g
        pinit = prof300small(i).p(1); %hPa
        
        %Calculate LCL 
        [T_LCL300(i),p_LCL300(i)] = atm.calculate_LCL(Tinit,rvinit,pinit*100); %T_LCL in K, p_LCL in Pa
        
        %Calculate gamma at LCL
        nz = length(prof300small(i).z);
        clear dqsdz
        prof300small(i).qs = smooth(prof300small(i).qs,5); % smooth before taking deri
        dqsdz(1) = (prof300small(i).qs(2) - prof300small(i).qs(1))/(prof300small(i).z(2) - prof300small(i).z(1)); %forward difference
        dqsdz(2:nz-1) = (prof300small(i).qs(3:nz) - prof300small(i).qs(1:nz-2))./(prof300small(i).z(3:nz) - prof300small(i).z(1:nz-2)); %centered difference
        dqsdz(nz) = (prof300small(i).qs(nz) - prof300small(i).qs(nz-1))/(prof300small(i).z(nz) - prof300small(i).z(nz-1)); %backward difference
        prof300small(i).gamma = -1*dqsdz'./prof300small(i).qs; %km^-1
        prof300small(i).gamma = prof300small(i).gamma/1000; %m^-1
        gammaLCL300(i) = interp1(prof300small(i).p*100,prof300small(i).gamma,p_LCL300(i));
        
        %Calcluate CAPE & LNB 
        [cape300small(i),p_LNB300(i),Tv_LNB300(i),T_LNB300(i)] = calculate_CAPE(Tinit,rvinit,pinit*100,prof300small(i).tv,prof300small(i).p*100,prof300small(i).ta,prof300small(i).z*1000);
        %CAPE in J/kg, p_LNB in Pa, Tv_LNB in K, T_LNB in K
        
        %calculate theory for CAPE based on LCL and LNB, don't plot
        [PE_R16,epsilon_R16,RH_R16_300,CAPE_R16_300] = fn_plot_CAPE_and_RH_RCE(0,T_LCL300(i),p_LCL300(i),T_LNB300(i),gammaLCL300(i));
        
        % Diagnose PE & epsilon implied by theory
        [PE_imp300small(i),eps_imp300small(i)] = find_PE_epsilon(PE_R16,epsilon_R16,double(CAPE_R16_300),double(RH_R16_300),double(cape300small(i)),double(hur300small_avg(i))/100);
        
        
        %%%%%%%%%%%%%%%%%
        %305%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%
        
        %calculate lower-troposphere indices
        [~,I2] = min(abs(prof305small(i).z-z2));
        [~,I1] = min(abs(prof305small(i).z-z1));
        ILT = false(length(prof305small(i).z),1); ILT(I1:I2) = true;
        
        % Calculate dz, ensuring it has same length as variables
        if prof305small(i).z(1) > 0 %if lowest model level is > 0, add 0 as the first level for dz calcluation
            dz = 1000*diff([0 ;prof305small(i).z]); %m
        else %otherwise add on a level to the top
            topdz = prof305small(i).z(end)-prof305small(i).z(end-1);
            dz = 1000*diff([prof305small(i).z ;prof305small(i).z(end)+topdz]); %m
        end
        
        % Calculate lower-tropospheric humidity
        hur305small_avg(i) = sum(prof305small(i).hur.*prof305small(i).rho.* ILT.*dz)./sum(prof305small(i).rho.* ILT.*dz);
        
        %Lift from lowest model level
        Tinit = prof305small(i).ta(1); %K
        rvinit = prof305small(i).rv(1); %g/g
        pinit = prof305small(i).p(1); %hPs
        
        %Calculate LCL 
        [T_LCL305(i),p_LCL305(i)] = atm.calculate_LCL(Tinit,rvinit,pinit*100); %T_LCL in K, p_LCL in Pa
        
        %Calculate gamma at LCL
        nz = length(prof305small(i).z);
        clear dqsdz
        prof305small(i).qs = smooth(prof305small(i).qs,5); % smooth before taking deri
        dqsdz(1) = (prof305small(i).qs(2) - prof305small(i).qs(1))/(prof305small(i).z(2) - prof305small(i).z(1)); %forward difference
        dqsdz(2:nz-1) = (prof305small(i).qs(3:nz) - prof305small(i).qs(1:nz-2))./(prof305small(i).z(3:nz) - prof305small(i).z(1:nz-2)); %centered difference
        dqsdz(nz) = (prof305small(i).qs(nz) - prof305small(i).qs(nz-1))/(prof305small(i).z(nz) - prof305small(i).z(nz-1)); %backward difference       
        prof305small(i).gamma = -1*dqsdz'./prof305small(i).qs; %km^-1
        prof305small(i).gamma = prof305small(i).gamma/1000; %m^-1
        gammaLCL305(i) = interp1(prof305small(i).p*100,prof305small(i).gamma,p_LCL305(i));
        
        %Calcluate CAPE & LNB 
        [cape305small(i),p_LNB305(i),Tv_LNB305(i),T_LNB305(i)] = calculate_CAPE(Tinit,rvinit,pinit*100,prof305small(i).tv,prof305small(i).p*100,prof305small(i).ta,prof305small(i).z*1000);
        %CAPE in J/kg, p_LNB in Pa, Tv_LNB in K, T_LNB in K
        
        %calculate theory for CAPE based on LCL and LNB, don't plot
        [PE_R16,epsilon_R16,RH_R16_305,CAPE_R16_305] = fn_plot_CAPE_and_RH_RCE(0,T_LCL305(i),p_LCL305(i),T_LNB305(i),gammaLCL305(i));
        
        % Diagnose PE & epsilon implied by theory
        [PE_imp305small(i),eps_imp305small(i)] = find_PE_epsilon(PE_R16,epsilon_R16,double(CAPE_R16_305),double(RH_R16_305),double(cape305small(i)),double(hur305small_avg(i))/100);     
        
    end
    
end

toc;

%% Calculate theoretical CAPE holding different things vary, from diagnosed PE & epsilon, and calculated LCL (cloud base) and LNB (cloud top)
for i = 1:length(small_model_list)
    
    %%295%%
    [CAPE_theory295(i),RH_theory295(i)] = calculate_CAPE_theory(T_LCL295(i),T_LNB295(i),p_LCL295(i),eps_imp295small(i),PE_imp295small(i),gammaLCL295(i),'gamma');
    [CAPE_theory295_vary_depth(i),RH_theory295_vary_depth(i)] = calculate_CAPE_theory(T_LCL295(i),T_LNB295(i),p_LCL295(i),nanmean(eps_imp295small),nanmean(PE_imp295small),nanmean(gammaLCL295),'gamma');
    [CAPE_theory295_vary_PE(i),RH_theory295_vary_PE(i)] = calculate_CAPE_theory(T_LCL295(i),nanmean(T_LNB295),p_LCL295(i),nanmean(eps_imp295small),PE_imp295small(i),nanmean(gammaLCL295),'gamma');
    [CAPE_theory295_vary_eps(i),RH_theory295_vary_eps(i)] = calculate_CAPE_theory(T_LCL295(i),nanmean(T_LNB295),p_LCL295(i),eps_imp295small(i),nanmean(PE_imp295small),nanmean(gammaLCL295),'gamma');
    [CAPE_theory295_vary_gamma(i),RH_theory295_vary_gamma(i)] = calculate_CAPE_theory(T_LCL295(i),nanmean(T_LNB295),p_LCL295(i),nanmean(eps_imp295small),nanmean(PE_imp295small),gammaLCL295(i),'gamma');
    
    %%300%%
    [CAPE_theory300(i),RH_theory300(i)] = calculate_CAPE_theory(T_LCL300(i),T_LNB300(i),p_LCL300(i),eps_imp300small(i),PE_imp300small(i),gammaLCL300(i),'gamma');
    [CAPE_theory300_vary_depth(i),RH_theory300_vary_depth(i)] = calculate_CAPE_theory(T_LCL300(i),T_LNB300(i),p_LCL300(i),nanmean(eps_imp300small),nanmean(PE_imp300small),nanmean(gammaLCL300),'gamma');
    [CAPE_theory300_vary_PE(i),RH_theory300_vary_PE(i)] = calculate_CAPE_theory(T_LCL300(i),nanmean(T_LNB300),p_LCL300(i),nanmean(eps_imp300small),PE_imp300small(i),nanmean(gammaLCL300),'gamma');
    [CAPE_theory300_vary_eps(i),RH_theory300_vary_eps(i)] = calculate_CAPE_theory(T_LCL300(i),nanmean(T_LNB300),p_LCL300(i),eps_imp300small(i),nanmean(PE_imp300small),nanmean(gammaLCL300),'gamma');
    [CAPE_theory300_vary_gamma(i),RH_theory300_vary_gamma(i)] = calculate_CAPE_theory(T_LCL300(i),nanmean(T_LNB300),p_LCL300(i),nanmean(eps_imp300small),nanmean(PE_imp300small),gammaLCL300(i),'gamma');
    
    %%305%%
    [CAPE_theory305(i),RH_theory305(i)] = calculate_CAPE_theory(T_LCL305(i),T_LNB305(i),p_LCL305(i),eps_imp305small(i),PE_imp305small(i),gammaLCL305(i),'gamma');
    [CAPE_theory305_vary_depth(i),RH_theory305_vary_depth(i)] = calculate_CAPE_theory(T_LCL305(i),T_LNB305(i),p_LCL305(i),nanmean(eps_imp305small),nanmean(PE_imp305small),nanmean(gammaLCL305),'gamma');
    [CAPE_theory305_vary_PE(i),RH_theory305_vary_PE(i)] = calculate_CAPE_theory(T_LCL305(i),nanmean(T_LNB305),p_LCL305(i),nanmean(eps_imp305small),PE_imp305small(i),nanmean(gammaLCL305),'gamma');
    [CAPE_theory305_vary_eps(i),RH_theory305_vary_eps(i)] = calculate_CAPE_theory(T_LCL305(i),nanmean(T_LNB305),p_LCL305(i),eps_imp305small(i),nanmean(PE_imp305small),nanmean(gammaLCL305),'gamma');
    [CAPE_theory305_vary_gamma(i),RH_theory305_vary_gamma(i)] = calculate_CAPE_theory(T_LCL305(i),nanmean(T_LNB305),p_LCL305(i),nanmean(eps_imp305small),nanmean(PE_imp305small),gammaLCL305(i),'gamma');
    
    
end

%% Save data
save(['RH_stab_' num2str(z1) '-' num2str(z2) '.mat'],'small_model_list', 'small_model_type', 'iPAR', 'iCRM', 'iVER', 'iLES', 'iEXP', 'symbs', 'sizes',...
    'colors_model_type', 'T1_295', 'p1_295', 'T1_300', 'p1_300','T1_305', 'p1_305',...
    'prof295small', 'prof300small','prof305small','hur295small_avg', 'hur300small_avg', 'hur305small_avg', ...
    'cape295small', 'cape300small', 'cape305small',...
    'PE_imp295small', 'PE_imp300small', 'PE_imp305small', 'eps_imp295small', 'eps_imp300small', 'eps_imp305small', ...
    'T_LCL295','T_LCL300','T_LCL305','p_LCL295','p_LCL300','p_LCL305','Tv_LNB295','Tv_LNB300','Tv_LNB305','T_LNB295','T_LNB300','T_LNB305','p_LNB295','p_LNB300','p_LNB305',...
    'gammaLCL295','gammaLCL300','gammaLCL305')


%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PLOTS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 2 -panel scatter plot
% figure('Position',[100 100 500 1600]) %for 3x1
figure('Position',[100 100 1000 400])

% Panel 1: Param Models, CAPE vs. RH on top of theory, colored by PE proxy
%plot theory
subplot(1,2,1)
fn_plot_CAPE_and_RH_RCE(1,nanmean(T_LCL300(iPAR)),nanmean(p_LCL300(iPAR)),nanmean(T_LNB300(iPAR)),nanmean(gammaLCL300(iPAR)));

cape_param = cape300small(iPAR);
hur_param = hur300small_avg(iPAR);
[r300small_param,p300small_param]=corrcoef(cape_param(~isnan(cape_param)),hur_param(~isnan(cape_param)))
pparam = polyfit(hur_param(~isnan(cape_param))/100,cape_param(~isnan(cape_param)),1);
PE300param = PE300(iPAR);

hold on
scatter(hur_param/100,cape_param,100,PE300(iPAR),'filled','HandleVisibility','off')
scatter(hur_param(isnan(PE300param))/100,cape_param(isnan(PE300param)),100,'k','filled','HandleVisibility','off')
plot(hur_param(~isnan(cape_param))/100,pparam(1)*hur_param(~isnan(cape_param))/100+pparam(2),'Color',[0.5 0 0.5],'LineWidth',3,...
    'DisplayName',sprintf('r = %3.2f; p = %3.2f',r300small_param(1,2),p300small_param(1,2)))

set(gca,'FontSize',20)
legend('location','northwest')
ylabel(['CAPE (J kg^{-1})'])
xlabel(['Lower-Tropospheric Relative Humidity'])
title('(a) Parameterized Convection')
colorbar
caxis([0 0.6e-3])
xlim([0.3 1])

% Panel 2: Explicit Models, CAPE vs. RH on top of theory, colored by PE proxy
%plot theory
subplot(1,2,2)
fn_plot_CAPE_and_RH_RCE(1,nanmean(T_LCL300(iEXP)),nanmean(p_LCL300(iEXP)),nanmean(T_LNB300(iEXP)),nanmean(gammaLCL300(iEXP)));

cape_expl = cape300small(iEXP);
hur_expl = hur300small_avg(iEXP);
[r300small_expl,p300small_expl]=corrcoef(cape_expl(~isnan(cape_expl)),hur_expl(~isnan(cape_expl)))
pexpl = polyfit(hur_expl(~isnan(cape_expl))/100,cape_expl(~isnan(cape_expl)),1);
PE300expl = PE300(iEXP);

hold on
scatter(hur_expl/100,cape_expl,100,PE300(iEXP),'filled','HandleVisibility','off')
scatter(hur_expl(isnan(PE300expl))/100,cape_expl(isnan(PE300expl)),100,'k','filled','HandleVisibility','off')
plot(hur_expl(~isnan(cape_expl))/100,pexpl(1)*hur_expl(~isnan(cape_expl))/100+pexpl(2),'Color',[0.5 0 0.5],'LineWidth',3,...
    'DisplayName',sprintf('r = %3.2f; p = %3.2f',r300small_expl(1,2),p300small_expl(1,2)))

set(gca,'FontSize',20)
legend('location','northwest')
ylabel(['CAPE (J kg^{-1})'])
xlabel(['Lower-Tropospheric Relative Humidity'])
title('(b) Explicit Convection')
colorbar
caxis([0 0.6e-3])
xlim([0.3 1])

gcfsavepdf(['Fig05-theory-hur-cape' num2str(z1) '-' num2str(z2) '-scatter-small300-2panel.pdf'])

%% Seperate out CRM-VER-LES for those models that have all 3. CAPE
% vs. RH on top of theory, colored by PE proxy. Different symbols for
% different models. Different size symbols for simulation type (large =
% CRM, medium = VER, small - LES)

figure('Position',[100 100 650 400])
% subplot(3,1,3)
fn_plot_CAPE_and_RH_RCE(1,nanmean(T_LCL300(iEXP)),nanmean(p_LCL300(iEXP)),nanmean(T_LNB300(iEXP)),nanmean(gammaLCL300(iEXP)));

hold on

%loop over each family of models
    %remove index corresponding to DALES-LES because it is NaN but its CRM
    %partner is present - want to ignore the latter. 
    iignore = find(strcmp(small_model_list(iLES),'DALES-LES')==1);
    iLES(iignore) = NaN;
    iLES = iLES(~isnan(iLES));
    
for i = 1:length(iLES)    
    ifam = [iLES(i) iLES(i)-1 iLES(i)-2]; %grab the LES, the VER (one before)  and the CRM (2 before) for that family
    sz = cell2mat(sizes(ifam)); %grab the sizes
    symb = symbs{ifam}; %grab the symbol (same for each familY)
    
    scatter(hur300small_avg(ifam)/100,cape300small(ifam),sz,PE300(ifam),'filled',symb,'HandleVisibility','off')
end

%legend labels
plot(nan,nan,'d','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0],'MarkerSize',7,'DisplayName','LES')
plot(nan,nan,'d','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0],'MarkerSize',10,'DisplayName','VER')
plot(nan,nan,'d','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0],'MarkerSize',15,'DisplayName','CRM')

%all CRM-VER-LES familes
ifam_all = [iLES iLES-1 iLES-2];
cape_fam = cape300small(ifam_all);
hur_fam = hur300small_avg(ifam_all);
[r300_fam,p300_fam] = corrcoef(cape_fam(~isnan(cape_fam)),hur_fam(~isnan(cape_fam)))
pfam = polyfit(hur_fam(~isnan(cape_fam))/100,cape_fam(~isnan(cape_fam)),1);
plot(hur_fam(~isnan(cape_fam))/100,pfam(1)*hur_fam(~isnan(cape_fam))/100+pfam(2),'Color',[1 0 1],'LineWidth',3,...
    'DisplayName',sprintf('ALL: r = %3.2f; p = %3.2f',r300_fam(1,2),p300_fam(1,2)))

%all LES
cape_LES = cape300small(iLES);
hur_LES = hur300small_avg(iLES);
[r300_LES,p300_LES] = corrcoef(cape_LES(~isnan(cape_LES)),hur_LES(~isnan(cape_LES)))
pLES = polyfit(hur_LES(~isnan(cape_LES))/100,cape_LES(~isnan(cape_LES)),1);
plot(hur_LES(~isnan(cape_LES))/100,pLES(1)*hur_LES(~isnan(cape_LES))/100+pLES(2),'Color',[0.25 0 0.25],'LineWidth',1,...
    'DisplayName',sprintf('LES: r = %3.2f; p = %3.2f',r300_LES(1,2),p300_LES(1,2)))

%all VER (with LES partner)
cape_VER = cape300small(iLES-1);
hur_VER = hur300small_avg(iLES-1);
[r300_VER,p300_VER] = corrcoef(cape_VER(~isnan(cape_VER)),hur_VER(~isnan(cape_VER)))
pVER = polyfit(hur_VER(~isnan(cape_VER))/100,cape_VER(~isnan(cape_VER)),1);
plot(hur_VER(~isnan(cape_VER))/100,pVER(1)*hur_VER(~isnan(cape_VER))/100+pVER(2),'Color',[0.50 0 0.50],'LineWidth',1,...
    'DisplayName',sprintf('VER: r = %3.2f; p = %3.2f',r300_VER(1,2),p300_VER(1,2)))

%all CRM with VER & LES partners
cape_Cles = cape300small(iLES-2);
hur_Cles = hur300small_avg(iLES-2);
[r300_Cles,p300_Cles] = corrcoef(cape_Cles(~isnan(cape_Cles)),hur_Cles(~isnan(cape_Cles)))
pCles = polyfit(hur_Cles(~isnan(cape_Cles))/100,cape_Cles(~isnan(cape_Cles)),1);
plot(hur_Cles(~isnan(cape_Cles))/100,pCles(1)*hur_Cles(~isnan(cape_Cles))/100+pCles(2),'Color',[0.75 0 0.75],'LineWidth',1,...
    'DisplayName',sprintf('CRM: r = %3.2f; p = %3.2f',r300_Cles(1,2),p300_Cles(1,2)))


%leegend labels
iLESs = iLES(~isnan(cape_LES));
ss = symbs(iLESs);
mname = small_model_list(iLESs-2);
for i = 1:length(iLESs)
    plot(nan,nan,ss{i},'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0],'MarkerSize',10,'DisplayName',mname{i})
end

set(gca,'FontSize',20)
legend('location','eastoutside')
ylabel(['CAPE (J kg^{-1})'])
xlabel(['Lower-Tropospheric Relative Humidity'])
title('CRM, VER, LES')
colorbar
caxis([0 0.6e-3])
xlim([0.3 1])


gcfsavepdf(['FigSX-theory-hur-cape' num2str(z1) '-' num2str(z2) '-scatter-small300-verles.pdf'])



%% Intermodel spread in theory-implied PE and proxy PE

PE300par = PE300(iPAR); PE_imp300par = PE_imp300small(iPAR);
PE300crm = PE300(iCRM); PE_imp300crm = PE_imp300small(iCRM);
PE300ver = PE300(iVER); PE_imp300ver = PE_imp300small(iVER);
PE300les = PE300(iLES); PE_imp300les = PE_imp300small(iLES);
PE300exp = PE300(iEXP); PE_imp300exp = PE_imp300small(iEXP);

%correlations
[rPE,pPE] = corrcoef(PE300(logical(~isnan(PE300).*~isnan(PE_imp300small'))),PE_imp300small(logical(~isnan(PE300).*~isnan(PE_imp300small'))))
[rPEparall,pPEparall] = corrcoef(PE300par(logical(~isnan(PE300par).*~isnan(PE_imp300par'))),PE_imp300par(logical(~isnan(PE300par).*~isnan(PE_imp300par'))));
[rPEcrm,pPEcrm] = corrcoef(PE300crm(logical(~isnan(PE300crm).*~isnan(PE_imp300crm'))),PE_imp300crm(logical(~isnan(PE300crm).*~isnan(PE_imp300crm'))));
[rPEver,pPEver] = corrcoef(PE300ver(logical(~isnan(PE300ver).*~isnan(PE_imp300ver'))),PE_imp300ver(logical(~isnan(PE300ver).*~isnan(PE_imp300ver'))));
[rPEles,pPEles] = corrcoef(PE300les(logical(~isnan(PE300les).*~isnan(PE_imp300les'))),PE_imp300les(logical(~isnan(PE300les).*~isnan(PE_imp300les'))));
[rPEexp,pPEexp] = corrcoef(PE300exp(logical(~isnan(PE300exp).*~isnan(PE_imp300exp'))),PE_imp300exp(logical(~isnan(PE300exp).*~isnan(PE_imp300exp'))))

% %correlation ignoring theory-implied param outlier
% PE_imp300par1 = PE_imp300par(PE_imp300par~=nanmax(PE_imp300par));
% PE300par1 = PE300par(PE_imp300par~=nanmax(PE_imp300par));
% [rPEpar,pPEpar] = corrcoef(PE300par1(logical(~isnan(PE300par1).*~isnan(PE_imp300par1'))),PE_imp300par1(logical(~isnan(PE300par1).*~isnan(PE_imp300par1'))));

%correlation ignoring 2 proxy PE outliers
PE300par1 = PE300par(PE300par<0.002);
PE_imp300par1 = PE_imp300par(PE300par<0.002);
[rPEpar,pPEpar] = corrcoef(PE300par1(logical(~isnan(PE300par1).*~isnan(PE_imp300par1'))),PE_imp300par1(logical(~isnan(PE300par1).*~isnan(PE_imp300par1'))));


figure;
% plot(PE300par1,PE_imp300par1,'o','MarkerSize',10,'MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[1 0 0],'DisplayName',sprintf('PAR: r = %3.2f; p = %3.2f',rPEpar(1,2),pPEpar(1,2))) %ignoring  outlier
% plot(PE300par,PE_imp300par,'o','MarkerSize',10,'MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[1 0 0],'DisplayName',sprintf('PAR: r = %3.2f; p = %3.2f',rPEpar(1,2),pPEpar(1,2))) %ignoring  outlier in corr
plot(PE300par,PE_imp300par,'o','MarkerSize',10,'MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[1 0 0],'DisplayName',sprintf('PAR: r = %3.2f; p = %3.2f',rPEparall(1,2),pPEparall(1,2))) 
hold on
plot(PE300(iCRM),PE_imp300small(iCRM),'o','MarkerSize',10,'MarkerEdgeColor',[0 0 1],'MarkerFaceColor',[0 0 1],'DisplayName',sprintf('CRM: r = %3.2f; p = %3.2f',rPEcrm(1,2),pPEcrm(1,2)))
plot(PE300(iVER),PE_imp300small(iVER),'o','MarkerSize',10,'MarkerEdgeColor',[0 180/255 1],'MarkerFaceColor',[0 180/255 1],'DisplayName',sprintf('VER: r = %3.2f; p = %3.2f',rPEver(1,2),pPEver(1,2)))
plot(PE300(iLES),PE_imp300small(iLES),'o','MarkerSize',10,'MarkerEdgeColor',[0 1 1],'MarkerFaceColor',[0 1 1],'DisplayName',sprintf('LES: r = %3.2f; p = %3.2f',rPEles(1,2),pPEles(1,2)))

set(gca,'FontSize',20)
xlabel('Proxy Precip. Efficiency')
ylabel('Theory-Implied Precip. Efficiency')
legend('location','southeast')
% title('Precipitation Efficiency')

% xlim([0 1.5e-3])

gcfsavepdf(['Fig-PE-PEproxy-' num2str(z1) '-' num2str(z2) '.pdf'])

%% change in PE proxy with warming vs. change in theory-implied PE with warming
%theory-implied PE
dPE = (100/10)*(PE_imp305small-PE_imp295small)./PE_imp300small;
dPEpar = dPE(iPAR);
dPEcrm = dPE(iCRM);
dPEver = dPE(iVER);
dPEles = dPE(iLES);
dPEexp = dPE(iEXP);

%proxy PE
dPEP = (100/10)*(PE305-PE295)./PE300;
dPEPpar = dPEP(iPAR);
dPEPcrm = dPEP(iCRM);
dPEPver = dPEP(iVER);
dPEPles = dPEP(iLES);
dPEPexp = dPEP(iEXP);

%correlations
[rPEPEP,pPEPEP] = corrcoef(dPEP(logical(~isnan(dPE').*~isnan(dPEP))),dPE(logical(~isnan(dPE').*~isnan(dPEP))))
[rPEPEPexp,pPEPEPexp] = corrcoef(dPEPexp(logical(~isnan(dPEexp').*~isnan(dPEPexp))),dPEexp(logical(~isnan(dPEexp').*~isnan(dPEPexp))))
[rPEPEPpar,pPEPEPpar] = corrcoef(dPEPpar(logical(~isnan(dPEpar').*~isnan(dPEPpar))),dPEpar(logical(~isnan(dPEpar').*~isnan(dPEPpar))));
[rPEPEPcrm,pPEPEPcrm] = corrcoef(dPEPcrm(logical(~isnan(dPEcrm').*~isnan(dPEPcrm))),dPEcrm(logical(~isnan(dPEcrm').*~isnan(dPEPcrm))));
[rPEPEPver,pPEPEPver] = corrcoef(dPEPver(logical(~isnan(dPEver').*~isnan(dPEPver))),dPEver(logical(~isnan(dPEver').*~isnan(dPEPver))));
[rPEPEPles,pPEPEPles] = corrcoef(dPEPles(logical(~isnan(dPEles').*~isnan(dPEPles))),dPEles(logical(~isnan(dPEles').*~isnan(dPEPles))));

figure;
plot(dPEPpar,dPEpar,'o','MarkerSize',10,'MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[1 0 0],'DisplayName',sprintf('PAR: r = %3.2f; p = %3.2f',rPEPEPpar(1,2),pPEPEPpar(1,2)))
hold on
plot(dPEPcrm,dPEcrm,'o','MarkerSize',10,'MarkerEdgeColor',[0 0 1],'MarkerFaceColor',[0 0 1],'DisplayName',sprintf('CRM: r = %3.2f; p = %3.2f',rPEPEPcrm(1,2),pPEPEPcrm(1,2)))
plot(dPEPver,dPEver,'o','MarkerSize',10,'MarkerEdgeColor',[0 180/255 1],'MarkerFaceColor',[0 180/255 1],'DisplayName',sprintf('VER: r = %3.2f; p = %3.2f',rPEPEPver(1,2),pPEPEPver(1,2)))
plot(dPEPles,dPEles,'o','MarkerSize',10,'MarkerEdgeColor',[0 1 1],'MarkerFaceColor',[0 1 1],'DisplayName',sprintf('LES: r = %3.2f; p = %3.2f',rPEPEPles(1,2),pPEPEPles(1,2)))
hline(0,'k--')
vline(0,'k--')

set(gca,'FontSize',20)
xlabel('Change in Proxy PE (% K^{-1})')
ylabel('Change in Theory-Implied PE (% K^{-1})')
legend('location','southeast')

gcfsavepdf(['Fig-PE-PEproxy-change' num2str(z1) '-' num2str(z2) '.pdf'])

%% Explanations for intermodel spread in CAPE

% figure('Position',[100 100 700 1600]) %for 4x1 subplot
figure('Position',[100 100 1000 1000])
CAPE = cape300small;
CAPEpar = CAPE(iPAR);
CAPEcrm = CAPE(iCRM);
CAPEver = CAPE(iVER);
CAPEles = CAPE(iLES);
CAPEexp = CAPE(iEXP);

% % %%%%%%%%%%%%%%%%
% % %Panel 1: CAPE  vs. theoretical CAPE
% % %%%%%%%%%%%%%%%%
% % CAPEtheory = CAPE_theory300;
% % CAPEtheory_par = CAPEtheory(iPAR);
% % CAPEtheory_crm = CAPEtheory(iCRM);
% % CAPEtheory_ver = CAPEtheory(iVER);
% % CAPEtheory_les = CAPEtheory(iLES);
% % CAPEtheory_exp = CAPEtheory(iEXP);
% % 
% % %correlations
% % [rCT,pCT] = corrcoef(CAPE(logical(~isnan(CAPE).*~isnan(CAPEtheory))),CAPEtheory(logical(~isnan(CAPE).*~isnan(CAPEtheory))))
% % [rCTexp,pCTexp] = corrcoef(CAPEexp(logical(~isnan(CAPEexp).*~isnan(CAPEtheory_exp))),CAPEtheory_exp(logical(~isnan(CAPEexp).*~isnan(CAPEtheory_exp))))
% % [rCTpar,pCTpar] = corrcoef(CAPEpar(logical(~isnan(CAPEpar).*~isnan(CAPEtheory_par))),CAPEtheory_par(logical(~isnan(CAPEpar).*~isnan(CAPEtheory_par))))
% % [rCTcrm,pCTcrm] = corrcoef(CAPEcrm(logical(~isnan(CAPEcrm).*~isnan(CAPEtheory_crm))),CAPEtheory_crm(logical(~isnan(CAPEcrm).*~isnan(CAPEtheory_crm))))
% % [rCTver,pCTver] = corrcoef(CAPEver(logical(~isnan(CAPEver).*~isnan(CAPEtheory_ver))),CAPEtheory_ver(logical(~isnan(CAPEver).*~isnan(CAPEtheory_ver))))
% % [rCTles,pCTles] = corrcoef(CAPEles(logical(~isnan(CAPEles).*~isnan(CAPEtheory_les))),CAPEtheory_les(logical(~isnan(CAPEles).*~isnan(CAPEtheory_les))))
% % 
% % subplot(2,2,1)
% % plot(CAPEtheory_par,CAPEpar,'o','MarkerSize',10,'MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[1 0 0],'DisplayName',sprintf('PAR: r = %3.2f; p = %3.2f',rCTpar(1,2),pCTpar(1,2)))
% % hold on
% % plot(CAPEtheory_crm,CAPEcrm,'o','MarkerSize',10,'MarkerEdgeColor',[0 0 1],'MarkerFaceColor',[0 0 1],'DisplayName',sprintf('CRM: r = %3.2f; p = %3.2f',rCTcrm(1,2),pCTcrm(1,2)))
% % plot(CAPEtheory_ver,CAPEver,'o','MarkerSize',10,'MarkerEdgeColor',[0 180/255 1],'MarkerFaceColor',[0 180/255 1],'DisplayName',sprintf('VER: r = %3.2f; p = %3.2f',rCTver(1,2),pCTver(1,2)))
% % plot(CAPEtheory_les,CAPEles,'o','MarkerSize',10,'MarkerEdgeColor',[0 1 1],'MarkerFaceColor',[0 1 1],'DisplayName',sprintf('LES: r = %3.2f; p = %3.2f',rCTles(1,2),pCTles(1,2)))
% % %plot 1:1 line
% % plot([0:100:4000],[0:100:4000],'k','HandleVisibility','off')
% % 
% % set(gca,'FontSize',20)
% % xlabel('Theoretical CAPE (J kg^{-1})')
% % ylabel('CAPE (J kg^{-1})')
% % title('(a) All Factors')
% % ylim([0 4000])
% % xlim([0 4000])
% % legend('location','northwest')

%%%%%%%%%%%%%%%%
%Panel 2: CAPE vs. theoretical CAPE with vary depth
%%%%%%%%%%%%%%%%
CAPEtheory = CAPE_theory300_vary_depth;
CAPEtheory_par = CAPEtheory(iPAR);
CAPEtheory_crm = CAPEtheory(iCRM);
CAPEtheory_ver = CAPEtheory(iVER);
CAPEtheory_les = CAPEtheory(iLES);
CAPEtheory_exp = CAPEtheory(iEXP);

%correlations
[rCT,pCT] = corrcoef(CAPE(logical(~isnan(CAPE).*~isnan(CAPEtheory))),CAPEtheory(logical(~isnan(CAPE).*~isnan(CAPEtheory))))
[rCTexp,pCTexp] = corrcoef(CAPEexp(logical(~isnan(CAPEexp).*~isnan(CAPEtheory_exp))),CAPEtheory_exp(logical(~isnan(CAPEexp).*~isnan(CAPEtheory_exp))))
[rCTpar,pCTpar] = corrcoef(CAPEpar(logical(~isnan(CAPEpar).*~isnan(CAPEtheory_par))),CAPEtheory_par(logical(~isnan(CAPEpar).*~isnan(CAPEtheory_par))))
[rCTcrm,pCTcrm] = corrcoef(CAPEcrm(logical(~isnan(CAPEcrm).*~isnan(CAPEtheory_crm))),CAPEtheory_crm(logical(~isnan(CAPEcrm).*~isnan(CAPEtheory_crm))))
[rCTver,pCTver] = corrcoef(CAPEver(logical(~isnan(CAPEver).*~isnan(CAPEtheory_ver))),CAPEtheory_ver(logical(~isnan(CAPEver).*~isnan(CAPEtheory_ver))))
[rCTles,pCTles] = corrcoef(CAPEles(logical(~isnan(CAPEles).*~isnan(CAPEtheory_les))),CAPEtheory_les(logical(~isnan(CAPEles).*~isnan(CAPEtheory_les))))

subplot(2,2,2)
plot(CAPEtheory_par,CAPEpar,'o','MarkerSize',10,'MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[1 0 0],'DisplayName',sprintf('PAR: r = %3.2f; p = %3.2f',rCTpar(1,2),pCTpar(1,2)))
hold on
plot(CAPEtheory_crm,CAPEcrm,'o','MarkerSize',10,'MarkerEdgeColor',[0 0 1],'MarkerFaceColor',[0 0 1],'DisplayName',sprintf('CRM: r = %3.2f; p = %3.2f',rCTcrm(1,2),pCTcrm(1,2)))
plot(CAPEtheory_ver,CAPEver,'o','MarkerSize',10,'MarkerEdgeColor',[0 180/255 1],'MarkerFaceColor',[0 180/255 1],'DisplayName',sprintf('VER: r = %3.2f; p = %3.2f',rCTver(1,2),pCTver(1,2)))
plot(CAPEtheory_les,CAPEles,'o','MarkerSize',10,'MarkerEdgeColor',[0 1 1],'MarkerFaceColor',[0 1 1],'DisplayName',sprintf('LES: r = %3.2f; p = %3.2f',rCTles(1,2),pCTles(1,2)))
%plot 1:1 line
plot([0:100:4000],[0:100:4000],'k','HandleVisibility','off')

set(gca,'FontSize',20)
xlabel('Theoretical CAPE (J kg^{-1})')
ylabel('CAPE (J kg^{-1})')
title('(b) Temperature of Convecting Top')
ylim([0 4000])
xlim([0 4000])
legend('location','northwest')

%%%%%%%%%%%%%%%%
%Panel 3 CAPE vs. theoretical CAPE with vary PE
%%%%%%%%%%%%%%%
CAPEtheory = CAPE_theory300_vary_PE;
CAPEtheory_par = CAPEtheory(iPAR);
CAPEtheory_crm = CAPEtheory(iCRM);
CAPEtheory_ver = CAPEtheory(iVER);
CAPEtheory_les = CAPEtheory(iLES);
CAPEtheory_exp = CAPEtheory(iEXP);

%correlations
[rCT,pCT] = corrcoef(CAPE(logical(~isnan(CAPE).*~isnan(CAPEtheory))),CAPEtheory(logical(~isnan(CAPE).*~isnan(CAPEtheory))))
[rCTexp,pCTexp] = corrcoef(CAPEexp(logical(~isnan(CAPEexp).*~isnan(CAPEtheory_exp))),CAPEtheory_exp(logical(~isnan(CAPEexp).*~isnan(CAPEtheory_exp))))
[rCTpar,pCTpar] = corrcoef(CAPEpar(logical(~isnan(CAPEpar).*~isnan(CAPEtheory_par))),CAPEtheory_par(logical(~isnan(CAPEpar).*~isnan(CAPEtheory_par))))
[rCTcrm,pCTcrm] = corrcoef(CAPEcrm(logical(~isnan(CAPEcrm).*~isnan(CAPEtheory_crm))),CAPEtheory_crm(logical(~isnan(CAPEcrm).*~isnan(CAPEtheory_crm))))
[rCTver,pCTver] = corrcoef(CAPEver(logical(~isnan(CAPEver).*~isnan(CAPEtheory_ver))),CAPEtheory_ver(logical(~isnan(CAPEver).*~isnan(CAPEtheory_ver))))
[rCTles,pCTles] = corrcoef(CAPEles(logical(~isnan(CAPEles).*~isnan(CAPEtheory_les))),CAPEtheory_les(logical(~isnan(CAPEles).*~isnan(CAPEtheory_les))))

subplot(2,2,3)
plot(CAPEtheory_par,CAPEpar,'o','MarkerSize',10,'MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[1 0 0],'DisplayName',sprintf('PAR: r = %3.2f; p = %3.2f',rCTpar(1,2),pCTpar(1,2)))
hold on
plot(CAPEtheory_crm,CAPEcrm,'o','MarkerSize',10,'MarkerEdgeColor',[0 0 1],'MarkerFaceColor',[0 0 1],'DisplayName',sprintf('CRM: r = %3.2f; p = %3.2f',rCTcrm(1,2),pCTcrm(1,2)))
plot(CAPEtheory_ver,CAPEver,'o','MarkerSize',10,'MarkerEdgeColor',[0 180/255 1],'MarkerFaceColor',[0 180/255 1],'DisplayName',sprintf('VER: r = %3.2f; p = %3.2f',rCTver(1,2),pCTver(1,2)))
plot(CAPEtheory_les,CAPEles,'o','MarkerSize',10,'MarkerEdgeColor',[0 1 1],'MarkerFaceColor',[0 1 1],'DisplayName',sprintf('LES: r = %3.2f; p = %3.2f',rCTles(1,2),pCTles(1,2)))
%plot 1:1 line
plot([0:100:4000],[0:100:4000],'k','HandleVisibility','off')

set(gca,'FontSize',20)
xlabel('Theoretical CAPE (J kg^{-1})')
ylabel('CAPE (J kg^{-1})')
title('(c) Precip. Efficiency')
ylim([0 4000])
xlim([0 4000])
legend('location','northwest')

%%%%%%%%%%%%%%%%
%Panel 4 CAPE vs. theoretical CAPE with vary PE
%%%%%%%%%%%%%%%
CAPEtheory = CAPE_theory300_vary_eps;
CAPEtheory_par = CAPEtheory(iPAR);
CAPEtheory_crm = CAPEtheory(iCRM);
CAPEtheory_ver = CAPEtheory(iVER);
CAPEtheory_les = CAPEtheory(iLES);
CAPEtheory_exp = CAPEtheory(iEXP);

%correlations
[rCT,pCT] = corrcoef(CAPE(logical(~isnan(CAPE).*~isnan(CAPEtheory))),CAPEtheory(logical(~isnan(CAPE).*~isnan(CAPEtheory))))
[rCTexp,pCTexp] = corrcoef(CAPEexp(logical(~isnan(CAPEexp).*~isnan(CAPEtheory_exp))),CAPEtheory_exp(logical(~isnan(CAPEexp).*~isnan(CAPEtheory_exp))))
[rCTpar,pCTpar] = corrcoef(CAPEpar(logical(~isnan(CAPEpar).*~isnan(CAPEtheory_par))),CAPEtheory_par(logical(~isnan(CAPEpar).*~isnan(CAPEtheory_par))))
[rCTcrm,pCTcrm] = corrcoef(CAPEcrm(logical(~isnan(CAPEcrm).*~isnan(CAPEtheory_crm))),CAPEtheory_crm(logical(~isnan(CAPEcrm).*~isnan(CAPEtheory_crm))))
[rCTver,pCTver] = corrcoef(CAPEver(logical(~isnan(CAPEver).*~isnan(CAPEtheory_ver))),CAPEtheory_ver(logical(~isnan(CAPEver).*~isnan(CAPEtheory_ver))))
[rCTles,pCTles] = corrcoef(CAPEles(logical(~isnan(CAPEles).*~isnan(CAPEtheory_les))),CAPEtheory_les(logical(~isnan(CAPEles).*~isnan(CAPEtheory_les))))

subplot(2,2,4)
plot(CAPEtheory_par,CAPEpar,'o','MarkerSize',10,'MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[1 0 0],'DisplayName',sprintf('PAR: r = %3.2f; p = %3.2f',rCTpar(1,2),pCTpar(1,2)))
hold on
plot(CAPEtheory_crm,CAPEcrm,'o','MarkerSize',10,'MarkerEdgeColor',[0 0 1],'MarkerFaceColor',[0 0 1],'DisplayName',sprintf('CRM: r = %3.2f; p = %3.2f',rCTcrm(1,2),pCTcrm(1,2)))
plot(CAPEtheory_ver,CAPEver,'o','MarkerSize',10,'MarkerEdgeColor',[0 180/255 1],'MarkerFaceColor',[0 180/255 1],'DisplayName',sprintf('VER: r = %3.2f; p = %3.2f',rCTver(1,2),pCTver(1,2)))
plot(CAPEtheory_les,CAPEles,'o','MarkerSize',10,'MarkerEdgeColor',[0 1 1],'MarkerFaceColor',[0 1 1],'DisplayName',sprintf('LES: r = %3.2f; p = %3.2f',rCTles(1,2),pCTles(1,2)))
%plot 1:1 line
plot([0:100:4000],[0:100:4000],'k','HandleVisibility','off')

set(gca,'FontSize',20)
xlabel('Theoretical CAPE (J kg^{-1})')
ylabel('CAPE (J kg^{-1})')
title('(d) Entrainment')
ylim([0 4000])
xlim([0 4000])
legend('location','northwest')

%%%%%%%%%%%%%%%%
%Panel 1 CAPE vs. theoretical CAPE with vary gamma
%%%%%%%%%%%%%%%
CAPEtheory = CAPE_theory300_vary_gamma;
CAPEtheory_par = CAPEtheory(iPAR);
CAPEtheory_crm = CAPEtheory(iCRM);
CAPEtheory_ver = CAPEtheory(iVER);
CAPEtheory_les = CAPEtheory(iLES);
CAPEtheory_exp = CAPEtheory(iEXP);

%correlations
[rCT,pCT] = corrcoef(CAPE(logical(~isnan(CAPE).*~isnan(CAPEtheory))),CAPEtheory(logical(~isnan(CAPE).*~isnan(CAPEtheory))))
[rCTexp,pCTexp] = corrcoef(CAPEexp(logical(~isnan(CAPEexp).*~isnan(CAPEtheory_exp))),CAPEtheory_exp(logical(~isnan(CAPEexp).*~isnan(CAPEtheory_exp))))
[rCTpar,pCTpar] = corrcoef(CAPEpar(logical(~isnan(CAPEpar).*~isnan(CAPEtheory_par))),CAPEtheory_par(logical(~isnan(CAPEpar).*~isnan(CAPEtheory_par))))
[rCTcrm,pCTcrm] = corrcoef(CAPEcrm(logical(~isnan(CAPEcrm).*~isnan(CAPEtheory_crm))),CAPEtheory_crm(logical(~isnan(CAPEcrm).*~isnan(CAPEtheory_crm))))
[rCTver,pCTver] = corrcoef(CAPEver(logical(~isnan(CAPEver).*~isnan(CAPEtheory_ver))),CAPEtheory_ver(logical(~isnan(CAPEver).*~isnan(CAPEtheory_ver))))
[rCTles,pCTles] = corrcoef(CAPEles(logical(~isnan(CAPEles).*~isnan(CAPEtheory_les))),CAPEtheory_les(logical(~isnan(CAPEles).*~isnan(CAPEtheory_les))))

subplot(2,2,1)
plot(CAPEtheory_par,CAPEpar,'o','MarkerSize',10,'MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[1 0 0],'DisplayName',sprintf('PAR: r = %3.2f; p = %3.2f',rCTpar(1,2),pCTpar(1,2)))
hold on
plot(CAPEtheory_crm,CAPEcrm,'o','MarkerSize',10,'MarkerEdgeColor',[0 0 1],'MarkerFaceColor',[0 0 1],'DisplayName',sprintf('CRM: r = %3.2f; p = %3.2f',rCTcrm(1,2),pCTcrm(1,2)))
plot(CAPEtheory_ver,CAPEver,'o','MarkerSize',10,'MarkerEdgeColor',[0 180/255 1],'MarkerFaceColor',[0 180/255 1],'DisplayName',sprintf('VER: r = %3.2f; p = %3.2f',rCTver(1,2),pCTver(1,2)))
plot(CAPEtheory_les,CAPEles,'o','MarkerSize',10,'MarkerEdgeColor',[0 1 1],'MarkerFaceColor',[0 1 1],'DisplayName',sprintf('LES: r = %3.2f; p = %3.2f',rCTles(1,2),pCTles(1,2)))
%plot 1:1 line
plot([0:100:4000],[0:100:4000],'k','HandleVisibility','off')

set(gca,'FontSize',20)
xlabel('Theoretical CAPE (J kg^{-1})')
ylabel('CAPE (J kg^{-1})')
title('(a) \gamma_{LCL}')
ylim([0 4000])
xlim([0 4000])
legend('location','northwest')


% gcfsavepdf('Fig-CAPE-Tlnb-interp.pdf')
gcfsavepdf(['Fig06-CAPE-' num2str(z1) '-' num2str(z2) '.pdf'])


%% Explanations for intermodel spread in CAPE change with warming

%Force cape at 305 to be NaN for DALES and DALES-damping to ignore them for
%changes with warming, because they are aggregated at 305
cape305small(strcmp(small_model_list,'DALES')==1)=NaN;
cape305small(strcmp(small_model_list,'DALES-damping')==1)=NaN;

%compute CC scaling as 100*1/es des/dT = 100*L/(RvT^2) ==> %/K where T = SST
CCscale = 100*c.Lv0/(c.Rv*300^2);

%change in CAPE with warming
dCAPE = (100/10)*(cape305small-cape295small)./cape300small; % percent/K over 10K range
dCAPEpar = dCAPE(iPAR);
dCAPEcrm = dCAPE(iCRM);
dCAPEver = dCAPE(iVER);
dCAPEles = dCAPE(iLES);
dCAPEexp = dCAPE(iEXP);

% figure('Position',[100 100 700 1600]) %for 4x1 subplot
figure('Position',[100 100 1000 1000])

% %%%%%%%%%%%%%%%%
% %Panel 1: CAPE  vs. theoretical CAPE
% %%%%%%%%%%%%%%%%
% dCAPEtheory = (100/10)*(CAPE_theory305-CAPE_theory295)./CAPE_theory300; % percent/K over 10K range
% dCAPEtheory_par = dCAPEtheory(iPAR);
% dCAPEtheory_crm = dCAPEtheory(iCRM);
% dCAPEtheory_ver = dCAPEtheory(iVER);
% dCAPEtheory_les = dCAPEtheory(iLES);
% dCAPEtheory_exp = dCAPEtheory(iEXP);
% 
% %correlations
% [rdCT,pdCT] = corrcoef(dCAPE(logical(~isnan(dCAPE).*~isnan(dCAPEtheory))),dCAPEtheory(logical(~isnan(dCAPE).*~isnan(dCAPEtheory))))
% [rdCTexp,pdCTexp] = corrcoef(dCAPEexp(logical(~isnan(dCAPEexp).*~isnan(dCAPEtheory_exp))),dCAPEtheory_exp(logical(~isnan(dCAPEexp).*~isnan(dCAPEtheory_exp))))
% [rdCTpar,pdCTpar] = corrcoef(dCAPEpar(logical(~isnan(dCAPEpar).*~isnan(dCAPEtheory_par))),dCAPEtheory_par(logical(~isnan(dCAPEpar).*~isnan(dCAPEtheory_par))));
% [rdCTcrm,pdCTcrm] = corrcoef(dCAPEcrm(logical(~isnan(dCAPEcrm).*~isnan(dCAPEtheory_crm))),dCAPEtheory_crm(logical(~isnan(dCAPEcrm).*~isnan(dCAPEtheory_crm))));
% [rdCTver,pdCTver] = corrcoef(dCAPEver(logical(~isnan(dCAPEver).*~isnan(dCAPEtheory_ver))),dCAPEtheory_ver(logical(~isnan(dCAPEver).*~isnan(dCAPEtheory_ver))));
% [rdCTles,pdCTles] = corrcoef(dCAPEles(logical(~isnan(dCAPEles).*~isnan(dCAPEtheory_les))),dCAPEtheory_les(logical(~isnan(dCAPEles).*~isnan(dCAPEtheory_les))));
% 
% subplot(2,2,1)
% plot(dCAPEtheory_par,dCAPEpar,'o','MarkerSize',10,'MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[1 0 0],'DisplayName',sprintf('PAR: r = %3.2f; p = %3.2f',rdCTpar(1,2),pdCTpar(1,2)))
% hold on
% plot(dCAPEtheory_crm,dCAPEcrm,'o','MarkerSize',10,'MarkerEdgeColor',[0 0 1],'MarkerFaceColor',[0 0 1],'DisplayName',sprintf('CRM: r = %3.2f; p = %3.2f',rdCTcrm(1,2),pdCTcrm(1,2)))
% plot(dCAPEtheory_ver,dCAPEver,'o','MarkerSize',10,'MarkerEdgeColor',[0 180/255 1],'MarkerFaceColor',[0 180/255 1],'DisplayName',sprintf('VER: r = %3.2f; p = %3.2f',rdCTver(1,2),pdCTver(1,2)))
% plot(dCAPEtheory_les,dCAPEles,'o','MarkerSize',10,'MarkerEdgeColor',[0 1 1],'MarkerFaceColor',[0 1 1],'DisplayName',sprintf('LES: r = %3.2f; p = %3.2f',rdCTles(1,2),pdCTles(1,2)))
% 
% set(gca,'FontSize',20)
% xlabel('Change in Theoretical CAPE (% K^{-1})')
% ylabel('Change in CAPE (% K^{-1})')
% title('(a) All Factors')
% ylim([0 20])
% xlim([0 20])
% plot([-10:10:30],[-10:10:30],'k','HandleVisibility','off')
% hline(CCscale,'k--')
% vline(0,'k--')
% legend('location','northwest')

%%%%%%%%%%%%%%%%
%Panel 2: CAPE change with warming vs. theoretical CAPE with vary depth change with warming
%%%%%%%%%%%%%%%
dCAPEtheory = (100/10)*(CAPE_theory305_vary_depth-CAPE_theory295_vary_depth)./CAPE_theory300_vary_depth; % percent/K over 10K range
dCAPEtheory_par = dCAPEtheory(iPAR);
dCAPEtheory_crm = dCAPEtheory(iCRM);
dCAPEtheory_ver = dCAPEtheory(iVER);
dCAPEtheory_les = dCAPEtheory(iLES);
dCAPEtheory_exp = dCAPEtheory(iEXP);

%correlations
[rdCT,pdCT] = corrcoef(dCAPE(logical(~isnan(dCAPE).*~isnan(dCAPEtheory))),dCAPEtheory(logical(~isnan(dCAPE).*~isnan(dCAPEtheory))))
[rdCTexp,pdCTexp] = corrcoef(dCAPEexp(logical(~isnan(dCAPEexp).*~isnan(dCAPEtheory_exp))),dCAPEtheory_exp(logical(~isnan(dCAPEexp).*~isnan(dCAPEtheory_exp))))
[rdCTpar,pdCTpar] = corrcoef(dCAPEpar(logical(~isnan(dCAPEpar).*~isnan(dCAPEtheory_par))),dCAPEtheory_par(logical(~isnan(dCAPEpar).*~isnan(dCAPEtheory_par))));
[rdCTcrm,pdCTcrm] = corrcoef(dCAPEcrm(logical(~isnan(dCAPEcrm).*~isnan(dCAPEtheory_crm))),dCAPEtheory_crm(logical(~isnan(dCAPEcrm).*~isnan(dCAPEtheory_crm))));
[rdCTver,pdCTver] = corrcoef(dCAPEver(logical(~isnan(dCAPEver).*~isnan(dCAPEtheory_ver))),dCAPEtheory_ver(logical(~isnan(dCAPEver).*~isnan(dCAPEtheory_ver))));
[rdCTles,pdCTles] = corrcoef(dCAPEles(logical(~isnan(dCAPEles).*~isnan(dCAPEtheory_les))),dCAPEtheory_les(logical(~isnan(dCAPEles).*~isnan(dCAPEtheory_les))));

subplot(2,2,2)
plot(dCAPEtheory_par,dCAPEpar,'o','MarkerSize',10,'MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[1 0 0],'DisplayName',sprintf('PAR: r = %3.2f; p = %3.2f',rdCTpar(1,2),pdCTpar(1,2)))
hold on
plot(dCAPEtheory_crm,dCAPEcrm,'o','MarkerSize',10,'MarkerEdgeColor',[0 0 1],'MarkerFaceColor',[0 0 1],'DisplayName',sprintf('CRM: r = %3.2f; p = %3.2f',rdCTcrm(1,2),pdCTcrm(1,2)))
plot(dCAPEtheory_ver,dCAPEver,'o','MarkerSize',10,'MarkerEdgeColor',[0 180/255 1],'MarkerFaceColor',[0 180/255 1],'DisplayName',sprintf('VER: r = %3.2f; p = %3.2f',rdCTver(1,2),pdCTver(1,2)))
plot(dCAPEtheory_les,dCAPEles,'o','MarkerSize',10,'MarkerEdgeColor',[0 1 1],'MarkerFaceColor',[0 1 1],'DisplayName',sprintf('LES: r = %3.2f; p = %3.2f',rdCTles(1,2),pdCTles(1,2)))

set(gca,'FontSize',20)
xlabel('Change in Theoretical CAPE (% K^{-1})')
ylabel('Change in CAPE (% K^{-1})')
title('(b) Temperature of Convecting Top')
ylim([0 20])
xlim([0 20])
plot([-10:10:30],[-10:10:30],'k','HandleVisibility','off')
hline(CCscale,'k--')
vline(0,'k--')
legend('location','northwest')

%%%%%%%%%%%%%%%%
%Panel 3: CAPE change with warming vs. theoretical CAPE with vary PE change with warming
%%%%%%%%%%%%%%%
dCAPEtheory = (100/10)*(CAPE_theory305_vary_PE-CAPE_theory295_vary_PE)./CAPE_theory300_vary_PE; % percent/K over 10K range
dCAPEtheory_par = dCAPEtheory(iPAR);
dCAPEtheory_crm = dCAPEtheory(iCRM);
dCAPEtheory_ver = dCAPEtheory(iVER);
dCAPEtheory_les = dCAPEtheory(iLES);
dCAPEtheory_exp = dCAPEtheory(iEXP);

%correlations
[rdCT,pdCT] = corrcoef(dCAPE(logical(~isnan(dCAPE).*~isnan(dCAPEtheory))),dCAPEtheory(logical(~isnan(dCAPE).*~isnan(dCAPEtheory))))
[rdCTexp,pdCTexp] = corrcoef(dCAPEexp(logical(~isnan(dCAPEexp).*~isnan(dCAPEtheory_exp))),dCAPEtheory_exp(logical(~isnan(dCAPEexp).*~isnan(dCAPEtheory_exp))))
[rdCTpar,pdCTpar] = corrcoef(dCAPEpar(logical(~isnan(dCAPEpar).*~isnan(dCAPEtheory_par))),dCAPEtheory_par(logical(~isnan(dCAPEpar).*~isnan(dCAPEtheory_par))));
[rdCTcrm,pdCTcrm] = corrcoef(dCAPEcrm(logical(~isnan(dCAPEcrm).*~isnan(dCAPEtheory_crm))),dCAPEtheory_crm(logical(~isnan(dCAPEcrm).*~isnan(dCAPEtheory_crm))));
[rdCTver,pdCTver] = corrcoef(dCAPEver(logical(~isnan(dCAPEver).*~isnan(dCAPEtheory_ver))),dCAPEtheory_ver(logical(~isnan(dCAPEver).*~isnan(dCAPEtheory_ver))));
[rdCTles,pdCTles] = corrcoef(dCAPEles(logical(~isnan(dCAPEles).*~isnan(dCAPEtheory_les))),dCAPEtheory_les(logical(~isnan(dCAPEles).*~isnan(dCAPEtheory_les))));

subplot(2,2,3)
plot(dCAPEtheory_par,dCAPEpar,'o','MarkerSize',10,'MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[1 0 0],'DisplayName',sprintf('PAR: r = %3.2f; p = %3.2f',rdCTpar(1,2),pdCTpar(1,2)))
hold on
plot(dCAPEtheory_crm,dCAPEcrm,'o','MarkerSize',10,'MarkerEdgeColor',[0 0 1],'MarkerFaceColor',[0 0 1],'DisplayName',sprintf('CRM: r = %3.2f; p = %3.2f',rdCTcrm(1,2),pdCTcrm(1,2)))
plot(dCAPEtheory_ver,dCAPEver,'o','MarkerSize',10,'MarkerEdgeColor',[0 180/255 1],'MarkerFaceColor',[0 180/255 1],'DisplayName',sprintf('VER: r = %3.2f; p = %3.2f',rdCTver(1,2),pdCTver(1,2)))
plot(dCAPEtheory_les,dCAPEles,'o','MarkerSize',10,'MarkerEdgeColor',[0 1 1],'MarkerFaceColor',[0 1 1],'DisplayName',sprintf('LES: r = %3.2f; p = %3.2f',rdCTles(1,2),pdCTles(1,2)))

set(gca,'FontSize',20)
xlabel('Change in Theoretical CAPE (% K^{-1})')
ylabel('Change in CAPE (% K^{-1})')
title('(c) Precip. Efficiency')
ylim([0 20])
xlim([0 20])
plot([-10:10:30],[-10:10:30],'k','HandleVisibility','off')
hline(CCscale,'k--')
vline(0,'k--')
legend('location','northwest')

%%%%%%%%%%%%%%%%
%Panel 4: CAPE change with warming vs. theoretical CAPE with vary entrainment change with warming
%%%%%%%%%%%%%%%
dCAPEtheory = (100/10)*(CAPE_theory305_vary_eps-CAPE_theory295_vary_eps)./CAPE_theory300_vary_eps; % percent/K over 10K range
dCAPEtheory_par = dCAPEtheory(iPAR);
dCAPEtheory_crm = dCAPEtheory(iCRM);
dCAPEtheory_ver = dCAPEtheory(iVER);
dCAPEtheory_les = dCAPEtheory(iLES);
dCAPEtheory_exp = dCAPEtheory(iEXP);

%correlations
[rdCT,pdCT] = corrcoef(dCAPE(logical(~isnan(dCAPE).*~isnan(dCAPEtheory))),dCAPEtheory(logical(~isnan(dCAPE).*~isnan(dCAPEtheory))))
[rdCTexp,pdCTexp] = corrcoef(dCAPEexp(logical(~isnan(dCAPEexp).*~isnan(dCAPEtheory_exp))),dCAPEtheory_exp(logical(~isnan(dCAPEexp).*~isnan(dCAPEtheory_exp))))
[rdCTpar,pdCTpar] = corrcoef(dCAPEpar(logical(~isnan(dCAPEpar).*~isnan(dCAPEtheory_par))),dCAPEtheory_par(logical(~isnan(dCAPEpar).*~isnan(dCAPEtheory_par))));
[rdCTcrm,pdCTcrm] = corrcoef(dCAPEcrm(logical(~isnan(dCAPEcrm).*~isnan(dCAPEtheory_crm))),dCAPEtheory_crm(logical(~isnan(dCAPEcrm).*~isnan(dCAPEtheory_crm))));
[rdCTver,pdCTver] = corrcoef(dCAPEver(logical(~isnan(dCAPEver).*~isnan(dCAPEtheory_ver))),dCAPEtheory_ver(logical(~isnan(dCAPEver).*~isnan(dCAPEtheory_ver))));
[rdCTles,pdCTles] = corrcoef(dCAPEles(logical(~isnan(dCAPEles).*~isnan(dCAPEtheory_les))),dCAPEtheory_les(logical(~isnan(dCAPEles).*~isnan(dCAPEtheory_les))));

subplot(2,2,4)
plot(dCAPEtheory_par,dCAPEpar,'o','MarkerSize',10,'MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[1 0 0],'DisplayName',sprintf('PAR: r = %3.2f; p = %3.2f',rdCTpar(1,2),pdCTpar(1,2)))
hold on
plot(dCAPEtheory_crm,dCAPEcrm,'o','MarkerSize',10,'MarkerEdgeColor',[0 0 1],'MarkerFaceColor',[0 0 1],'DisplayName',sprintf('CRM: r = %3.2f; p = %3.2f',rdCTcrm(1,2),pdCTcrm(1,2)))
plot(dCAPEtheory_ver,dCAPEver,'o','MarkerSize',10,'MarkerEdgeColor',[0 180/255 1],'MarkerFaceColor',[0 180/255 1],'DisplayName',sprintf('VER: r = %3.2f; p = %3.2f',rdCTver(1,2),pdCTver(1,2)))
plot(dCAPEtheory_les,dCAPEles,'o','MarkerSize',10,'MarkerEdgeColor',[0 1 1],'MarkerFaceColor',[0 1 1],'DisplayName',sprintf('LES: r = %3.2f; p = %3.2f',rdCTles(1,2),pdCTles(1,2)))

set(gca,'FontSize',20)
xlabel('Change in Theoretical CAPE (% K^{-1})')
ylabel('Change in CAPE (% K^{-1})')
title('(d) Entrainment')
ylim([0 20])
xlim([0 20])
plot([-10:10:30],[-10:10:30],'k','HandleVisibility','off')
hline(CCscale,'k--')
vline(0,'k--')
legend('location','northwest')

%%%%%%%%%%%%%%%%
%Panel 1: CAPE change with warming vs. theoretical CAPE with vary gamma change with warming
%%%%%%%%%%%%%%%
dCAPEtheory = (100/10)*(CAPE_theory305_vary_gamma-CAPE_theory295_vary_gamma)./CAPE_theory300_vary_gamma; % percent/K over 10K range
dCAPEtheory_par = dCAPEtheory(iPAR);
dCAPEtheory_crm = dCAPEtheory(iCRM);
dCAPEtheory_ver = dCAPEtheory(iVER);
dCAPEtheory_les = dCAPEtheory(iLES);
dCAPEtheory_exp = dCAPEtheory(iEXP);

%correlations
[rdCT,pdCT] = corrcoef(dCAPE(logical(~isnan(dCAPE).*~isnan(dCAPEtheory))),dCAPEtheory(logical(~isnan(dCAPE).*~isnan(dCAPEtheory))))
[rdCTexp,pdCTexp] = corrcoef(dCAPEexp(logical(~isnan(dCAPEexp).*~isnan(dCAPEtheory_exp))),dCAPEtheory_exp(logical(~isnan(dCAPEexp).*~isnan(dCAPEtheory_exp))))
[rdCTpar,pdCTpar] = corrcoef(dCAPEpar(logical(~isnan(dCAPEpar).*~isnan(dCAPEtheory_par))),dCAPEtheory_par(logical(~isnan(dCAPEpar).*~isnan(dCAPEtheory_par))));
[rdCTcrm,pdCTcrm] = corrcoef(dCAPEcrm(logical(~isnan(dCAPEcrm).*~isnan(dCAPEtheory_crm))),dCAPEtheory_crm(logical(~isnan(dCAPEcrm).*~isnan(dCAPEtheory_crm))));
[rdCTver,pdCTver] = corrcoef(dCAPEver(logical(~isnan(dCAPEver).*~isnan(dCAPEtheory_ver))),dCAPEtheory_ver(logical(~isnan(dCAPEver).*~isnan(dCAPEtheory_ver))));
[rdCTles,pdCTles] = corrcoef(dCAPEles(logical(~isnan(dCAPEles).*~isnan(dCAPEtheory_les))),dCAPEtheory_les(logical(~isnan(dCAPEles).*~isnan(dCAPEtheory_les))));

subplot(2,2,1)
plot(dCAPEtheory_par,dCAPEpar,'o','MarkerSize',10,'MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[1 0 0],'DisplayName',sprintf('PAR: r = %3.2f; p = %3.2f',rdCTpar(1,2),pdCTpar(1,2)))
hold on
plot(dCAPEtheory_crm,dCAPEcrm,'o','MarkerSize',10,'MarkerEdgeColor',[0 0 1],'MarkerFaceColor',[0 0 1],'DisplayName',sprintf('CRM: r = %3.2f; p = %3.2f',rdCTcrm(1,2),pdCTcrm(1,2)))
plot(dCAPEtheory_ver,dCAPEver,'o','MarkerSize',10,'MarkerEdgeColor',[0 180/255 1],'MarkerFaceColor',[0 180/255 1],'DisplayName',sprintf('VER: r = %3.2f; p = %3.2f',rdCTver(1,2),pdCTver(1,2)))
plot(dCAPEtheory_les,dCAPEles,'o','MarkerSize',10,'MarkerEdgeColor',[0 1 1],'MarkerFaceColor',[0 1 1],'DisplayName',sprintf('LES: r = %3.2f; p = %3.2f',rdCTles(1,2),pdCTles(1,2)))

set(gca,'FontSize',20)
xlabel('Change in Theoretical CAPE (% K^{-1})')
ylabel('Change in CAPE (% K^{-1})')
title('(a) \gamma_{LCL}')
ylim([0 20])
xlim([0 20])
plot([-10:10:30],[-10:10:30],'k','HandleVisibility','off')
hline(CCscale,'k--')
vline(0,'k--')
legend('location','northwest')

% gcfsavepdf('Fig-CAPE-Tlnb-interp.pdf')
gcfsavepdf(['Fig08-CAPEchange-' num2str(z1) '-' num2str(z2) '.pdf'])

%% COMBINED PLOT FOR Explanations for intermodel spread in RH AND ITS CHANGES WITH WARMING

% figure('Position',[100 100 700 1600]) %for 4x1 subplot
figure('Position',[100 100 1000 1000])
RH = hur300small_avg;
RHpar = hur300small_avg(iPAR);
RHcrm = hur300small_avg(iCRM);
RHver = hur300small_avg(iVER);
RHles = hur300small_avg(iLES);
RHexp = hur300small_avg(iEXP);

%Explanations for intermodel spread in RH
%%%%%%%%%%%%%%%%
%Panel 1: RH  vs. theoretical RH with vary PE
%%%%%%%%%%%%%%%%
RHtheory = RH_theory300_vary_PE;
RHtheory_par = RHtheory(iPAR);
RHtheory_crm = RHtheory(iCRM);
RHtheory_ver = RHtheory(iVER);
RHtheory_les = RHtheory(iLES);
RHtheory_exp = RHtheory(iEXP);

%correlations
[rRT,pRT] = corrcoef(RH(logical(~isnan(RH).*~isnan(RHtheory))),RHtheory(logical(~isnan(RH).*~isnan(RHtheory))))
[rRTexp,pRTexp] = corrcoef(RHexp(logical(~isnan(RHexp).*~isnan(RHtheory_exp))),RHtheory_exp(logical(~isnan(RHexp).*~isnan(RHtheory_exp))))
[rRTpar,pRTpar] = corrcoef(RHpar(logical(~isnan(RHpar).*~isnan(RHtheory_par))),RHtheory_par(logical(~isnan(RHpar).*~isnan(RHtheory_par))));
[rRTcrm,pRTcrm] = corrcoef(RHcrm(logical(~isnan(RHcrm).*~isnan(RHtheory_crm))),RHtheory_crm(logical(~isnan(RHcrm).*~isnan(RHtheory_crm))));
[rRTver,pRTver] = corrcoef(RHver(logical(~isnan(RHver).*~isnan(RHtheory_ver))),RHtheory_ver(logical(~isnan(RHver).*~isnan(RHtheory_ver))));
[rRTles,pRTles] = corrcoef(RHles(logical(~isnan(RHles).*~isnan(RHtheory_les))),RHtheory_les(logical(~isnan(RHles).*~isnan(RHtheory_les))));

subplot(2,2,1)
plot(RHtheory_par*100,RHpar,'o','MarkerSize',10,'MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[1 0 0],'DisplayName',sprintf('PAR: r = %3.2f; p = %3.2f',rRTpar(1,2),pRTpar(1,2)))
hold on
plot(RHtheory_crm*100,RHcrm,'o','MarkerSize',10,'MarkerEdgeColor',[0 0 1],'MarkerFaceColor',[0 0 1],'DisplayName',sprintf('CRM: r = %3.2f; p = %3.2f',rRTcrm(1,2),pRTcrm(1,2)))
plot(RHtheory_ver*100,RHver,'o','MarkerSize',10,'MarkerEdgeColor',[0 180/255 1],'MarkerFaceColor',[0 180/255 1],'DisplayName',sprintf('VER: r = %3.2f; p = %3.2f',rRTver(1,2),pRTver(1,2)))
plot(RHtheory_les*100,RHles,'o','MarkerSize',10,'MarkerEdgeColor',[0 1 1],'MarkerFaceColor',[0 1 1],'DisplayName',sprintf('LES: r = %3.2f; p = %3.2f',rRTles(1,2),pRTles(1,2)))

set(gca,'FontSize',20)
xlabel('Theoretical RH')
ylabel('RH')
title('(a) Precip. Efficiency')
xlim([30 100])
ylim([30 100])
plot([30:10:100],[30:10:100],'k','HandleVisibility','off')
legend('location','northwest')

%%%%%%%%%%%%%%%%
%Panel 2: RH  vs. theoretical RH with vary entrainment
%%%%%%%%%%%%%%%%
RHtheory = RH_theory300_vary_eps;
RHtheory_par = RHtheory(iPAR);
RHtheory_crm = RHtheory(iCRM);
RHtheory_ver = RHtheory(iVER);
RHtheory_les = RHtheory(iLES);
RHtheory_exp = RHtheory(iEXP);

%correlations
[rRT,pRT] = corrcoef(RH(logical(~isnan(RH).*~isnan(RHtheory))),RHtheory(logical(~isnan(RH).*~isnan(RHtheory))))
[rRTexp,pRTexp] = corrcoef(RHexp(logical(~isnan(RHexp).*~isnan(RHtheory_exp))),RHtheory_exp(logical(~isnan(RHexp).*~isnan(RHtheory_exp))))
[rRTpar,pRTpar] = corrcoef(RHpar(logical(~isnan(RHpar).*~isnan(RHtheory_par))),RHtheory_par(logical(~isnan(RHpar).*~isnan(RHtheory_par))));
[rRTcrm,pRTcrm] = corrcoef(RHcrm(logical(~isnan(RHcrm).*~isnan(RHtheory_crm))),RHtheory_crm(logical(~isnan(RHcrm).*~isnan(RHtheory_crm))));
[rRTver,pRTver] = corrcoef(RHver(logical(~isnan(RHver).*~isnan(RHtheory_ver))),RHtheory_ver(logical(~isnan(RHver).*~isnan(RHtheory_ver))));
[rRTles,pRTles] = corrcoef(RHles(logical(~isnan(RHles).*~isnan(RHtheory_les))),RHtheory_les(logical(~isnan(RHles).*~isnan(RHtheory_les))));

subplot(2,2,2)
plot(RHtheory_par*100,RHpar,'o','MarkerSize',10,'MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[1 0 0],'DisplayName',sprintf('PAR: r = %3.2f; p = %3.2f',rRTpar(1,2),pRTpar(1,2)))
hold on
plot(RHtheory_crm*100,RHcrm,'o','MarkerSize',10,'MarkerEdgeColor',[0 0 1],'MarkerFaceColor',[0 0 1],'DisplayName',sprintf('CRM: r = %3.2f; p = %3.2f',rRTcrm(1,2),pRTcrm(1,2)))
plot(RHtheory_ver*100,RHver,'o','MarkerSize',10,'MarkerEdgeColor',[0 180/255 1],'MarkerFaceColor',[0 180/255 1],'DisplayName',sprintf('VER: r = %3.2f; p = %3.2f',rRTver(1,2),pRTver(1,2)))
plot(RHtheory_les*100,RHles,'o','MarkerSize',10,'MarkerEdgeColor',[0 1 1],'MarkerFaceColor',[0 1 1],'DisplayName',sprintf('LES: r = %3.2f; p = %3.2f',rRTles(1,2),pRTles(1,2)))

set(gca,'FontSize',20)
xlabel('Theoretical RH')
ylabel('RH')
title('(b) Entrainment')
xlim([30 100])
ylim([30 100])
plot([30:10:100],[30:10:100],'k','HandleVisibility','off')
legend('location','northwest')


% Explanations for intermodel spread in RH change with warming
%Force RH at 305 to be NaN for DALES and DALES-damping to ignore them for
%changes with warming, because they are aggregated at 305
hur305small_avg(strcmp(small_model_list,'DALES')==1)=NaN;
hur305small_avg(strcmp(small_model_list,'DALES-damping')==1)=NaN;

%change in RH with warming
dRH = (100/10)*(hur305small_avg-hur295small_avg)./hur300small_avg; % percent/K over 10K range
dRHpar = dRH(iPAR);
dRHcrm = dRH(iCRM);
dRHver = dRH(iVER);
dRHles = dRH(iLES);
dRHexp = dRH(iEXP);


%%%%%%%%%%%%%%%%
%Panel 3: RH change with warming vs. theoretical RH with vary PE change with warming
%%%%%%%%%%%%%%%%
dRHtheory = (100/10)*(RH_theory305_vary_PE-RH_theory295_vary_PE)./RH_theory300_vary_PE; % percent/K over 10K range;
dRHtheory_par = dRHtheory(iPAR);
dRHtheory_crm = dRHtheory(iCRM);
dRHtheory_ver = dRHtheory(iVER);
dRHtheory_les = dRHtheory(iLES);
dRHtheory_exp = dRHtheory(iEXP);

%correlations
[rdRT,pdRT] = corrcoef(dRH(logical(~isnan(dRH).*~isnan(dRHtheory))),dRHtheory(logical(~isnan(dRH).*~isnan(dRHtheory))))
[rdRTexp,pdRTexp] = corrcoef(dRHexp(logical(~isnan(dRHexp).*~isnan(dRHtheory_exp))),dRHtheory_exp(logical(~isnan(dRHexp).*~isnan(dRHtheory_exp))))
[rdRTpar,pdRTpar] = corrcoef(dRHpar(logical(~isnan(dRHpar).*~isnan(dRHtheory_par))),dRHtheory_par(logical(~isnan(dRHpar).*~isnan(dRHtheory_par))));
[rdRTcrm,pdRTcrm] = corrcoef(dRHcrm(logical(~isnan(dRHcrm).*~isnan(dRHtheory_crm))),dRHtheory_crm(logical(~isnan(dRHcrm).*~isnan(dRHtheory_crm))));
[rdRTver,pdRTver] = corrcoef(dRHver(logical(~isnan(dRHver).*~isnan(dRHtheory_ver))),dRHtheory_ver(logical(~isnan(dRHver).*~isnan(dRHtheory_ver))));
[rdRTles,pdRTles] = corrcoef(dRHles(logical(~isnan(dRHles).*~isnan(dRHtheory_les))),dRHtheory_les(logical(~isnan(dRHles).*~isnan(dRHtheory_les))));

subplot(2,2,3)
plot(dRHtheory_par,dRHpar,'o','MarkerSize',10,'MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[1 0 0],'DisplayName',sprintf('PAR: r = %3.2f; p = %3.2f',rdRTpar(1,2),pdRTpar(1,2)))
hold on
plot(dRHtheory_crm,dRHcrm,'o','MarkerSize',10,'MarkerEdgeColor',[0 0 1],'MarkerFaceColor',[0 0 1],'DisplayName',sprintf('CRM: r = %3.2f; p = %3.2f',rdRTcrm(1,2),pdRTcrm(1,2)))
plot(dRHtheory_ver,dRHver,'o','MarkerSize',10,'MarkerEdgeColor',[0 180/255 1],'MarkerFaceColor',[0 180/255 1],'DisplayName',sprintf('VER: r = %3.2f; p = %3.2f',rdRTver(1,2),pdRTver(1,2)))
plot(dRHtheory_les,dRHles,'o','MarkerSize',10,'MarkerEdgeColor',[0 1 1],'MarkerFaceColor',[0 1 1],'DisplayName',sprintf('LES: r = %3.2f; p = %3.2f',rdRTles(1,2),pdRTles(1,2)))

set(gca,'FontSize',20)
xlabel('Change in Theoretical RH (% K^{-1})')
ylabel('Change in RH (% K^{-1})')
title('(c) Precipitation Efficiency')
xlim([-1 3])
ylim([-1 3])
% xlim([-5 5])
% ylim([-5 5])
vline(0,'k--')
hline(0,'k--')
plot([-4:1:4],[-4:1:4],'k','HandleVisibility','off')
legend('location','southeast')

%%%%%%%%%%%%%%%%
%Panel 4: RH change with warming vs. theoretical RH with vary entrainment change with warming
%%%%%%%%%%%%%%%%
dRHtheory = (100/10)*(RH_theory305_vary_eps-RH_theory295_vary_eps)./RH_theory300_vary_eps; % percent/K over 10K range;
dRHtheory_par = dRHtheory(iPAR);
dRHtheory_crm = dRHtheory(iCRM);
dRHtheory_ver = dRHtheory(iVER);
dRHtheory_les = dRHtheory(iLES);
dRHtheory_exp = dRHtheory(iEXP);

%correlations
[rdRT,pdRT] = corrcoef(dRH(logical(~isnan(dRH).*~isnan(dRHtheory))),dRHtheory(logical(~isnan(dRH).*~isnan(dRHtheory))))
[rdRTexp,pdRTexp] = corrcoef(dRHexp(logical(~isnan(dRHexp).*~isnan(dRHtheory_exp))),dRHtheory_exp(logical(~isnan(dRHexp).*~isnan(dRHtheory_exp))))
[rdRTpar,pdRTpar] = corrcoef(dRHpar(logical(~isnan(dRHpar).*~isnan(dRHtheory_par))),dRHtheory_par(logical(~isnan(dRHpar).*~isnan(dRHtheory_par))));
[rdRTcrm,pdRTcrm] = corrcoef(dRHcrm(logical(~isnan(dRHcrm).*~isnan(dRHtheory_crm))),dRHtheory_crm(logical(~isnan(dRHcrm).*~isnan(dRHtheory_crm))));
[rdRTver,pdRTver] = corrcoef(dRHver(logical(~isnan(dRHver).*~isnan(dRHtheory_ver))),dRHtheory_ver(logical(~isnan(dRHver).*~isnan(dRHtheory_ver))));
[rdRTles,pdRTles] = corrcoef(dRHles(logical(~isnan(dRHles).*~isnan(dRHtheory_les))),dRHtheory_les(logical(~isnan(dRHles).*~isnan(dRHtheory_les))));

subplot(2,2,4)
plot(dRHtheory_par,dRHpar,'o','MarkerSize',10,'MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[1 0 0],'DisplayName',sprintf('PAR: r = %3.2f; p = %3.2f',rdRTpar(1,2),pdRTpar(1,2)))
hold on
plot(dRHtheory_crm,dRHcrm,'o','MarkerSize',10,'MarkerEdgeColor',[0 0 1],'MarkerFaceColor',[0 0 1],'DisplayName',sprintf('CRM: r = %3.2f; p = %3.2f',rdRTcrm(1,2),pdRTcrm(1,2)))
plot(dRHtheory_ver,dRHver,'o','MarkerSize',10,'MarkerEdgeColor',[0 180/255 1],'MarkerFaceColor',[0 180/255 1],'DisplayName',sprintf('VER: r = %3.2f; p = %3.2f',rdRTver(1,2),pdRTver(1,2)))
plot(dRHtheory_les,dRHles,'o','MarkerSize',10,'MarkerEdgeColor',[0 1 1],'MarkerFaceColor',[0 1 1],'DisplayName',sprintf('LES: r = %3.2f; p = %3.2f',rdRTles(1,2),pdRTles(1,2)))

set(gca,'FontSize',20)
xlabel('Change in Theoretical RH (% K^{-1})')
ylabel('Change in RH (% K^{-1})')
title('(d) Entrainment')
xlim([-1 3])
ylim([-1 3])
% xlim([-5 5])
% ylim([-5 5])
vline(0,'k--')
hline(0,'k--')
plot([-4:1:4],[-4:1:4],'k','HandleVisibility','off')
legend('location','southeast')


gcfsavepdf(['Fig07-RHRHchange-' num2str(z1) '-' num2str(z2) '.pdf'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %% Explanations for intermodel spread in RH
% % figure('Position',[100 100 700 1600]) %for 4x1 subplot
% figure('Position',[100 100 1000 1000])
% RH = hur300small_avg;
% RHpar = hur300small_avg(iPAR);
% RHcrm = hur300small_avg(iCRM);
% RHver = hur300small_avg(iVER);
% RHles = hur300small_avg(iLES);
% RHexp = hur300small_avg(iEXP);
% 
% %%%%%%%%%%%%%%%%
% %Panel 1: RH  vs. theoretical RH with vary gamma
% %%%%%%%%%%%%%%%%
% RHtheory = RH_theory300_vary_gamma;
% RHtheory_par = RHtheory(iPAR);
% RHtheory_crm = RHtheory(iCRM);
% RHtheory_ver = RHtheory(iVER);
% RHtheory_les = RHtheory(iLES);
% RHtheory_exp = RHtheory(iEXP);
% 
% %correlations
% [rRT,pRT] = corrcoef(RH(logical(~isnan(RH).*~isnan(RHtheory))),RHtheory(logical(~isnan(RH).*~isnan(RHtheory))))
% [rRTexp,pRTexp] = corrcoef(RHexp(logical(~isnan(RHexp).*~isnan(RHtheory_exp))),RHtheory_exp(logical(~isnan(RHexp).*~isnan(RHtheory_exp))))
% [rRTpar,pRTpar] = corrcoef(RHpar(logical(~isnan(RHpar).*~isnan(RHtheory_par))),RHtheory_par(logical(~isnan(RHpar).*~isnan(RHtheory_par))));
% [rRTcrm,pRTcrm] = corrcoef(RHcrm(logical(~isnan(RHcrm).*~isnan(RHtheory_crm))),RHtheory_crm(logical(~isnan(RHcrm).*~isnan(RHtheory_crm))));
% [rRTver,pRTver] = corrcoef(RHver(logical(~isnan(RHver).*~isnan(RHtheory_ver))),RHtheory_ver(logical(~isnan(RHver).*~isnan(RHtheory_ver))));
% [rRTles,pRTles] = corrcoef(RHles(logical(~isnan(RHles).*~isnan(RHtheory_les))),RHtheory_les(logical(~isnan(RHles).*~isnan(RHtheory_les))));
% 
% subplot(2,2,1)
% plot(RHtheory_par*100,RHpar,'o','MarkerSize',10,'MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[1 0 0],'DisplayName',sprintf('PAR: r = %3.2f; p = %3.2f',rRTpar(1,2),pRTpar(1,2)))
% hold on
% plot(RHtheory_crm*100,RHcrm,'o','MarkerSize',10,'MarkerEdgeColor',[0 0 1],'MarkerFaceColor',[0 0 1],'DisplayName',sprintf('CRM: r = %3.2f; p = %3.2f',rRTcrm(1,2),pRTcrm(1,2)))
% plot(RHtheory_ver*100,RHver,'o','MarkerSize',10,'MarkerEdgeColor',[0 180/255 1],'MarkerFaceColor',[0 180/255 1],'DisplayName',sprintf('VER: r = %3.2f; p = %3.2f',rRTver(1,2),pRTver(1,2)))
% plot(RHtheory_les*100,RHles,'o','MarkerSize',10,'MarkerEdgeColor',[0 1 1],'MarkerFaceColor',[0 1 1],'DisplayName',sprintf('LES: r = %3.2f; p = %3.2f',rRTles(1,2),pRTles(1,2)))
% 
% set(gca,'FontSize',20)
% xlabel('Theoretical RH')
% ylabel('RH')
% title('(a) \gamma')
% xlim([30 100])
% ylim([30 100])
% plot([30:10:100],[30:10:100],'k','HandleVisibility','off')
% legend('location','northwest')
% 
% %%%%%%%%%%%%%%%%
% %Panel 2: RH  vs. theoretical RH with vary depth
% %%%%%%%%%%%%%%%%
% RHtheory = RH_theory300_vary_depth;
% RHtheory_par = RHtheory(iPAR);
% RHtheory_crm = RHtheory(iCRM);
% RHtheory_ver = RHtheory(iVER);
% RHtheory_les = RHtheory(iLES);
% RHtheory_exp = RHtheory(iEXP);
% 
% %correlations
% [rRT,pRT] = corrcoef(RH(logical(~isnan(RH).*~isnan(RHtheory))),RHtheory(logical(~isnan(RH).*~isnan(RHtheory))))
% [rRTexp,pRTexp] = corrcoef(RHexp(logical(~isnan(RHexp).*~isnan(RHtheory_exp))),RHtheory_exp(logical(~isnan(RHexp).*~isnan(RHtheory_exp))))
% [rRTpar,pRTpar] = corrcoef(RHpar(logical(~isnan(RHpar).*~isnan(RHtheory_par))),RHtheory_par(logical(~isnan(RHpar).*~isnan(RHtheory_par))));
% [rRTcrm,pRTcrm] = corrcoef(RHcrm(logical(~isnan(RHcrm).*~isnan(RHtheory_crm))),RHtheory_crm(logical(~isnan(RHcrm).*~isnan(RHtheory_crm))));
% [rRTver,pRTver] = corrcoef(RHver(logical(~isnan(RHver).*~isnan(RHtheory_ver))),RHtheory_ver(logical(~isnan(RHver).*~isnan(RHtheory_ver))));
% [rRTles,pRTles] = corrcoef(RHles(logical(~isnan(RHles).*~isnan(RHtheory_les))),RHtheory_les(logical(~isnan(RHles).*~isnan(RHtheory_les))));
% %set NaNs and tiny to 0 correlation
% rRTpar(rRTpar<10^-4)=0;
% rRTpar(isnan(rRTpar))=0;
% pRTpar(isnan(pRTpar))=1;
% rRTcrm(rRTcrm<10^-4)=0;
% rRTcrm(isnan(rRTcrm))=0;
% pRTcrm(isnan(pRTcrm))=1;
% rRTver(rRTver<10^-4)=0;
% rRTver(isnan(rRTver))=0;
% pRTver(isnan(pRTver))=1;
% rRTles(rRTles<10^-4)=0;
% rRTles(isnan(rRTles))=0;
% pRTles(isnan(pRTles))=1;
% 
% subplot(2,2,2)
% plot(RHtheory_par*100,RHpar,'o','MarkerSize',10,'MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[1 0 0],'DisplayName',sprintf('PAR: r = %3.2f; p = %3.2f',rRTpar(1,2),pRTpar(1,2)))
% hold on
% plot(RHtheory_crm*100,RHcrm,'o','MarkerSize',10,'MarkerEdgeColor',[0 0 1],'MarkerFaceColor',[0 0 1],'DisplayName',sprintf('CRM: r = %3.2f; p = %3.2f',rRTcrm(1,2),pRTcrm(1,2)))
% plot(RHtheory_ver*100,RHver,'o','MarkerSize',10,'MarkerEdgeColor',[0 180/255 1],'MarkerFaceColor',[0 180/255 1],'DisplayName',sprintf('VER: r = %3.2f; p = %3.2f',rRTver(1,2),pRTver(1,2)))
% plot(RHtheory_les*100,RHles,'o','MarkerSize',10,'MarkerEdgeColor',[0 1 1],'MarkerFaceColor',[0 1 1],'DisplayName',sprintf('LES: r = %3.2f; p = %3.2f',rRTles(1,2),pRTles(1,2)))
% 
% set(gca,'FontSize',20)
% xlabel('Theoretical RH')
% ylabel('RH')
% title('(b) Temperature of Convecting Top')
% xlim([30 100])
% ylim([30 100])
% plot([30:10:100],[30:10:100],'k','HandleVisibility','off')
% legend('location','northwest')
% 
% %%%%%%%%%%%%%%%%
% %Panel 3: RH  vs. theoretical RH with vary PE
% %%%%%%%%%%%%%%%%
% RHtheory = RH_theory300_vary_PE;
% RHtheory_par = RHtheory(iPAR);
% RHtheory_crm = RHtheory(iCRM);
% RHtheory_ver = RHtheory(iVER);
% RHtheory_les = RHtheory(iLES);
% RHtheory_exp = RHtheory(iEXP);
% 
% %correlations
% [rRT,pRT] = corrcoef(RH(logical(~isnan(RH).*~isnan(RHtheory))),RHtheory(logical(~isnan(RH).*~isnan(RHtheory))))
% [rRTexp,pRTexp] = corrcoef(RHexp(logical(~isnan(RHexp).*~isnan(RHtheory_exp))),RHtheory_exp(logical(~isnan(RHexp).*~isnan(RHtheory_exp))))
% [rRTpar,pRTpar] = corrcoef(RHpar(logical(~isnan(RHpar).*~isnan(RHtheory_par))),RHtheory_par(logical(~isnan(RHpar).*~isnan(RHtheory_par))));
% [rRTcrm,pRTcrm] = corrcoef(RHcrm(logical(~isnan(RHcrm).*~isnan(RHtheory_crm))),RHtheory_crm(logical(~isnan(RHcrm).*~isnan(RHtheory_crm))));
% [rRTver,pRTver] = corrcoef(RHver(logical(~isnan(RHver).*~isnan(RHtheory_ver))),RHtheory_ver(logical(~isnan(RHver).*~isnan(RHtheory_ver))));
% [rRTles,pRTles] = corrcoef(RHles(logical(~isnan(RHles).*~isnan(RHtheory_les))),RHtheory_les(logical(~isnan(RHles).*~isnan(RHtheory_les))));
% 
% subplot(2,2,3)
% plot(RHtheory_par*100,RHpar,'o','MarkerSize',10,'MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[1 0 0],'DisplayName',sprintf('PAR: r = %3.2f; p = %3.2f',rRTpar(1,2),pRTpar(1,2)))
% hold on
% plot(RHtheory_crm*100,RHcrm,'o','MarkerSize',10,'MarkerEdgeColor',[0 0 1],'MarkerFaceColor',[0 0 1],'DisplayName',sprintf('CRM: r = %3.2f; p = %3.2f',rRTcrm(1,2),pRTcrm(1,2)))
% plot(RHtheory_ver*100,RHver,'o','MarkerSize',10,'MarkerEdgeColor',[0 180/255 1],'MarkerFaceColor',[0 180/255 1],'DisplayName',sprintf('VER: r = %3.2f; p = %3.2f',rRTver(1,2),pRTver(1,2)))
% plot(RHtheory_les*100,RHles,'o','MarkerSize',10,'MarkerEdgeColor',[0 1 1],'MarkerFaceColor',[0 1 1],'DisplayName',sprintf('LES: r = %3.2f; p = %3.2f',rRTles(1,2),pRTles(1,2)))
% 
% set(gca,'FontSize',20)
% xlabel('Theoretical RH')
% ylabel('RH')
% title('(c) Precip. Efficiency')
% xlim([30 100])
% ylim([30 100])
% plot([30:10:100],[30:10:100],'k','HandleVisibility','off')
% legend('location','northwest')
% 
% %%%%%%%%%%%%%%%%
% %Panel 4: RH  vs. theoretical RH with vary entrainment
% %%%%%%%%%%%%%%%%
% RHtheory = RH_theory300_vary_eps;
% RHtheory_par = RHtheory(iPAR);
% RHtheory_crm = RHtheory(iCRM);
% RHtheory_ver = RHtheory(iVER);
% RHtheory_les = RHtheory(iLES);
% RHtheory_exp = RHtheory(iEXP);
% 
% %correlations
% [rRT,pRT] = corrcoef(RH(logical(~isnan(RH).*~isnan(RHtheory))),RHtheory(logical(~isnan(RH).*~isnan(RHtheory))))
% [rRTexp,pRTexp] = corrcoef(RHexp(logical(~isnan(RHexp).*~isnan(RHtheory_exp))),RHtheory_exp(logical(~isnan(RHexp).*~isnan(RHtheory_exp))))
% [rRTpar,pRTpar] = corrcoef(RHpar(logical(~isnan(RHpar).*~isnan(RHtheory_par))),RHtheory_par(logical(~isnan(RHpar).*~isnan(RHtheory_par))));
% [rRTcrm,pRTcrm] = corrcoef(RHcrm(logical(~isnan(RHcrm).*~isnan(RHtheory_crm))),RHtheory_crm(logical(~isnan(RHcrm).*~isnan(RHtheory_crm))));
% [rRTver,pRTver] = corrcoef(RHver(logical(~isnan(RHver).*~isnan(RHtheory_ver))),RHtheory_ver(logical(~isnan(RHver).*~isnan(RHtheory_ver))));
% [rRTles,pRTles] = corrcoef(RHles(logical(~isnan(RHles).*~isnan(RHtheory_les))),RHtheory_les(logical(~isnan(RHles).*~isnan(RHtheory_les))));
% 
% subplot(2,2,4)
% plot(RHtheory_par*100,RHpar,'o','MarkerSize',10,'MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[1 0 0],'DisplayName',sprintf('PAR: r = %3.2f; p = %3.2f',rRTpar(1,2),pRTpar(1,2)))
% hold on
% plot(RHtheory_crm*100,RHcrm,'o','MarkerSize',10,'MarkerEdgeColor',[0 0 1],'MarkerFaceColor',[0 0 1],'DisplayName',sprintf('CRM: r = %3.2f; p = %3.2f',rRTcrm(1,2),pRTcrm(1,2)))
% plot(RHtheory_ver*100,RHver,'o','MarkerSize',10,'MarkerEdgeColor',[0 180/255 1],'MarkerFaceColor',[0 180/255 1],'DisplayName',sprintf('VER: r = %3.2f; p = %3.2f',rRTver(1,2),pRTver(1,2)))
% plot(RHtheory_les*100,RHles,'o','MarkerSize',10,'MarkerEdgeColor',[0 1 1],'MarkerFaceColor',[0 1 1],'DisplayName',sprintf('LES: r = %3.2f; p = %3.2f',rRTles(1,2),pRTles(1,2)))
% 
% set(gca,'FontSize',20)
% xlabel('Theoretical RH')
% ylabel('RH')
% title('(d) Entrainment')
% xlim([30 100])
% ylim([30 100])
% plot([30:10:100],[30:10:100],'k','HandleVisibility','off')
% legend('location','northwest')
% 
% % gcfsavepdf(['FigXX-RH-' num2str(z1) '-' num2str(z2) '.pdf'])
% 
% 
% %% Explanations for intermodel spread in RH change with warming
% %Force RH at 305 to be NaN for DALES and DALES-damping to ignore them for
% %changes with warming, because they are aggregated at 305
% hur305small_avg(strcmp(small_model_list,'DALES')==1)=NaN;
% hur305small_avg(strcmp(small_model_list,'DALES-damping')==1)=NaN;
% 
% %change in RH with warming
% dRH = (100/10)*(hur305small_avg-hur295small_avg)./hur300small_avg; % percent/K over 10K range
% dRHpar = dRH(iPAR);
% dRHcrm = dRH(iCRM);
% dRHver = dRH(iVER);
% dRHles = dRH(iLES);
% dRHexp = dRH(iEXP);
% 
% % figure('Position',[100 100 700 1600]) %for 4x1 subplot
% figure('Position',[100 100 1000 1000])
% 
% %%%%%%%%%%%%%%%%
% %Panel 1: RH change with warming vs. theoretical RH with vary gamma change with warming
% %%%%%%%%%%%%%%%%
% dRHtheory = (100/10)*(RH_theory305_vary_gamma-RH_theory295_vary_gamma)./RH_theory300_vary_gamma; % percent/K over 10K range;
% dRHtheory_par = dRHtheory(iPAR);
% dRHtheory_crm = dRHtheory(iCRM);
% dRHtheory_ver = dRHtheory(iVER);
% dRHtheory_les = dRHtheory(iLES);
% dRHtheory_exp = dRHtheory(iEXP);
% 
% %correlations
% [rdRT,pdRT] = corrcoef(dRH(logical(~isnan(dRH).*~isnan(dRHtheory))),dRHtheory(logical(~isnan(dRH).*~isnan(dRHtheory))))
% [rdRTexp,pdRTexp] = corrcoef(dRHexp(logical(~isnan(dRHexp).*~isnan(dRHtheory_exp))),dRHtheory_exp(logical(~isnan(dRHexp).*~isnan(dRHtheory_exp))))
% [rdRTpar,pdRTpar] = corrcoef(dRHpar(logical(~isnan(dRHpar).*~isnan(dRHtheory_par))),dRHtheory_par(logical(~isnan(dRHpar).*~isnan(dRHtheory_par))));
% [rdRTcrm,pdRTcrm] = corrcoef(dRHcrm(logical(~isnan(dRHcrm).*~isnan(dRHtheory_crm))),dRHtheory_crm(logical(~isnan(dRHcrm).*~isnan(dRHtheory_crm))));
% [rdRTver,pdRTver] = corrcoef(dRHver(logical(~isnan(dRHver).*~isnan(dRHtheory_ver))),dRHtheory_ver(logical(~isnan(dRHver).*~isnan(dRHtheory_ver))));
% [rdRTles,pdRTles] = corrcoef(dRHles(logical(~isnan(dRHles).*~isnan(dRHtheory_les))),dRHtheory_les(logical(~isnan(dRHles).*~isnan(dRHtheory_les))));
% 
% subplot(2,2,1)
% plot(dRHtheory_par,dRHpar,'o','MarkerSize',10,'MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[1 0 0],'DisplayName',sprintf('PAR: r = %3.2f; p = %3.2f',rdRTpar(1,2),pdRTpar(1,2)))
% hold on
% plot(dRHtheory_crm,dRHcrm,'o','MarkerSize',10,'MarkerEdgeColor',[0 0 1],'MarkerFaceColor',[0 0 1],'DisplayName',sprintf('CRM: r = %3.2f; p = %3.2f',rdRTcrm(1,2),pdRTcrm(1,2)))
% plot(dRHtheory_ver,dRHver,'o','MarkerSize',10,'MarkerEdgeColor',[0 180/255 1],'MarkerFaceColor',[0 180/255 1],'DisplayName',sprintf('VER: r = %3.2f; p = %3.2f',rdRTver(1,2),pdRTver(1,2)))
% plot(dRHtheory_les,dRHles,'o','MarkerSize',10,'MarkerEdgeColor',[0 1 1],'MarkerFaceColor',[0 1 1],'DisplayName',sprintf('LES: r = %3.2f; p = %3.2f',rdRTles(1,2),pdRTles(1,2)))
% 
% set(gca,'FontSize',20)
% xlabel('Change in Theoretical RH (% K^{-1})')
% ylabel('Change in RH (% K^{-1})')
% title('(a) \gamma')
% xlim([-1 3])
% ylim([-1 3])
% vline(0,'k--')
% hline(0,'k--')
% plot([-4:1:4],[-4:1:4],'k','HandleVisibility','off')
% legend('location','northwest')
% 
% %%%%%%%%%%%%%%%%
% %Panel 2: RH change with warming vs. theoretical RH with vary depth change with warming
% %%%%%%%%%%%%%%%%
% dRHtheory = (100/10)*(RH_theory305_vary_depth-RH_theory295_vary_depth)./RH_theory300_vary_depth; % percent/K over 10K range;
% dRHtheory_par = dRHtheory(iPAR);
% dRHtheory_crm = dRHtheory(iCRM);
% dRHtheory_ver = dRHtheory(iVER);
% dRHtheory_les = dRHtheory(iLES);
% dRHtheory_exp = dRHtheory(iEXP);
% 
% %correlations
% [rdRT,pdRT] = corrcoef(dRH(logical(~isnan(dRH).*~isnan(dRHtheory))),dRHtheory(logical(~isnan(dRH).*~isnan(dRHtheory))))
% [rdRTexp,pdRTexp] = corrcoef(dRHexp(logical(~isnan(dRHexp).*~isnan(dRHtheory_exp))),dRHtheory_exp(logical(~isnan(dRHexp).*~isnan(dRHtheory_exp))))
% [rdRTpar,pdRTpar] = corrcoef(dRHpar(logical(~isnan(dRHpar).*~isnan(dRHtheory_par))),dRHtheory_par(logical(~isnan(dRHpar).*~isnan(dRHtheory_par))));
% [rdRTcrm,pdRTcrm] = corrcoef(dRHcrm(logical(~isnan(dRHcrm).*~isnan(dRHtheory_crm))),dRHtheory_crm(logical(~isnan(dRHcrm).*~isnan(dRHtheory_crm))));
% [rdRTver,pdRTver] = corrcoef(dRHver(logical(~isnan(dRHver).*~isnan(dRHtheory_ver))),dRHtheory_ver(logical(~isnan(dRHver).*~isnan(dRHtheory_ver))));
% [rdRTles,pdRTles] = corrcoef(dRHles(logical(~isnan(dRHles).*~isnan(dRHtheory_les))),dRHtheory_les(logical(~isnan(dRHles).*~isnan(dRHtheory_les))));
% %set NaNs and tiny to 0 correlation
% rdRTpar(rdRTpar<10^-4)=0;
% rdRTpar(isnan(rdRTpar))=0;
% pdRTpar(isnan(pdRTpar))=1;
% rdRTcrm(rdRTcrm<10^-4)=0;
% rdRTcrm(isnan(rdRTcrm))=0;
% pdRTcrm(isnan(pdRTcrm))=1;
% rdRTver(rdRTver<10^-4)=0;
% rdRTver(isnan(rdRTver))=0;
% pdRTver(isnan(pdRTver))=1;
% rdRTles(rdRTles<10^-4)=0;
% rdRTles(isnan(rdRTles))=0;
% pdRTles(isnan(pdRTles))=1;
% 
% subplot(2,2,2)
% plot(dRHtheory_par,dRHpar,'o','MarkerSize',10,'MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[1 0 0],'DisplayName',sprintf('PAR: r = %3.2f; p = %3.2f',rdRTpar(1,2),pdRTpar(1,2)))
% hold on
% plot(dRHtheory_crm,dRHcrm,'o','MarkerSize',10,'MarkerEdgeColor',[0 0 1],'MarkerFaceColor',[0 0 1],'DisplayName',sprintf('CRM: r = %3.2f; p = %3.2f',rdRTcrm(1,2),pdRTcrm(1,2)))
% plot(dRHtheory_ver,dRHver,'o','MarkerSize',10,'MarkerEdgeColor',[0 180/255 1],'MarkerFaceColor',[0 180/255 1],'DisplayName',sprintf('VER: r = %3.2f; p = %3.2f',rdRTver(1,2),pdRTver(1,2)))
% plot(dRHtheory_les,dRHles,'o','MarkerSize',10,'MarkerEdgeColor',[0 1 1],'MarkerFaceColor',[0 1 1],'DisplayName',sprintf('LES: r = %3.2f; p = %3.2f',rdRTles(1,2),pdRTles(1,2)))
% 
% set(gca,'FontSize',20)
% xlabel('Change in Theoretical RH (% K^{-1})')
% ylabel('Change in RH (% K^{-1})')
% title('(b) Temperature of Convecting Top')
% xlim([-1 3])
% ylim([-1 3])
% vline(0,'k--')
% hline(0,'k--')
% plot([-4:1:4],[-4:1:4],'k','HandleVisibility','off')
% legend('location','northwest')
% 
% %%%%%%%%%%%%%%%%
% %Panel 3: RH change with warming vs. theoretical RH with vary PE change with warming
% %%%%%%%%%%%%%%%%
% dRHtheory = (100/10)*(RH_theory305_vary_PE-RH_theory295_vary_PE)./RH_theory300_vary_PE; % percent/K over 10K range;
% dRHtheory_par = dRHtheory(iPAR);
% dRHtheory_crm = dRHtheory(iCRM);
% dRHtheory_ver = dRHtheory(iVER);
% dRHtheory_les = dRHtheory(iLES);
% dRHtheory_exp = dRHtheory(iEXP);
% 
% %correlations
% [rdRT,pdRT] = corrcoef(dRH(logical(~isnan(dRH).*~isnan(dRHtheory))),dRHtheory(logical(~isnan(dRH).*~isnan(dRHtheory))))
% [rdRTexp,pdRTexp] = corrcoef(dRHexp(logical(~isnan(dRHexp).*~isnan(dRHtheory_exp))),dRHtheory_exp(logical(~isnan(dRHexp).*~isnan(dRHtheory_exp))))
% [rdRTpar,pdRTpar] = corrcoef(dRHpar(logical(~isnan(dRHpar).*~isnan(dRHtheory_par))),dRHtheory_par(logical(~isnan(dRHpar).*~isnan(dRHtheory_par))));
% [rdRTcrm,pdRTcrm] = corrcoef(dRHcrm(logical(~isnan(dRHcrm).*~isnan(dRHtheory_crm))),dRHtheory_crm(logical(~isnan(dRHcrm).*~isnan(dRHtheory_crm))));
% [rdRTver,pdRTver] = corrcoef(dRHver(logical(~isnan(dRHver).*~isnan(dRHtheory_ver))),dRHtheory_ver(logical(~isnan(dRHver).*~isnan(dRHtheory_ver))));
% [rdRTles,pdRTles] = corrcoef(dRHles(logical(~isnan(dRHles).*~isnan(dRHtheory_les))),dRHtheory_les(logical(~isnan(dRHles).*~isnan(dRHtheory_les))));
% 
% subplot(2,2,3)
% plot(dRHtheory_par,dRHpar,'o','MarkerSize',10,'MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[1 0 0],'DisplayName',sprintf('PAR: r = %3.2f; p = %3.2f',rdRTpar(1,2),pdRTpar(1,2)))
% hold on
% plot(dRHtheory_crm,dRHcrm,'o','MarkerSize',10,'MarkerEdgeColor',[0 0 1],'MarkerFaceColor',[0 0 1],'DisplayName',sprintf('CRM: r = %3.2f; p = %3.2f',rdRTcrm(1,2),pdRTcrm(1,2)))
% plot(dRHtheory_ver,dRHver,'o','MarkerSize',10,'MarkerEdgeColor',[0 180/255 1],'MarkerFaceColor',[0 180/255 1],'DisplayName',sprintf('VER: r = %3.2f; p = %3.2f',rdRTver(1,2),pdRTver(1,2)))
% plot(dRHtheory_les,dRHles,'o','MarkerSize',10,'MarkerEdgeColor',[0 1 1],'MarkerFaceColor',[0 1 1],'DisplayName',sprintf('LES: r = %3.2f; p = %3.2f',rdRTles(1,2),pdRTles(1,2)))
% 
% set(gca,'FontSize',20)
% xlabel('Change in Theoretical RH (% K^{-1})')
% ylabel('Change in RH (% K^{-1})')
% title('(c) Precipitation Efficiency')
% xlim([-1 3])
% ylim([-1 3])
% vline(0,'k--')
% hline(0,'k--')
% plot([-4:1:4],[-4:1:4],'k','HandleVisibility','off')
% legend('location','northwest')
% 
% %%%%%%%%%%%%%%%%
% %Panel 4: RH change with warming vs. theoretical RH with vary entrainment change with warming
% %%%%%%%%%%%%%%%%
% dRHtheory = (100/10)*(RH_theory305_vary_eps-RH_theory295_vary_eps)./RH_theory300_vary_eps; % percent/K over 10K range;
% dRHtheory_par = dRHtheory(iPAR);
% dRHtheory_crm = dRHtheory(iCRM);
% dRHtheory_ver = dRHtheory(iVER);
% dRHtheory_les = dRHtheory(iLES);
% dRHtheory_exp = dRHtheory(iEXP);
% 
% %correlations
% [rdRT,pdRT] = corrcoef(dRH(logical(~isnan(dRH).*~isnan(dRHtheory))),dRHtheory(logical(~isnan(dRH).*~isnan(dRHtheory))))
% [rdRTexp,pdRTexp] = corrcoef(dRHexp(logical(~isnan(dRHexp).*~isnan(dRHtheory_exp))),dRHtheory_exp(logical(~isnan(dRHexp).*~isnan(dRHtheory_exp))))
% [rdRTpar,pdRTpar] = corrcoef(dRHpar(logical(~isnan(dRHpar).*~isnan(dRHtheory_par))),dRHtheory_par(logical(~isnan(dRHpar).*~isnan(dRHtheory_par))));
% [rdRTcrm,pdRTcrm] = corrcoef(dRHcrm(logical(~isnan(dRHcrm).*~isnan(dRHtheory_crm))),dRHtheory_crm(logical(~isnan(dRHcrm).*~isnan(dRHtheory_crm))));
% [rdRTver,pdRTver] = corrcoef(dRHver(logical(~isnan(dRHver).*~isnan(dRHtheory_ver))),dRHtheory_ver(logical(~isnan(dRHver).*~isnan(dRHtheory_ver))));
% [rdRTles,pdRTles] = corrcoef(dRHles(logical(~isnan(dRHles).*~isnan(dRHtheory_les))),dRHtheory_les(logical(~isnan(dRHles).*~isnan(dRHtheory_les))));
% 
% subplot(2,2,4)
% plot(dRHtheory_par,dRHpar,'o','MarkerSize',10,'MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[1 0 0],'DisplayName',sprintf('PAR: r = %3.2f; p = %3.2f',rdRTpar(1,2),pdRTpar(1,2)))
% hold on
% plot(dRHtheory_crm,dRHcrm,'o','MarkerSize',10,'MarkerEdgeColor',[0 0 1],'MarkerFaceColor',[0 0 1],'DisplayName',sprintf('CRM: r = %3.2f; p = %3.2f',rdRTcrm(1,2),pdRTcrm(1,2)))
% plot(dRHtheory_ver,dRHver,'o','MarkerSize',10,'MarkerEdgeColor',[0 180/255 1],'MarkerFaceColor',[0 180/255 1],'DisplayName',sprintf('VER: r = %3.2f; p = %3.2f',rdRTver(1,2),pdRTver(1,2)))
% plot(dRHtheory_les,dRHles,'o','MarkerSize',10,'MarkerEdgeColor',[0 1 1],'MarkerFaceColor',[0 1 1],'DisplayName',sprintf('LES: r = %3.2f; p = %3.2f',rdRTles(1,2),pdRTles(1,2)))
% 
% set(gca,'FontSize',20)
% xlabel('Change in Theoretical RH (% K^{-1})')
% ylabel('Change in RH (% K^{-1})')
% title('(d) Entrainment')
% xlim([-1 3])
% ylim([-1 3])
% vline(0,'k--')
% hline(0,'k--')
% plot([-4:1:4],[-4:1:4],'k','HandleVisibility','off')
% legend('location','northwest')
% 
% 
% % gcfsavepdf(['FigXX-RHchange-' num2str(z1) '-' num2str(z2) '.pdf'])
% 
% 
% 
