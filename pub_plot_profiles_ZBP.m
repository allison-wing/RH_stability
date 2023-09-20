%% Plots model profiles (T deviation from moist adiabat) and compares them to ZBP model
% lift parcel from lowest model level
clear all

load('RH_stab_2-5.mat')

%Concatenate SSTs together (SST, models)
prof = cat(1,prof295small,prof300small,prof305small);
eps_imp = cat(1,eps_imp295small,eps_imp300small,eps_imp305small);
PE_imp = cat(1,PE_imp295small,PE_imp300small,PE_imp305small);
T_LCL = cat(1,T_LCL295,T_LCL300,T_LCL305);
p_LCL = cat(1,p_LCL295,p_LCL300,p_LCL305);
T_LNB = cat(1,T_LNB295,T_LNB300,T_LNB305);
p_LNB = cat(1,p_LNB295,p_LNB300,p_LNB305);

c = atm.load_constants;
gamma = 0; %reversible adiabat
deltaT = 40;
ice = 1;

SSTs = [295 300 305];
lw = [1 2 3];

for i = 1:length(small_model_list) %which model
    % for i = 1:3
    
    
    if strcmp(small_model_list(i),'DALES-VER')==1 || strcmp(small_model_list(i),'DALES-LES')==1 || ...
            strcmp(small_model_list(i),'DALES-damping-VER')==1 || ...
            strcmp(small_model_list(i) ,'UKMO-RA1-T')==1
        fprintf('skip model')
    else
        figure
        for k = 1:length(SSTs) %which SST
            
            %Calculate moist adiabat, lifting from lowest model levels
            [T_par,r_par,rl_par,ri_par,prof(k,i).tmarho] = atm.calculate_adiabat(prof(k,i).ta(1),prof(k,i).rv(1),prof(k,i).p(1)*100,prof(k,i).p*100);
                  
            %Calculate theoretical buoyancy profile
            % Find the pressures we want to solve for
            Itrop = prof(k,i).p*100<p_LCL(k,i) & prof(k,i).p*100>p_LNB(k,i);
            iLCL = find(prof(k,i).p*100>p_LCL(k,i),1,'last');
            iLNB = find(prof(k,i).p*100>p_LNB(k,i),1,'last');
            p_ZBP = [p_LCL(k,i); prof(k,i).p(Itrop)*100; p_LNB(k,i)];
            
            % Now do a slightly weird calculation to find z
            z_LCL = prof(k,i).z(iLCL)*1000 + c.Rd./c.g.*prof(k,i).tv(iLCL).*(log(prof(k,i).p(iLCL)*100)-log(p_LCL(k,i)));
            z_LNB = prof(k,i).z(iLNB)*1000 + c.Rd./c.g.*prof(k,i).tv(iLNB).*(log(prof(k,i).p(iLNB)*100)-log(p_LNB(k,i)));
            z_ZBP = [z_LCL; prof(k,i).z(Itrop)*1000; z_LNB];
            
            % Calculate ZBP model with a given set of pressures
            [T_ZBP,p_ZBP,z_ZBP,RH_ZBP,Tm_ZBP] = calculate_ZBP(z_ZBP,p_ZBP,T_LCL(k,i),eps_imp(k,i),PE_imp(k,i),'gamma');
                     
            %%%%%%%%%%%%%%%%
            % Plot %
            %%%%%%%%%%%%%%%%
            %Plot model profiles
            hold on
            plot(prof(k,i).tmarho-prof(k,i).tv,prof(k,i).z,'Color',prof295small(i).color,'LineWidth',lw(k),'DisplayName',small_model_list{i})
            
            %Plot ZBP profiles
            plot(Tm_ZBP-T_ZBP,z_ZBP/1000,'k','LineWidth',lw(k),'DisplayName',small_model_list{i})
            
        end %SST loop
        
        %set plot properties
        set(gca,'FontSize',18)
        %vline(0,'k--')
        xlabel('T_{m,\rho} - T_\rho (K)')
        ylabel('Z (km)')
        title(small_model_list{i})
        ylim([0 18])
        xlim([-5 15])
        %         hline(2,'k--')
        %         hline(5,'k--')
                %gcfsavepdf(['./proftest/Fig_prof_' small_model_list{i} '.pdf'])
    end %exclude models with missing data
end %model loop