function [PE_vec,ent_vec,RH_mat,CAPE_mat,gammab] = fn_plot_CAPE_and_RH_RCE(doplot,Tb,pb,Tt,gammaLCL)

%
% Calculate CAPE as a function of RH for the zero-buoyancy plume model
%
% Performs calculation following analytic theory of Romps (2016; R16)
%
%Input:
% doplot = 1 to make the figure, 0 to not plot
% Tb = Cloud base temperature (K)
% pb = Cloud base pressure (pa)
% Tt = Temperature of top of convecting layer (K)
%   gammaLCL = gamma at the LCL (m^-1)


%% Input parameters

if nargin < 1
    %defaults
    doplot = 0;
    Tb = 295;                   % Cloud base temperature (K)
    pb = 95000;                 % Cloud base pressure (pa)
    Tt = 220;                   % Temperature of top of convecting layer (K)% FAT assumption
elseif nargin<2
    %defaults
    Tb = 295;                   % Cloud base temperature (K)
    pb = 95000;                 % Cloud base pressure (pa)
    Tt = 220;                   % Temperature of top of convecting layer (K)% FAT assumption
elseif nargin<3
    %defaults
    pb = 95000;                 % Cloud base pressure (pa)
    Tt = 220;                   % Temperature of top of convecting layer (K)% FAT assumption
elseif nargin<4
    %defaults
    Tt = 220;                   % Temperature of top of convecting layer (K)% FAT assumption
end

% Now calculate the phase space
ent_vec = [0:0.01:10].*1e-3;
PE_vec = [0:0.01:1];
[ent_mat,PE_mat] = meshgrid(ent_vec,PE_vec);

%% Calculate CAPE according to the analytic theory of Romps (2016)
[CAPE_mat,RH_mat,CAPE_simple_mat,gammab] = calculate_CAPE_theory(Tb,Tt,pb,ent_mat,PE_mat,gammaLCL,'gamma');


%% Plot figure showing CAPE vs RH for R16 solution

if doplot==1
    % Which values should we plot
    PE_plot = [0.2 0.4 0.6 0.8 1];
    eps_plot = [0.05 0.1 0.2 0.4 0.8].*1e-3;
    
    % Which values do we label
    PE_lab = [0.2 0.6  1];
    PE_pos = [2500 3500 2000];
    eps_lab = [0.05 0.1 0.2 0.4 0.8].*1e-3;
    eps_pos = [0.7 0.7 0.7 0.7 0.9];
    
    % Make the figure
    %fig.bfig(16,12)
    % figure
    
    set(gca,'ylim',[0 6000])
    set(gca,'xlim',[0 1])
    hold on
    
    % Plot the lines of constant precipitation efficiency
    for i = 1:length(PE_plot)
        
        % Find the precipitation efficiency we want
        I = find(PE_vec==PE_plot(i));
        
        % Plot the line
        pp = plot(RH_mat(I,:),CAPE_mat(I,:),'k','LineWidth',2,'HandleVisibility','off');
        
        % Label the line
        if ismember(PE_plot(i),PE_lab)
            j = find(PE_lab==PE_plot(i));
            ll = fig.inline_label(pp,['PE = ' num2str(PE_plot(i))],[],PE_pos(j));
        end
        
    end
    
    
    % Plot the lines of constant entrainment rate
    for i = 1:length(eps_plot)
        
        % Find the precipitation efficiency we want
        I = find(ent_vec==eps_plot(i));
        
        % Plot the line
        pp = plot(RH_mat(:,I),CAPE_mat(:,I),'color',[0.75 0.75 0.75],'LineWidth',2,'HandleVisibility','off');
        
        % Label the line
        if ismember(eps_plot(i),eps_lab)
            j = find(eps_lab==eps_plot(i));
            ll = fig.inline_label(pp,['\epsilon = ' num2str(eps_plot(i).*1000)],eps_pos(j),[],'color',[0.5 0.5 0.5]);
        end
        
    end
    
    set(gca,'FontSize',16)
    xlabel('Relative Humidity')
    ylabel('CAPE(J kg^{-1})')
    box off
end

return
end




