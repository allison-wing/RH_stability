function [PE_imp,eps_imp] = find_PE_epsilon(PE_R16,epsilon_R16,CAPE_R16,RH_R16,capeval,hurval)
% [PE_imp,eps_imp] = find_PE_epsilon(PE_R16,epsilon_R16,CAPE_R16,RH_R16,capeval,hurval)
% Function to find value of precipitation efficiency and entrainment implied
% by Romps (2016) (R16) theory using simulated value of CAPE and RH
% 
% Input:
% PE_R16: array of precipitation efficiency values the theory is defined on
% epsilon_R16: array of epsilon values the theory is defined on
% CAPE_R16: matrix of CAPE values predicted by R16 given PE and epsilon (PE x epsilon)
% RH_R16: matrix of RH values predicted by R16 given PE and epsilon (PE x epsilon)
% capeval: CAPE value from simulation
% hurval: RH value from simulation
% 
% Output:
% PE_imp: theory-implied value of precipitation efficiency
% eps_imp: theory-implied value of entrainment (epsilon)


% Method:
% Calculate contours of CAPE and RH for the simulated values and
% see where they intersect.

% The contours
C1 = contourc(epsilon_R16,PE_R16,CAPE_R16,[capeval capeval]);
C2 = contourc(epsilon_R16,PE_R16,RH_R16,[hurval hurval]);

% Find intersection
[eps_imp,PE_imp] = polyxpoly(C1(1,2:end),C1(2,2:end),C2(1,2:end),C2(2,2:end));

% %find contour as function of PE and entrainment for specific value of CAPE
% %and hur
% figure;
% [Mcape]=contour(PE_R16,epsilon_R16,CAPE_R16',[capeval capeval]);
% hold on
% [Mhur]=contour(PE_R16,epsilon_R16,RH_R16',[hurval hurval]);
% 
% %ignore first column of M
% %first row is x points - PE
% %second row is y points - epsilon
% %each column is a coordinate for the contour of that value
% 
% %plot curves against each other
% % figure;
% % plot(Mcape(1,2:end),Mcape(2,2:end))
% % hold on
% % plot(Mhur(1,2:end),Mhur(2,2:end))
% xlabel('Precipitation Efficiency (PE)')
% ylabel('Entrainment (\epsilon)')
% 
% %find where the curves intersect. This is the theory-implied value of PE
% %and epsilon given the CAPE and hur values
% [PE_imp,eps_imp] = intersections(Mcape(1,2:end),Mcape(2,2:end),Mhur(1,2:end),Mhur(2,2:end),true);
% 
% % %what to do if multiple values are returned
% % %remove any that are outside our theoretical range
% % PE_imp(PE_imp<0)=NaN;
% % PE_imp(PE_imp>1)=NaN;
% % eps_imp(eps_imp<0)=NaN;
% % eps_imp(eps_imp>3e-3)=NaN;
% % 
% % %remove any that are outside the range for that contour
% % PE_imp(PE_imp>max(Mcape(1,:)))=NaN;
% % PE_imp(PE_imp<min(Mcape(1,:)))=NaN;
% % PE_imp(PE_imp>max(Mhur(1,:)))=NaN;
% % PE_imp(PE_imp<min(Mhur(1,:)))=NaN;
% % 
% % eps_imp(eps_imp>max(Mcape(2,:)))=NaN;
% % eps_imp(eps_imp<min(Mcape(2,:)))=NaN;
% % eps_imp(eps_imp>max(Mhur(2,:)))=NaN;
% % eps_imp(eps_imp<min(Mhur(2,:)))=NaN;
% 
% hold on
% plot(PE_imp,eps_imp,'k*')


return
end
