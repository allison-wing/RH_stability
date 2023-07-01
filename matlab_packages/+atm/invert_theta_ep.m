function [T,r] = invert_theta_ep(theta_ep,r,Tstar,p)

% Invert the pseudo-equivelant potetial temperature
% using iteraction procedure. Should converge within a few iterations
% Use Bolton formula, taken from Emanuel (1994), p 132

c = atm.load_constants('bolton');
c.ice = 0;

% No condensation: simply calculate T assuming Tstar and r are unchanged.
T = theta_ep./ (  (c.p00./p).^(0.2854.*(1-0.28.*r)).*exp( r .* (1+0.81.*r) .* (3376./Tstar -2.54) )   );  

% Find out where the super saturation is
rsat = atm.r_sat(T,p,c.type,c.ice);
saturated = rsat<r;
Tsat = T(saturated);
psat = p(saturated);
rsat = rsat(saturated);
Tstarsat = Tstar(saturated);

pisat = (c.p00./psat).^(0.2854);

% Begin iteration to fix supersaturation
i_iter = 0;
dT = 1;
while dT > 0.01 % Assume we are converged if error in theta_ep < 0.01 K.

% Iteration number
i_iter = i_iter+1;


% Calculate error in theta_ep
[theta_ep_shoot,Tstarsat] = atm.calculate_theta_ep(Tsat,rsat,psat);
dtheta_ep = theta_ep(saturated) - theta_ep_shoot;

% Calculate partial derivatives (simply)
dthepdT = (  pisat.*exp( rsat .* (3376./Tstarsat -2.54) )   );  
dthepdr = 3376./Tstarsat.*Tsat.*(  pisat.*exp( rsat .*(3376./Tstarsat -2.54) )   );  
drdT = rsat.*c.Lv0./(c.Rv.*Tstarsat.^2);

dTdthep = 1./(dthepdT + dthepdr.*drdT);

% Find new Temperature in saturated regions
Tnew = Tsat + 0.8.*dtheta_ep.*dTdthep;


% Calculate difference
dT = max(abs(dtheta_ep(:)));

if i_iter > 12 
    Tnew(abs(dtheta_ep)>100) = nan;
end
Tsat = Tnew;
rsat = atm.r_sat(Tsat,psat,c.type,c.ice);

if i_iter > 12; disp(['iter: ' num2str(i_iter) ', error: ' num2str(dT,4) ' K']); end

end

T(saturated) = Tsat;
r(saturated) = rsat;

