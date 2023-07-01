function [Tw,Tw_simple] = calculate_isobaric_wetbulb(T,p,rt,varargin)
%
% Function to calculate the isobaric wetbulb temperatre
%
% [Tw,Tw_simple] = calculate_isobaric_wetbulb(T,p,rt,[type,ice,deltaT])
%
% Calculates wetbulb temperature (Tw) based on Temperature (T) pressure (p)
% and total water mixing ratio (rt). If parcel is supersaturated, Tw = T.
%
% Two different methods to calculate the wetbulb temperature
%
% 1) simple: the traditional method. Saturate the parcel by adding vapour
% at constant enthalpy.
%
% 2) default: method that takes into account the fact that enthalpy
% changes when you add condensed water to a parcel. The physics of this
% method is as follows:
%
%      STEP 1: Add the required amount of condensed water to the parcel 
%              isothermally and isobarically
%
%      STEP 2: Evaporate this water until the parcel is saturated
%
%
% The challenge is that we do not know the amount of water to add at the
% beginning, so we need to iterate to find this.
%
% In testing, the two methods do not differ by more than 0.25 K, and
% typically the difference is much smaller. As long as one is consistent,
% it does not particularly matter which method is used.
%

%% Calculate thermodynamics of parcel

% Load constants
thermo_opts = {varargin{1:end}};
c = atm.load_constants(thermo_opts{:});

% Ensure the parcel is not supersaturated
[rvp,~,~] = atm.saturation_adjustment(p,T,rt,thermo_opts{:});

% calculate enthalpy per unit dry air mass:
h = atm.calculate_enthalpy(T,p,rvp,thermo_opts{1:end});

%% Calculate the wetbulb using the simple method

% Now calculate the temperature assuming fixed enthalpy and saturation
% Use the input temperature as a first guess.
[Tw,rv,~,~] = atm.invert_enthalpy(h,p,'sat',T);

Tw_simple = Tw;


%% Calulate the wetbulb using the comprehensive method

% This requires an iteration until the change in temperature is reduced
% below a threshold
max_dT = 1;
i_iter = 0;
while max_dT > 0.0001
    i_iter = i_iter+1;

   %% STEP 1: Add water to the parcel isothermally and isobarically

   % Use the previous iteration to estimate the amount of water needed
   rc = rv - rvp;

   % Calculate the phase(s) of the condensed water based on the parcel temp
   [fliq,fice] = atm.calculate_frac_ice(T,thermo_opts{1:end});

   % Recalculate the enthalphy.
   h = atm.calculate_enthalpy(T,p,rv,fliq.*rc,fice.*rc,thermo_opts{1:end});

   %% STEP 2: Evaporate the water at constant enthalpy
   % Now calculate the temperature assuming fixed enthalpy and saturation
   % Use the last iteration wetbulb temperature as a first guess.
   [Tw_new,rv,rl,ri] = atm.invert_enthalpy(h,p,'sat',Tw);
   max_dT = max(abs(Tw_new - Tw));
   Tw = Tw_new;
   
   if i_iter > 20 && i_iter <= 50
       % Now we are really likely to be not converged. Print some output to warn the user
       disp(['iter: ' num2str(i_iter) ', error = ' num2str(max_dT) ' K'])
   elseif i_iter > 50
       % Pull the plug
       error('vapour pressure inversion did not converge. This may be because the input data is corrupt.')
   end
   
   

end
    

    
    
