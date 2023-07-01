function const = thermo_constants(varargin)
%
% Load a table of constants for thermodynamic calculations
% Includes physical constants, default microphysical parameters
%
% Optional input argument chooses the set of constants to use
%
%
% 'default' : Exact thermodynamics under the assumption that isobaric specific heat capacities are constant
% 'bolton'  : Use Bolton's formulas for saturation vapor pressure
% 'teten'   : Use Teten's formulas for saturation vapor pressure
% 'sam'     : Use constants consistent with SAM cloud-resolving model
% 'fms'     : Use constants consistent with Simplified GCM based on the GFDL FMS model 
%		(as used by O'Gorman & Schneider etc.)
%

    % Assume default constants unless input variable
    const.type = 'default';
    if nargin >= 1; const.type = varargin{1}; end
    


    %% Default microphysical parameters

    const.gamma = 0;            % Default precipitation fallout parameter (0 is reversible, 1 is pseudo-adiabatic)
    const.ice = 1;              % include ice?
    const.deltaT = 40;		% mixed-phase range (K)
    

    %% Physical constants (based on constants in CM1)

    % Isobaric heat capacities (J/kg/K)
    const.cp        = 1005.7;    	% dry air
    const.cpv       = 1870.0;		% water vapor
    const.cpl       = 4190.0;		% liquid water
    const.cpi       = 2106.0;		% solid water

    % Gas constants (J/kg/K)
    const.Rd        = 287.04;		% dry air
    const.Rv        = 461.5;

    % Latent heats (J/kg)
    const.Lv0       = 2501000.0;	% liquid-gas
    const.Ls0       = 2834000.0;	% solid-gas

    % gravitational acceleration (m/s^2)
    const.g         = 9.81;

    % Liquid water density (kg/m^3)
    const.rhol = 1000;


    %% Reference values

    const.T0        = 273.16; 		% temperature (K)
    const.p00 	    = 100000;		% Pressure (Pa)
    const.e0        = 611.2;		% vapor pressure (Pa)


    %% Non thermodynamic constants
    const.sigma = 5.67e-8;

    %% Altered constants for different thermodynamics
    switch lower(const.type)

       case 'default'

       case 'bolton'

       case 'teten'

       case 'fms'

         % choose the following for consistency with the FMS code
          const.Rd                  = 287.04;                 % J / kg / K
          const.eps                 = 0.622;                  % gas constant water vapor [J/kg/K]
          const.Rv                  = const.Rd/const.eps;

          const.kappa             = 2.0/7.0;
          const.cp                = const.Rd/const.kappa;     % specific heat at constant pressure
          const.cpv               = const.cp;                 % specific heat water vapor [J/kg/K] - SAME AS FOR DRY AIR
          const.cpl               = const.cp;                 % specific heat liquid water [J/kg/K] - SAME AS FOR DRY AIR
          const.cp_ocean          = 3989.24495292815;         % heat capacity ocean water [J/kg/K]
          const.rhol              = 1000;                     % density of liquid water [kg/m^3]
          const.Lv0               = 2.5e6;                    % latent heat of evaporation [J/kg]
          const.ice         	  = 0;
          const.g                 = 9.80665; 

        case 'fmsi' 
          
                   % choose the following for consistency with the FMS code
          const.Rd                  = 287.04;                 % J / kg / K
          const.eps                 = 0.622;                  % gas constant water vapor [J/kg/K]
          const.Rv                  = const.Rd/const.eps;

          const.kappa             = 2.0/7.0;
          const.cp                = const.Rd/const.kappa;     % specific heat at constant pressure
          const.cpv               = const.cp;                 % specific heat water vapor [J/kg/K] - SAME AS FOR DRY AIR
          const.cpl               = const.cp;                 % specific heat liquid water [J/kg/K] - SAME AS FOR DRY AIR
          const.cp_ocean          = 3989.24495292815;         % heat capacity ocean water [J/kg/K]
          const.rhol              = 1000;                     % density of liquid water [kg/m^3]
          const.Lv0               = 2.5e6;                    % latent heat of evaporation [J/kg]
          const.Ls0               = 2.83e6;                    % latent heat of sublimation [J/kg]
          const.ice         	  = 1;
          const.g                 = 9.80665; 

          
          
       case 'sam'

          const.cp = 1004.0;      %                 ! Specific heat of air, J/kg/K
          const.Lv0 = 2.5104e+06; %                 ! Latent heat of condensation, J/kg
          const.Ls0 = 2.8440e+06; %                 ! Latent heat of sublimation, J/kg

          const.Rv = 461.;        %                 ! Gas constant for water vapor, J/kg/K
          const.Rd = 287.;        %                 ! Gas constant for dry air, J/kg/K

          const.deltaT = 20;
          const.tbgmin = 253.16;  %                 ! Minimum temperature for cloud water., K
          const.tbgmax = 273.16;  %                 ! Maximum temperature for cloud ice, K
          const.tprmin = 268.16;  %                 ! Minimum temperature for rain, K
          const.tprmax = 283.16;  %                 ! Maximum temperature for snow+graupel, K
          const.tgrmin = 223.16;  %                 ! Minimum temperature for snow, K
          const.tgrmax = 283.16;  %                 ! Maximum temperature for graupel, K

          const.deltaT = 20;			    % Consistent with tbgmin/tbgmax above

       otherwise
          error('Cannot find the constants requested')

    end

    %% Numerical parameters
    const.epsT = 1e-3;                        % Precision of temperature iterations (K)
    
    %% Override constants if input
    if nargin >= 2; const.ice  = varargin{2}; end 
    if nargin >= 3; const.deltaT  = varargin{3}; end 
    
    %% Derived parameters

    const.cv        = const.cp-const.Rd;
    const.cvv       = const.cpv-const.Rv;
    const.eps       = const.Rd/const.Rv;
    
