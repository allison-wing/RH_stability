function [T,rv,rl,ri,T_rho,varargout] = calculate_adiabat(Ts,rts,ps,p_in,varargin)
%
% Function to calculate an adiabatic parcel ascent with various thermodynamic assumptions
%
% [T,rv,rl,ri,T_rho,[T_LCL,p_LCL,p_out]] = calculate_adiabat(Ts,rts,ps,p,[gamma,type,ice,deltaT,N_factor])
%
% INPUTS:
%
% Ts = parcel temperature (K)
% rts = parcel total mixing ratio (kg/kg)
% ps = Parcel pressure (Pa)
% p = pressure levels to integrate parcel (Pa)
%
% gamma = fraction of condensate removed from parcel
% type, ice, deltaT are microphysical parameters
% See e_sat.m for details.
%
% OUTPUTS:
%
% T = temperature of parcel at all pressure levels (K)
% rv = water vapor mixing raio of parcel (kg/kg)
% rl = liquid water mixing ratio of parcel (kg/kg)
% ri = solid water mixing ratio of parcel (kg/kg)
% T_rho = density temperature of parcel (K)
%
% p must either be a vector or of size [Np size(Ts)] and should be sorted largest to smallest
% Ts, ps, and rts may be up to 3 dimensions
% The outputs are of size [Np size(Ts)]
%
% Parcel is lifted from ps to the top of the pressure matrix
% all values of p greater than  ps are filled with Nans
%
% This is helpful to vectorize the parcel ascent calculation over a 
% matrix of different columns which may have different surface pressures
%
%
% The adiabat is calculated in three steps.
%
% 1) The LCL pressure is calculated using the method of Romps (2017).
%    This method provides an analytic exact solution for the LCL in terms
%    of the Lambert W function. The MATLAB library for the W function must
%    be loaded, so this part of the calculation can be slow the first time
%    the adiabat calculator is called. 
%
%    If the LCL lies in the mixed-phase range, the calculation is no longer
%    analytic, and an iteration is required. See calculate_LCL.m for more
%    details
%   
% 2) For levels below the LCL, the parcel is assumed to conserve its
%    potential temperature theta = T*(p/p_LCL)^(Rm/cpm) and its mixing
%    ratio r. This gives the temperature by simple rearrangement
%
% 3) Above the LCL, the parcel temperature is calculated by by integrating 
%    the equation for the lapse rate:
%
%                       dT/dp = f(T,rt,p) 
%
%    See the function calculate_dTdp_adiabatic.m for more details
%
%    The ascent is performed by calculating dT/dp assuming no fallout and
%    then removing a fraction gamma of the water condensed during that step
%    isothermally. In the case of gamma == 1, dT/dp is calculated assuming 
%    no condensed water for pseudoadiabatic ascent.
%
% 4) For the case of truly reversible ascent (gamma = 0, deltaT = 0). An
%    alternative method is used in which entropy is assumed to be conserved,
%    and at each level an inversion is performed to calculate the
%    temperature. Special care is taken to correctly calculate the
%    isothermal freezing region in which the condensate freezes into ice. 
%    The inversion of entropy requires an iteration, and so can be slow, 
%    but it is very accurate. 
%
% The accuracy of these calculations may be checked using the script
% test_script.m. For the case of no ice and the default thermodynamics,
% an exact solution for a reversible parcel ascent is conservation of
% entropy. This script is accurate to within 0.5 K for dp = 50 hPa, and to
% within 0.1 K for dp = 20 hPa.
% 

%% Optional arguments

% Read in the thermodyanmic arguments

% irreversibilty, thermodynamics, and treatment of ice:
% These modify the default characteristics set in load_constants.m
c = atm.load_constants(varargin{2:end});
gamma = c.gamma;

if nargin >= 5; gamma  = varargin{1};  end


% Read in N_factor determining the vertical grid

% N_factor:
% N_factor = 1   % use the input pressure array 
% N_factor > 1   % Increase the number of pressure levels by N_factor
% N_factor < 0   % Make the mean pressure interval smaller than N_factor (Pa)

% Default setting is to ensure ~5 hPa intervals unless ascent is fully
% reversible (in which case no integration is performed and dp should not
% influence accuracy)
if gamma == 0 && c.deltaT == 0
    N_factor = 1;
else
    N_factor = -500;
end

if nargin >= 9; N_factor  = varargin{1};  end


%% Set up the grids

% Arrangement of the columns
xygrid = size(Ts);

% Resmaple the pressure matrix to ensure the integration can be performed accurately. 
[p,Ip,dp] = resample_pressure_matrix(p_in,xygrid,N_factor);

if c.gamma > 0 || c.deltaT > 0
    persistent big_dp
    if max(abs(dp(:))) >= 5000 & isempty(big_dp)
         warning(['calculate_adiabat.m: You have pressure intervals ' ...
                  'greater than or equal to 50 hPa. ' ...
                  'It is recommened that you subsample the pressure ' ...
                  'matrix for better accuracy of the integration.'])
         big_dp = 1;
    end
end

% Initialize profiles
T  = nan([size(p,1) xygrid ]); % temperature
rv = nan([size(p,1) xygrid ]); % vapor
rl = nan([size(p,1) xygrid ]); % liquid
ri = nan([size(p,1) xygrid ]); % solid
rt = nan([size(p,1) xygrid ]); % total water


% Initialize surface values
[rvs,~,~] = atm.saturation_adjustment(ps,Ts,rts,c.type,c.ice,c.deltaT);
Rm = c.Rd + c.Rv.*rvs;
cpm = c.cp + c.cpv.*rvs;


if gamma ==1
    rts = rvs; 
end
   
% Calculate LCL level
[T_LCL,p_LCL] = atm.calculate_LCL(Ts,rts,ps,varargin{2:end});

% If LCL is below surface, set it to surface
T_LCL(p_LCL>ps) = Ts(p_LCL>ps);
p_LCL(p_LCL>ps) = ps(p_LCL>ps);


  
if c.deltaT > 0 || gamma > 0

    %% non reversible case    
    
    % Initialise diagnostics
    I_above_LCL_prev = zeros(size(Ts));
    
    % Loop from lowest level upwards (from highest pressure to lowest)
    for k = 1:size(p,1)

        % Pressure at this level
        pk = reshape(p(k,:,:,:),xygrid);

        % Find regions where the surface pressure is larger than the current level p(k)
        I_above_surface = ps>=pk;

        % Find regions where the LCL pressure is smaller than the current level p(k)
        I_subcloud = I_above_surface & p_LCL<=reshape(p(k,:,:,:),xygrid);

        % In subcloud layer, total water is conserved
        rt(k,I_subcloud) = rts(I_subcloud);
        
        % Invert the entropy in the subcloud region
        [T(k,I_subcloud),rv(k,I_subcloud),rl(k,I_subcloud),ri(k,I_subcloud)] ...
            = invert_subcloud(rts(I_subcloud),pk(I_subcloud),T_LCL(I_subcloud),p_LCL(I_subcloud),c);

        % Find regions where the LCL pressure is larger than the current level p(k)
        I_above_LCL = I_above_surface & ~I_subcloud;

        % Find regions where the LCL is between levels (k-1) and k
        I_at_LCL = I_above_LCL & ~I_above_LCL_prev;

        %% Integrate upwards (towards lower pressure) from the LCL to the current model level for the regions where the current level p(k) is just above the LCL

        % Calculate the dp
        if sum(I_at_LCL(:))>0
            dps = reshape(p(k,:,:,:),xygrid) - p_LCL;
            [T(k,I_at_LCL),rt(k,I_at_LCL),rv(k,I_at_LCL),rl(k,I_at_LCL),ri(k,I_at_LCL)] = integrate_upwards(T_LCL(I_at_LCL),rts(I_at_LCL),p_LCL(I_at_LCL),dps(I_at_LCL),gamma,c.type,c.ice,c.deltaT);
        end


        %% Now integrate from level k to level k+1
        if k < size(p,1)
            [T(k+1,I_above_LCL),rt(k+1,I_above_LCL),rv(k+1,I_above_LCL),rl(k+1,I_above_LCL),ri(k+1,I_above_LCL)] = integrate_upwards(T(k,I_above_LCL),rt(k,I_above_LCL),p(k,I_above_LCL),dp(k,I_above_LCL),gamma,c.type,c.ice,c.deltaT);
        end

        % Set the previous levels diagnostics
        I_above_LCL_prev = I_above_LCL;
        
    
    end

else
    
    %% For truly reversible ascent

    
    % Calculate entropy of parcel. This is conserved for reversible ascent
    % The rest of this section is deadicated to inverting the entropy for
    % temperature at each pressur level
    s_parcel = atm.calculate_entropy(Ts,ps,rts,c.type,c.ice,c.deltaT);

    % Initialise entropies at freezing
    s_frz_min = nan(size(Ts));
    s_frz_max = nan(size(Ts));
    
    
    for k = 1:size(p,1)

        % Pressure at this level
        pk = reshape(p(k,:,:,:),xygrid);
        
        % Find regions where the surface pressure is larger than the current level p(k)
        I_above_surface = ps>=pk;
        
        % total water is conserved
        rt(k,I_above_surface) = rts(I_above_surface);
        
        % total water at this level
        rtk = reshape(rt(k,:,:,:),xygrid);
                
        % Find regions where the LCL pressure is smaller than the current level p(k)
        I_subcloud = I_above_surface & p_LCL<=pk;

        % Invert the entropy in the subcloud region
        [T(k,I_subcloud),rv(k,I_subcloud),rl(k,I_subcloud),ri(k,I_subcloud)] ...
            = invert_subcloud(rts(I_subcloud),pk(I_subcloud),T_LCL(I_subcloud),p_LCL(I_subcloud),c);

        % Find regions where the LCL pressure is larger than the current level p(k)
        I_above_LCL = I_above_surface & ~I_subcloud;
    
        if sum(I_above_LCL > 0)

            
            % Calculate maximum and minimum entropy for air parcel at freezing
            s_frz_max(I_above_LCL) = atm.calculate_entropy(c.T0.*(1+2.*eps),p(k,I_above_LCL),rt(k,I_above_LCL),c.type,c.ice,c.deltaT);
            s_frz_min(I_above_LCL) = atm.calculate_entropy(c.T0.*(1-2.*eps),p(k,I_above_LCL),rt(k,I_above_LCL),c.type,c.ice,c.deltaT);
       

            %% Solve for regions within the freezing zone

            % Find regions within the freezing zone
            I_frz = s_parcel > s_frz_min & s_parcel < s_frz_max;

            [T(k,I_frz),rv(k,I_frz),rl(k,I_frz),ri(k,I_frz)] = invert_freezing(s_parcel(I_frz),rtk(I_frz),pk(I_frz),s_frz_min(I_frz),s_frz_max(I_frz),c);
            
           
            %% Solve for regions outside the freezing zone

            % Find regions not in freezing zone
            I_nfrz = I_above_LCL & ~I_frz;

            % Construct two estimates of the temperature at level k
            
            % First estimate: assume temperature is the same as level below
            % This is always greater than the actual T(k) value
            if k == 1
               T_pos = T_LCL(I_nfrz);
               p_lower = p_LCL(I_nfrz);
            else
               T_pos = shiftdim(T(k-1,I_nfrz),1);
               p_lower = shiftdim(p(k-1,I_nfrz),1);
            end
            
            % Second estimate: calculate temperature assuming parcel conserves (approximate) potential temperature
            % This is always lower than the actual T(k) value
            T_neg = T_pos.*(shiftdim(p(k,I_nfrz),1)./p_lower).^(c.Rd./c.cp);
            
            [T(k,I_nfrz),rv(k,I_nfrz),rl(k,I_nfrz),ri(k,I_nfrz)] = invert_saturated(s_parcel(I_nfrz),rtk(I_nfrz),pk(I_nfrz),T_neg,T_pos,c);
            
        end
    end
    
end


% Resample back to the original pressures
T = T(Ip); T = reshape(T,[size(p_in,1) xygrid]);
rv = rv(Ip); rv = reshape(rv,[size(p_in,1) xygrid]);
rl = rl(Ip); rl = reshape(rl,[size(p_in,1) xygrid]);
ri = ri(Ip); ri = reshape(ri,[size(p_in,1) xygrid]);



T_rho = T.*(1 + rv./c.eps)./(1+rv+rl+ri);

if nargout > 5; varargout{1} = T_LCL; end
if nargout > 6; varargout{2} = p_LCL; end
if nargout > 7; varargout{3} = p; end

end


function [Tkp1,rtkp1,rvkp1,rlkp1,rikp1] = integrate_upwards(Tk,rtk,pk,dp,gamma,type,ice,deltaT)
%
% Integrate the equation for enthalpy upwards. 
% 1) 2nd order Runge-Kutta while conserving total water
% 2) precipitation fall out by fraction gamma
%
% If gamma == 1, then the calculation of dT/dp is done assuming no water,
% and no precipitation fallout is required.
%
    % calculate predictor
    [dTdp,rvk] = atm.calculate_dTdp_adiabatic(Tk,rtk,pk,gamma,type,ice,deltaT);
    
    T1 = Tk + dTdp.*dp./2;
    p1 = (pk+dp./2);
    [rv1,rl1,ri1] = atm.saturation_adjustment(p1,T1,rtk,type,ice,deltaT);
    if rl1+ri1 > 0
       [rl1,ri1] = atm.simple_fallout(T1,rvk-rv1,rl1,ri1,gamma,type,ice,deltaT);
    end
    rt1 = rv1+rl1+ri1;
    
    % calculate corrector 
    dTdp = atm.calculate_dTdp_adiabatic(T1,rt1,p1,gamma,type,ice,deltaT);
    
    Tkp1 = Tk + dTdp.*dp;
    pkp1 = pk+dp;
    [rvkp1,rlkp1,rikp1] = atm.saturation_adjustment(pkp1,Tkp1,rtk,type,ice,deltaT);
    if rlkp1+rikp1 > 0
       [rlkp1,rikp1] = atm.simple_fallout(Tkp1,rvk-rvkp1,rlkp1,rikp1,gamma,type,ice,deltaT);
    end
    rtkp1 = rvkp1+rlkp1+rikp1;
    
    
end
    


function [T,rv,rl,ri] = invert_subcloud(rt,p,T_LCL,p_LCL,c)
% Invert the subcloud region
        
     Rm = c.Rd + c.Rv.*rt;
     cpm = c.cp + c.cpv.*rt;


     % Find the temperature in the regions below the LCL
     T = T_LCL.*(p./p_LCL).^(Rm./cpm);

     rv = rt;
     rl = zeros(size(rv));
     ri = zeros(size(rv));



end

function [T,rv,rl,ri] = invert_freezing(s,rt,p,s_frz_min,s_frz_max,c)
% Invert the freezing region
        
    % Set temperature and mixing ratio at freezing
    T = c.T0.*ones(size(p));
    rv = atm.r_sat(T,p,c.type,c.ice,c.deltaT);

    % Solve for the condensate breakdown in this region
    fliq = ( s - s_frz_min )./( s_frz_max- s_frz_min );


    rl = fliq     .* (rt - rv);
    ri = (1-fliq) .* (rt - rv);


end

function [T,rv,rl,ri] = invert_saturated(s,rt,p,T_neg,T_pos,c)
% Invert the saturated region   
    

    % Use Bisection to solve for the temperature
    % This is robust, but somewhat slow. Might be useful to coed up a
    % better method (e.g., Brent's method)
    %
    % This method is vectorised, but it iterates until every individual
    % calculation has converged. This is a little inefficient, as ideally
    % we should do different numbers of iterations for each calculation.
    % Given MATLABs efficiency at array multiplication, this is probably
    % only a minor consideration.
    %
    
    iter = 0;
    while max(abs(T_pos - T_neg)) > c.epsT

       % Calculate error in entropy for the initial range
       ds_pos = atm.calculate_entropy(T_pos,p,rt,c.type,c.ice,c.deltaT) - s;
       ds_neg = atm.calculate_entropy(T_neg,p,rt,c.type,c.ice,c.deltaT) - s;

       if any(ds_pos.*ds_neg > 0 )
          disp(['ds+ = ' num2str(ds_pos) ', ds- = ' num2str(ds_neg)])
          error(['Entropy inversion function: zero not contained in interval given. Inversion failed on iter: ' num2str(iter)])
       end
        
        
       T_bisect = (T_pos+T_neg)./2;    % Bisection method
       %T_bisect = ( T_neg.*ds_pos - T_pos.*ds_neg ) ./ ( ds_pos - ds_neg ); % False position method
       
       % Calculate error in new estimate
       ds_bi = atm.calculate_entropy(T_bisect,p,rt,c.type,c.ice,c.deltaT) - s;
       
       % Calculate new interval
       if ds_neg ==0
           T_pos = T_neg;
       end
       I = ds_bi.*ds_neg < 0;

       T_pos(I) = T_bisect(I);
       T_neg(~I) = T_bisect(~I);

       


        % Update counter
        iter = iter + 1;
        
        if iter > 30 
            disp(num2str(max(abs(T_pos - T_neg))))
            error('Entropy inversion function: method not converging')
        end
    end

    T = (T_pos+T_neg)./2;
    [rv,rl,ri] = atm.saturation_adjustment(p,T,rt,c.type,c.ice,c.deltaT);



end


function [p,Ip,dp] = resample_pressure_matrix(p_in,xygrid,N_in)
% Create appropriate pressure matrix
% This routine ensures the pressure matrix is the right size
% Additionally, depending on the varible N_factor, it resamples the
% pressure matrix to ensure that the parcel ascent calculated accurately
%
% [p,Ip,dp] = resample_pressure_matrix(p_in,xygrid,N_factor)
%
% Inputs:
%           p_in:   input pressure matrix, either np-vector or of size [np xygrid] (Pa)
%           xygrid: input size of parcel array 1, 2 or 3D (int)
%           N_in: factor by which to decrease pressure intervals (see below)
% 
% Outputs: 
%           p: new pressure matrix of size [np_new xygrid]
%           Ip: indices of old pressures in new pressure matrix
%           dp: pressure intervals
%
%
% If N_in > 0, we simply increase the size of the pressure vector by a
% factor N_in uniformly. This means that the pressure intervals are reduced
% by the factor N_in everywhere.
%
% If N_in < 0, it represents the target pressure interval. We increase the
% number of levels such that the average value of the pressure intervals at
% each level equal -N_in (in Pa)


% Ensure first dimension of pressure is the vertical
if sum(size(p_in)>1) == 1 	% If the pressure is a vector
  p_in = p_in(:);
end


if N_in < 0
    % If N_in < 0 calculate the N values required
    
    % Calculate average pressure interval of the input
    dp = diff(p_in,1,1);
    dp_mean = mean(dp(:,:),2);

    N = ceil(dp_mean./N_in);   
  
else
    % If N_in > 0, simply take input value
    N = repmat(N_in,[1 size(p_in,1)-1]);
end

if N_in > 1 || N_in < 0
    % Resample pressure levels
    
    % Initialize the higher-resolution arrays
    p = zeros([sum(N)+1 size(p_in,2) size(p_in,3) size(p_in,4)]);
    Ip = false([sum(N)+1 size(p_in,2) size(p_in,3) size(p_in,4)]);

    % Now loop over the levels
    k_out = 1;
    for k_in = 1:size(p_in,1)-1
        
        % For each level, set the next N values in the high-resolution arrays        
        Ip(k_out,:,:,:) = true;
        for j = 1:N(k_in)
           p(k_out+j-1,:,:,:) = p_in(k_in,:,:,:)+dp(k_in,:,:,:).*(j-1)./N(k_in);
        end

        k_out = k_out+N(k_in);
    end
    
    % Set the final value
    Ip(end,:,:,:) = true;
    p(end,:,:,:) = p_in(end,:,:,:);
    

else
    
    % Keep pressure the same
    Ip = true([sum(N)+1 size(p_in,2) size(p_in,3) size(p_in,4)]);
    p = p_in;
    
end
  
% Ensure pressure is the right size
if sum(size(p)>1) == 1 	
  p = repmat(p,[1 xygrid]);
end

% Calculate pressure intervals
dp = diff(p,1,1);




end


