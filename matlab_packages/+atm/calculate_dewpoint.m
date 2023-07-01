function Td = calculate_dewpoint(T,p,r,varargin)
%
% Calculate the dew point temperature Td such that
%
% e_sat(Td) = e
%
% where e is the vapour pressure.
%
% Td = calculate_dewpoint(T,p,r,[type,ice,deltaT]) 
%
% For Bolton thermodynamics with no ice, there is a closed form solution,
% otherwise an iteration is used with the Bolton solution as the first
% guess
%
% In principle this function can deal with a mixed phase, but it is
% probably better to use the no ice case.
%

c = atm.load_constants(varargin{1:end});

% Calculate the vapour pressure
e = r.*p./(c.eps + r);


% Use the Bolton approximation as a first guess
% The liquid-ice breakdown is not exact
Tdl = 243.5.*log(e./611.2)./(17.67 - log(e./611.2)) +273.15;
Tdi = 265.49.*log(e./611.2)./(21.8745584 - log(e./611.2)) +273.15;
[fliq,fice] = atm.calculate_frac_ice(Tdl,varargin{1:end});
Td = Tdl.*fliq + Tdi.*fice;


% Find the error in this approximation
e_test = atm.e_sat(Td,varargin{1:end});
e_err_max = max(abs(e_test - e));

% Iterate to find the actual dew point temperature
i_iter = 0;
while e_err_max > 1e-2
    i_iter = i_iter+1;
 
    desdT = atm.desatdT(T,varargin{1:end});
    Td = Td - (e_test - e)./desdT;
    
    e_test = atm.e_sat(Td,varargin{1:end});
    
    e_err_max = max(abs(e_test - e));
    
    if i_iter > 20 && i_iter <= 50
       % Now we are really likely to be not converged. Print some output to warn the user
       disp(['iter: ' num2str(i_iter) ', error = ' num2str(e_err_max) ' K'])
    elseif i_iter > 50
       % Pull the plug
       error('vapour pressure inversion did not converge. This may be because the input data is corrupt.')
    end
    
end