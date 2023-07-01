function [y_out,t_out] = rk_step(fun,y_in,t_in,dt,varargin)
%
% Performs a single step of a Runge Kutta integration
% to solve an initial value problem. 
%
% The MATLAB builtin functions ODE45 and ODE23 and their relatives
% generally provide similar functionality. However, they are generally not 
% designed to perform the calculation one step at a time, which is useful
% in situations where the equation coefficients themselves are defined on a
% grid.
%
% y_out = rk_step(fun,t_in,y_in,dt [,method])
%
%
% Steps forward the equation
%
% dy/dt = fun(t,y)
%
% from the initial conditoin (t_in,y_in) by the step dt. 
% Outputs are the value t_out = t_in + dt and y_out, an estimate of y at t = t_out.
%
% Note that y can be a vector or ND array.
%
% INPUTS
% 
% fun       function handle representing dy/dt
% y_in      value of y at t = t_in
% t_in      initial time
% dt        timestep
%
% method    [optional]
%	    either a string determining the RK method to use
%           availale methods are
%		rk4      - standard 4th order Runge-Kutta (default)
%		euler    - simple Euler first order method
%		midpoint - 2nd order using midpoint rule
%
%	    otherwise a structure containing the vectors a,b and 
%	    matrix c that define an RK method.
%	    see the Wikipedia page on RK methods for more details
%
% OUTPUTS
%
% y_out    estimate of y_out at t = t_out
% t_out    t_in + dt
%
%
%
%


%% Set the integration method %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get the method we require
method = 'rk4';
if nargin > 4; method = lower(varargin{1}); end


% Set the matrices for the method
% This is not very ideal as we are recalculating the matrices at each step
% If you really need speed, might want to code up a bespoke solution
switch method

   case {'rk1','euler'}
      c = 0;
      b = 1;
      a = 0;

   case {'rk2','midpoint'}
      c = [0 0.5];
      b = [0 1];
      a = [0 0; 0.5 0];
      
   case {'rk3','ralston'}
      c = [0 0.5 3/4];
      b = [2 3 4]./9;
      a = [0 0 0; 0.5 0 0; 0 3/4 0];
       
   case {'rk4','classic','default'}     
      c = [0 0.5 0.5 1];
      b = [1 2 2 1]./6;
      a = [0 0 0 0; 0.5 0 0 0; 0 0.5 0 0; 0 0 1 0];

    case {'rk5','dormand'}
      % From Dormand & Prince (1980): A family of embedded Runge-Kutta formulae
      % This uses 7 levels and is designed to use an adaptive RK method.
      % Same algorithm as ODE45 matlab function.
      c = [    0           1/5       3/10        4/5         8/9         1       1   ];
      b = [  5179/57600	    0	    7571/16695 393/640	-92097/339200 187/2100  1/40 ];
      a = [    0            0         0           0           0          0       0; ...
              1/5           0         0           0           0          0       0; ...
              3/40         9/40       0           0           0          0       0; ...
             44/45       -56/15	    32/9          0           0          0       0; ...
            19372/6561 -25360/2187 64448/6561 -212/729        0          0       0; ...
            9017/3168  -355/33  46732/5247	    49/176   -5103/18656     0       0
             35/384	      0        500/1113	   125/192   -2187/6784	   11/84     0];
              
 
    otherwise

      if ischar(method)
          error('Unknown Runge-Kutta method. Perhaps you need to define your own matrices')
      elseif isstruct(method)
          c = method.c;
          b = method.b;
          a = method.a;
      else
          error('Malformed Runge-Kutta method.')
      end

end

% Number of stages in RK step
N = length(c);


%% Support for multidimensional functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Transpose the RK matrices
a = a';
b = b';
c = c';

% Set some dimension variables
ydims = size(y_in);
Ndims = length(ydims);
rem_first_dim = [2:Ndims+1 1];
Nsub = repmat({':'}, [1 Ndims]);

%% Integrate one RK step %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k = zeros([N ydims]);

% Calculate the intemerdiate steps
for i = 1:N
    
    % Calculate increment and ensure it is the same size as y_in 
    increm = permute(sum(a(1:i-1,i).*k(1:i-1,Nsub{:}),1),rem_first_dim);
    
    % Evaluate function at the intermediate point
    k(i,:) = fun( t_in + c(i).*dt , y_in + dt.*increm );
end

% Calculate the variables at the end of the timestep
y_out = y_in + permute(dt.*sum(k.*b,1),rem_first_dim);
t_out = t_in + dt;



end      

