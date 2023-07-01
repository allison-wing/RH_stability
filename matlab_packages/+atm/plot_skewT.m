function [Tx,Tdx,py] = plot_skewT(varargin)
%
% Plot a skew-T log-p diagram
% Plot data onto a skew-T log-p diagram
%

%% Defaults %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Example input data
p = [100 150 200 250 300 350 400 450 500 550 600 650 700 750 800 825 850 875 900    925 950  975 1000].*100;
T = [-28 -31 -30 -37 -32 -25 -20 -10  0   3   7   10  15  17  14  16  20  22  24    26   28   30   33];
Td =[NaN NaN NaN NaN NaN NaN -40 -43 -20 -10  -5 -12  1   5   11  15  14  13  13.5  14   14.5 15  16];



% Plot markers as well as lines
plot_markers = 0;
bw = 0;

%% Read in inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin > 0;    T = varargin{1}; end
if nargin > 1;    Td = varargin{2}; end
if nargin > 2;    p = varargin{3}; end
if nargin > 3;    bw = varargin{4}; end
if nargin > 4;    plot_markers = varargin{5}; end    


%% Plotting Options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% axis limits
temp_lim = [-30 50];                    % Limit to temperature axis in deg C
pres_lim = [100 1050];                  % Limit to pressure axis in hPa


% Lines to plot
pres_ticks = 100:100:1000;              % Pressure tick marks in hPa

temp_ticks = -100:10:60;                % Temperature tick marks in deg C
thta_ticks = 200:10:600;                % Potential temperature tik marks in deg C
thte_ticks = 200:20:600;                % Equivalent potential temperature tick marks in deg C

rsat_ticks = [0.4 1 2 3 5 8 12 20 30 40 60];  % Mixing ratio ticks in g/kg
rsat_pmin  = 500;                       % Maximum pressure to plot mixing ratio lines (hPa)



% Colours
if bw
    col_Tp = [0 0 0.5];
    col_Td = [0.3 0.3 1];

    col_T = [0 0 0];
    col_th= [0.5 0.5 0.5];
    col_the= [0.35 0.35 0.35];
    col_r = [0.1 0.1 0.1];
    col_p = [0.25 0.25 0.25];
else
    col_Tp = [0 0 0];
    col_Td = [0.4 0.4 0.4];

    col_T = [1 0 0];
    col_th= [0.6 0.4 0.2];
    col_the= [0.2 0.6 0.2];
    col_r = [0.2 0.4 0.6];
    col_p = [0.25 0.25 0.25];
end


%% Setup the axis and some vectors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load the constants for thermodynamic calculations: Bolton without ice 
c = atm.load_constants('bolton',0,40);

fig.set_default('paper')

fig.bfig(20,20)
set(gcf,'defaultaxesfontsize',10)

axbl = axes('position',[0.1300    0.1100    0.7750    0.8150]);

% Limits to axes in axes coordinates
x_lim = temp_lim;           
y_lim = log(pres_lim.*100);   


% Make a pressure vector (Pa)
pres = (1100:-10:10).*100;

% Make a Temperature vector (K)
temp = 183.15:1:363.15;

% Axes
X = temp-273.15;   % Temperature in deg C 
Y = log(pres);     % Log of pressure 


%% Calculate the positions of the profiles to plot

py = log(p);
Tx = T + ( x_lim(2)-x_lim(1) )./( y_lim(1)-y_lim(2) ).* (py-y_lim(2));
Tdx = Td + ( x_lim(2)-x_lim(1) )./( y_lim(1)-y_lim(2) ).* (py-y_lim(2));




%% Populate the matrices of thermodynamic variables %%%%%%%%%%%%%%%%%%%%%%%
temp_m = zeros(length(pres),length(temp));
pres_m = zeros(length(pres),length(temp));

% Calculate temperature and pressure
for i = 1:length(pres)
    pres_m(i,:)  = pres(i);
    temp_m(i,:)  = temp - ( x_lim(2)-x_lim(1) )./( y_lim(1)-y_lim(2) ).* (Y(i)-y_lim(2));
end
   
% Calculate potential temperatures and saturation mixing ratio
thta_m  = temp_m.*(c.p00./pres_m).^(c.Rd./c.cp);
rsat_m  = atm.r_sat(temp_m,pres_m,'bolton',0,40);
thte_m = atm.calculate_theta_ep(temp_m,rsat_m,pres_m);



%% Convert the matrices to the units we want %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pres_m = pres_m./100;       % convert to hPa
temp_m = temp_m - 273.15;   % convert to deg C

rsat_m = rsat_m.*1000;      % convert to g/kg


%% Contour the matrices of thermodynamic variables %%%%%%%%%%%%%%%%%%%%%%%%

% Pressure
contour(X,Y,pres_m,pres_ticks,'color',col_p);

hold on

% Temperature
[cT,hT] = contour(X,Y,temp_m,temp_ticks,'color',col_T);
%clabel(cT,hT,'color','r','LabelSpacing',1008)

% Potential temperature
[cTh,hTh] =  contour(X,Y,thta_m,thta_ticks,'color',col_th);
clabel(cTh,hTh,[260 280 300 320 340 360 380 400],'color',col_th,'LabelSpacing',500)

% Saturation mixing ratio
Ir = find(pres <= rsat_pmin.*100,1,'first');
[cr,hr] =  contour(X,Y(1:Ir),rsat_m(1:Ir,:),rsat_ticks,'--','color',col_r,'labelspacing',1000);
clabel(cr,hr,'color',col_r)

% Equivalent potential temperature
[cthe,hthe] = contour(X,Y,thte_m,thte_ticks,'color',col_the);
clabel(cthe,hthe,[260 300 340 380 420],'color',col_the,'LabelSpacing',400)


%% Set axis limits etc %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(gca,'xlim',x_lim,'ylim',y_lim,'ydir','reverse')

for i = 1:length(pres_ticks)
   Y_labels{i} = num2str(pres_ticks(i));
end

set(gca,'xtick',temp_ticks,'ytick',log(pres_ticks.*100),'yticklabel',Y_labels);
set(gca,'xcolor',col_T)

ylabel('pressure (hPa)')
xlabel('temperature (\circ{}C)')
box off


%% Plot the profile on top %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes(axbl)

plot(Tx,py,'k-','color',col_Tp,'linewidth',2,'markersize',4,'markerfacecolor',col_Tp)
hold on
plot(Tdx,py,'-','color',col_Td,'linewidth',2,'markersize',4,'markerfacecolor',col_Td)

if plot_markers
  plot(Tx,py,'ko','color',col_Tp,'linewidth',2,'markersize',4,'markerfacecolor',col_Tp)
  plot(Tdx,py,'o','color',col_Td,'linewidth',2,'markersize',4,'markerfacecolor',col_Td)
end


%% Secondary axis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

axtr = axes('position',[0.1300    0.1100    0.7750    0.8150])
set(gca, 'Color', 'none')
set(gca,'yaxislocation','right','xaxislocation','top')

% Make right hand axis 
rh_ticks = temp_ticks(temp_ticks >temp_lim(1) & temp_ticks<temp_lim(2));

rh_tick_label = {};
for i = 1:length(rh_ticks)
    rh_tick_label{i} = num2str(rh_ticks(i));
end

rh_ticks = 1 - (rh_ticks - temp_lim(1))./(temp_lim(2)-temp_lim(1));

rh_ticks = flip(rh_ticks);
rh_tick_label = flip(rh_tick_label);

set(gca,'ycolor',col_T,'ytick',rh_ticks,'yticklabel',rh_tick_label)

% Make top axis
% Axis liits, starts with highest of top
top_lim = [x_lim(1)- (x_lim(2)-x_lim(1)) x_lim(1)];



set(gca,'xlim',top_lim,'xtick',temp_ticks);

set(gca,'xcolor',col_T)
top_tick_label = get(gca,'xticklabel');
%top_tick_label{1} = '';
%set(gca,'xticklabel',top_tick_label);
