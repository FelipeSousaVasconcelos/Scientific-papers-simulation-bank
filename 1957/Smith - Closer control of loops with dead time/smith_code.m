% @article{smith1957closer,
%   title={Closer control of loops with dead time},
%   author={Smith, Otto JM},
%   journal={Chemical engineering progress},
%   volume={53},
%   pages={217--219},
%   year={1957}
% }

% Felipe José de Sousa Vasconcelos
% First modification: 22/06/2024
% Last modification: 23/06/2024

%% Initialization
clear all
close all
clc

warning off

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot,'defaultTextInterpreter','latex');
rng('default');

%% simulation paramters
Tsim = 300;
Tref = 10;
ref = 1;
Tdist = 150;
dist = -0.1;

s = tf('s');
fig = 1;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% II. Process Transient
%
% Identification methods
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Process and curves:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K_real = 1.5;
T1_real = 0.2;
T2_real = 4;
L_real = 5;
G = zpk((K_real/((s+T1_real)*(s+T2_real)))*exp(-L_real*s));

% Get the step response
% [y, t] = step(G);
t = (0:0.01:35)';
u = ones(length(t),1);
y = lsim(G,u,t);

% Calculate the numerical derivative of the step response
dy = diff(y) ./ diff(t);

% Find the point where the slope is maximum
[~, idx_max_slope] = max(dy);

% The point of maximum slope
t_max_slope = t(idx_max_slope);
y_max_slope = y(idx_max_slope);

% Equation of the tangent line: y = slope*(t - t_tang) + y_tang
slope = dy(idx_max_slope);
tangent = slope * (t - t_max_slope) + y_max_slope;  % Tangent line

% Points to draw the tangent line
y1_tang = 0;
x1_tang = (y1_tang - y_max_slope)/slope + t_max_slope;
y2_tang = max(y);
x2_tang = (y2_tang - y_max_slope)/slope + t_max_slope;
tangent_lineX = linspace(x1_tang,x2_tang,length(t));
tangent_lineY = linspace(y1_tang,y2_tang,length(t));

% Points to draw the b line
y1_b = 0;
x1_b = (y1_b - y_max_slope)/slope + t_max_slope - tangent_lineX(1) + L_real;
y2_b = max(y);
x2_b = (y2_b - y_max_slope)/slope + t_max_slope - tangent_lineX(1) + L_real;
b_lineX = linspace(x1_b,x2_b,length(t));
b_lineY = linspace(y1_b,y2_b,length(t));

% Plot the step response and the tangent line
offset1 = max(y)*0.1;
hf = figure(fig); fig = fig + 1;
plot(t, y, 'b', 'LineWidth', 1.5); hold on;
plot(tangent_lineX, tangent_lineY, 'r--', 'LineWidth', 1.5);
plot(b_lineX, b_lineY, 'g--', 'LineWidth', 1.5);
plot(t_max_slope, y_max_slope, 'ko', 'MarkerFaceColor', 'k'); % Point of maximum slope
plot(tangent_lineX(end), tangent_lineY(end), 'ko', 'MarkerFaceColor', 'k');

line([0, 0], [min(y)-offset1, 0], 'Color', 'black', 'LineStyle', '--');
line([L_real, L_real], [min(y)-3*offset1, 0], 'Color', 'black', 'LineStyle', '--');
line([tangent_lineX(1), tangent_lineX(1)], [tangent_lineY(1),tangent_lineY(1)-4*offset1], 'Color', 'black', 'LineStyle', '--');
line([t_max_slope, t_max_slope], [min(y)-3*offset1, y_max_slope], 'Color', 'black', 'LineStyle', '--');
line([tangent_lineX(end), tangent_lineX(end)], [tangent_lineY(1)-4*offset1, tangent_lineY(end)], 'Color', 'black', 'LineStyle', '--');

% T_c = T_1 + T_2
ha = annotation('doublearrow');
ha.Parent = hf.CurrentAxes;  % associate annotation with current axes
ha.X = [t_max_slope, tangent_lineX(end)];
ha.Y = [tangent_lineY(1)-offset1, tangent_lineY(1)-offset1];
text(ha.X(1),ha.Y(1)+ha.Y(1)*0.5,'$T_c = T_1 + T_2$','interpreter','Latex');

% T_d
plot(0, 0, 'ko', 'MarkerFaceColor', 'k');
plot(L_real, 0, 'ko', 'MarkerFaceColor', 'k');
ha = annotation('doublearrow');
ha.Parent = hf.CurrentAxes;  % associate annotation with current axes
ha.X = [0, L_real];
ha.Y = [min(y)-offset1, min(y)-offset1];
text((ha.X(1)+ha.X(end))/2,ha.Y(1)+ha.Y(1)*0.5,'$T_d$','interpreter','Latex');

% T_b
ha = annotation('doublearrow');
ha.Parent = hf.CurrentAxes;  % associate annotation with current axes
ha.X = [L_real,tangent_lineX(1)];
ha.Y = [min(y)-3*offset1, min(y)-3*offset1];
text((ha.X(1)+ha.X(end))/2,ha.Y(1)+ha.Y(1)*0.2,'$T_b$','interpreter','Latex');

% T_f
ha = annotation('doublearrow');
ha.Parent = hf.CurrentAxes;  % associate annotation with current axes
ha.X = [tangent_lineX(1),t_max_slope];
ha.Y = [min(y)-3*offset1, min(y)-3*offset1];
text((ha.X(1)+ha.X(end))/2,ha.Y(1)+ha.Y(1)*0.2,'$T_f$','interpreter','Latex');

% T_a
ha = annotation('doublearrow');
ha.Parent = hf.CurrentAxes;  % associate annotation with current axes
ha.X = [tangent_lineX(1),tangent_lineX(end)];
ha.Y = [min(y)-4*offset1, min(y)-4*offset1];
text((ha.X(1)+ha.X(end))/2,ha.Y(1)+ha.Y(1)*0.2,'$T_a$','interpreter','Latex');

xlabel('Time (s)','interpreter','Latex');
ylabel('Output','interpreter','Latex');
title('Output transient curve or recorded controller','interpreter','Latex');
% text(t_max_slope, y_max_slope, 'b')
grid on;
axis([t(1) t(end) min(y)-5.5*offset1 y(end)])

Ta = tangent_lineX(end) - tangent_lineX(1);
Tb = tangent_lineX(1) - L_real;
Tc = tangent_lineX(end) - t_max_slope;
Td = L_real;
Tf = t_max_slope - tangent_lineX(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Method from Rudolf Carl Oldenbourg, Hans Sartorius
% The Dynamics of Automatic Controls - 1948
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fun = @(x) (1 + x) * x^(exp(1) / (1 - exp(1))) - (Tc / Ta); % Eq.(4)
x0 = 0.1;
options = optimoptions('fsolve', 'Display', 'iter'); % Display iteration info
x_sol = fsolve(fun, x0, options);

T1_ = Tc/(1+x_sol); % Eq.(2)
T1_Oldenbourg = Ta/((1/x_sol)^(exp(1) / (1 - exp(1)))); % Eq.(3)

T2_Oldenbourg = x_sol*T1_Oldenbourg;

% Function with delay
f_Oldenbourg_delayed = @(t) (t >= L_real) .* (1 - exp(-(t-L_real)./T1_Oldenbourg)/(1-x_sol) + (x_sol*exp(-(t-L_real)./T2_Oldenbourg))/(1-x_sol));

% Evaluate the delayed function
f_Oldenbourg = evalfr(G,0)* f_Oldenbourg_delayed(t);

G_Oldenbourg = zpk((1/((s+T1_Oldenbourg)*(s+T2_Oldenbourg)))*exp(-L_real*s));
y_Oldenbourg = lsim(G_Oldenbourg,u,t);
K_Oldenbourg = y(end)/y_Oldenbourg(end);
G_Oldenbourg = zpk((K_Oldenbourg/((s+T1_Oldenbourg)*(s+T2_Oldenbourg)))*exp(-L_real*s));
y_Oldenbourg = lsim(G_Oldenbourg,u,t);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Method from Smith
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculating a
% first one is y at Td+Tb
a1_index = find(t>=Td+Tb,1,'first');
a1 = y(a1_index)
g1 = a1*exp(1)
% second one is b at Td+Tb, which gives a*exp(1)
a2_index = find(b_lineX>=Td+Tb,1,'first');
g2 = b_lineY(a2_index) 
a2 = g2/exp(1)


a = a2;
g = a*(exp(1) + 0.53/(1+(150*a)^(-exp(1)))) % I am not sure what is written in Eq. (6)

Tb1 = g*Ta
Tb

if a >= 0.005
    T2 = Tb*(1+10*a+(exp(1)-1)*(30*a)^2)
else
    T2 = (Tb+Tf)*(1-200*(0.032-a)*(1+(0.086+0.0015/(0.032-a)^(-1)))^(-1))
end
T1 = Tc-T2

G_smith = zpk((1/((s+T1)*(s+T2)))*exp(-L_real*s));
y_smith = lsim(G_smith,u,t);
K_smith = y(end)/y_smith(end);
G_smith = zpk((K_smith/((s+T1)*(s+T2)))*exp(-L_real*s));
y_smith = lsim(G_smith,u,t);

% Comparison of step responses
figure(fig); fig = fig + 1;

plot(t, y, 'k', 'LineWidth', 1); hold on;
plot(t,f_Oldenbourg, '--b', 'LineWidth', 1);
plot(t,y_Oldenbourg, '--m', 'LineWidth', 1);
plot(t, y_smith, '-.r', 'LineWidth', 1);
xlabel('Time (s)','interpreter','Latex');
ylabel('Output','interpreter','Latex');
grid
leg1 = legend({'Real output','Oldenbourg and Sartorius (1948) f','Oldenbourg and Sartorius (1948) tf','Smith (1957)'},'FontName','Times New Roman','FontSize',12,'location','southeast');
set(leg1(1),'Interpreter','latex');
legend boxoff
