% @article{manitius1979finite,
%   title={Finite spectrum assignment problem for systems with delays},
%   author={Manitius, AWOAZ and Olbrot, A},
%   journal={IEEE transactions on Automatic Control},
%   volume={24},
%   number={4},
%   pages={541--552},
%   year={1979},
%   publisher={IEEE}
% }

% Felipe José de Sousa Vasconcelos
% First modification: 28/06/2024
% Last modification: 30/06/2024

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
Tsim = 30;
Tref = 0;
ref = 1;

s = tf('s');
fig = 1;

%% System
a = 2;
b0 = 0;
b1 = 1;
h = 2;

%% f
% lamda = a+f*exp(-a*h)
p = -1;
ba = b0 + expm(-a*h)*b1;
f = (p-a)/ba;

%% Closed loop system
CL = tf(b1,[1 -(a+ba*f)]);

%% Simulation
w = f*(exp(-a*h)-exp(-h*s))/(s-a);

% W integral
S1 = exp(-a*h)/(s-a);

% Model SS S1
dssP1 = ss(S1, 'min');
A1z = dssP1.A; %Aj=T^-1*A*T
B1z = dssP1.B; %Bj=T^-1*B;
C1z = dssP1.C; %Cj=C*T;

% Compute the integral 
N=2000;
[tal1,Q1,KSoma1,Ad1,Bd1,Cd1,Dd1] = paramW(N,h,A1z,B1z);

sim('manitius_ex2_simu.slx');

%% Figures
figure(fig); fig = fig + 1;

subplot(2,1,1)
plot(t, y_open_loop, 'k', 'LineWidth', 2);
xlabel('Time (s)','interpreter','Latex');
ylabel('Output','interpreter','Latex');
title('Open loop response','interpreter','Latex');
grid
subplot(2,1,2)
plot(t,y_closed_loop, 'k', 'LineWidth', 2); 
hold on
plot(t,y_closed_loop_tf, '--r', 'LineWidth', 2);
xlabel('Time (s)','interpreter','Latex');
ylabel('Output','interpreter','Latex');
leg1 = legend({'Closed loop response','Closed loop response tf'},'FontName','Times New Roman','FontSize',12,'location','southeast');
set(leg1(1),'Interpreter','latex');
legend boxoff
grid
