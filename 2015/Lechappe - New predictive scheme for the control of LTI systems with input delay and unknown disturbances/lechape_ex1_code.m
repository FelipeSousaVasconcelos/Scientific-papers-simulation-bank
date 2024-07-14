% @article{lechappe2015new,
%   title={New predictive scheme for the control of LTI systems with input delay and unknown disturbances},
%   author={L{\'e}chapp{\'e}, Vincent and Moulay, Emmanuel and Plestan, Franck and Glumineau, Alain and Chriette, Abdelhamid},
%   journal={Automatica},
%   volume={52},
%   pages={179--184},
%   year={2015},
%   publisher={Elsevier}
% }


% Felipe José de Sousa Vasconcelos
% First modification: 12/07/2024
% Last modification: 12/07/2024

%% Initialization
clear all
close all
clc

warning off

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot,'defaultTextInterpreter','latex');
rng('default');

%% simulation parameters
Tsim = 25;
fig = 1;

%% System
A = [0 1; -9 3];
B = [0 1]';
h = 0.5;

% I am not sure what initial conditions the authors actually used!
x0_artstein = [1.5 -23]';
x0_lechappe = [1.5 -23]';

%% controller
kp = 45;
kd = 18;
ki = 60;

% Zhong integral
N = 500;
[tal,Q,KSoma,Ad,Bd,Cd,Dd] = paramZhong(N,h,A,B);

%% Simulation
sim('lechape_ex1_simu.slx');

%% Figures
figure(fig); fig = fig + 1;
plot(t, x_artstein(:,1), 'k' , 'LineWidth', 2);
hold on
plot(t, xp_artstein(:,1), 'b', 'LineWidth', 2);
plot(t, ref, '--', 'Color', "#ff6600", 'LineWidth', 1);
xlabel('time (s)','interpreter','Latex');
leg1 = legend({'$x_1$\hspace{50pt}','$x_{\hat{p}1}$\hspace{50pt}','$d\hspace{50pt}$'},'FontName','Times New Roman','FontSize',12,'location','northeast');
set(leg1(1),'Interpreter','latex');
axis([0 25 -4.7 4.7])
grid

figure(fig); fig = fig + 1;
plot(t, x_lechappe(:,1), 'k' , 'LineWidth', 2);
hold on
plot(t, xp_lechappe(:,1), 'r', 'LineWidth', 2);
plot(t, ref, '--', 'Color', "#ff6600", 'LineWidth', 1);
xlabel('time (s)','interpreter','Latex');
leg1 = legend({'$x_1$\hspace{50pt}','$X_{\hat{p}1}$\hspace{50pt}','$d\hspace{50pt}$'},'FontName','Times New Roman','FontSize',12,'location','northeast');
set(leg1(1),'Interpreter','latex');
axis([0 25 -4.7 4.7])
grid

