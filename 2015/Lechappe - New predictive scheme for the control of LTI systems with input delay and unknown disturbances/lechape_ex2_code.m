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
x0 = [1.5; 1];

%% controller
K = [-45 -18];

% Zhong integral
N = 500;
[tal,Q,KSoma,Ad,Bd,Cd,Dd] = paramZhong(N,h,A,B);

%% Criteria (24) and (25) 
norm_expm_Ah = norm(expm(A*h));
dmax = 3; 
Dmax = 1.5;

c24 = dmax/Dmax > h*(norm(expm(A*h))/(norm(expm(A*h))-1));
c25 = dmax/Dmax > h;

if c24 == 1 && c25 == 1
    disp('Criteria checked!')
end

%% Simulation
sim('lechape_ex2_simu.slx');

%% Figures
figure(fig); fig = fig + 1;
subplot(2,1,1)
plot(t, ref, '--', 'Color', "#ff6600", 'LineWidth', 1);
xlabel('time (s)','interpreter','Latex');
ylabel('$d(t)$','interpreter','Latex');
grid
subplot(2,1,2)
plot(t, vecnorm(x_artstein,2,2), '-.b' , 'LineWidth', 2);
hold on
plot(t, vecnorm(x_lechappe,2,2), 'r', 'LineWidth', 2);
xlabel('time (s)','interpreter','Latex');
ylabel('$||x(t)||$','interpreter','Latex');
leg1 = legend({'$||x(t)||$ with (27)','$||x(t)||$ with (28)'},'FontName','Times New Roman','FontSize',12,'location','northeast');
set(leg1(1),'Interpreter','latex');
axis([0 25 0 4.7])
grid

