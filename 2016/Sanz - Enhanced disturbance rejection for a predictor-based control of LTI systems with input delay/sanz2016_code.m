% @article{sanz2016enhanced,
%   title={Enhanced disturbance rejection for a predictor-based control of LTI systems with input delay},
%   author={Sanz, Ricardo and Garcia, Pedro and Albertos, Pedro},
%   journal={Automatica},
%   volume={72},
%   pages={205--208},
%   year={2016},
%   publisher={Elsevier}
% }

% Felipe José de Sousa Vasconcelos
% First modification: 13/07/2024
% Last modification: 13/07/2024

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

%% Check condition (31)

w = 0.5;
if w<2/h
    disp('Condition (31) is checked!')
end

%% Observer matrix
w1 = 1.49;
r = 1;
CH = [1 h];
for j = 0:r
    c(j+1) = (factorial(r+1)/(factorial(j+1)))*w1^(j+1);
end
AH_w1 = [-c(1) 1
      -c(2) 0];
BH_w1 = [c(1) c(2)]';

Aw1 = AH_w1;
B1w1 = AH_w1*BH_w1*pinv(B)-BH_w1*pinv(B)*A;
B2w1 = BH_w1;
C1 = CH;
C2w1 = CH*BH_w1*pinv(B);

w2 = 10;
for j = 0:r
    c(j+1) = (factorial(r+1)/(factorial(j+1)))*w2^(j+1);
end
AH_w2 = [-c(1) 1
      -c(2) 0];
BH_w2 = [c(1) c(2)]';

Aw2 = AH_w2;
B1w2 = AH_w2*BH_w2*pinv(B)-BH_w2*pinv(B)*A;
B2w2 = BH_w2;
C2w2 = CH*BH_w2*pinv(B);

%% Simulation
sim('sanz2016_simu.slx');

%% Figures
figure(fig); fig = fig + 1;
subplot(2,1,1)
plot(t, ref, '--k', 'LineWidth', 2);
xlabel('time (s)','interpreter','Latex');
ylabel('$d(t)$','interpreter','Latex');
axis([0 25 -4 4])
grid
subplot(2,1,2)
plot(t, vecnorm(x_sanz_w1,2,2), ':b' , 'LineWidth', 1);
hold on
plot(t, vecnorm(x_sanz_w2,2,2), 'b' , 'LineWidth', 1);
plot(t, vecnorm(x_lechappe,2,2), '--r', 'LineWidth', 1);
xlabel('time (s)','interpreter','Latex');
ylabel('$||x(t)||$','interpreter','Latex');
leg1 = legend({'Proposed $\omega_0=1.49$','Proposed $\omega_0=10$','Lechappe'},'FontName','Times New Roman','FontSize',12,'location','northeast');
set(leg1(1),'Interpreter','latex');
axis([0 25 0 2])
grid






