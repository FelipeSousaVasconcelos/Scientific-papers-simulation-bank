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
Tsim = 20;
Tref = 0;
ref = 1;
fig = 1;

%% System
syms mu
A = [0 0; 1 1];
B0 = [0; 0];
B1 = [1; 0];
C = [1 0];

Bhat = B0+B1*mu;
det([Bhat A*Bhat])

%% B(A)
syms h 

BA = B0 + expm(-A*h)*B1

h = 2; BA = double(subs(BA));
% Controllability
Co = ctrb(A,BA)
unco = length(A) - rank(Co)
if unco == 0
    disp('The pair (A,B(A)) is controllable');
else
    disp('The system is not controllable');
end

%% F
syms h s f1 f2
BA = B0 + expm(-A*h)*B1
F = [f1 f2];

M = eye(size(A,1))*s - A - BA*F;
detM = simplify(det(M));
eqn = detM - (s + 1)^2;

coeffsMatrix = [coeffs(eqn,s)];

solutions = solve(coeffsMatrix, [f1 f2]);

F = [simplify(solutions.f1) solutions.f2]

%% Closed loop transfer function
simplify(F*inv(eye(2)*s-A-BA*F)*B1)

%% Simulation
h = 2; F = double(subs(F)); BA = B0 + expm(-A*h)*B1;
s = tf('s');

[b,a] = ss2tf(A+BA*F,B1,F,0);

YU = tf(b,a);

% W1 integral
S1 = 1/s;

% Model SS S1
dssP1 = ss(S1, 'min');
A1z = dssP1.A; %Aj=T^-1*A*T
B1z = dssP1.B; %Bj=T^-1*B;
C1z = dssP1.C; %Cj=C*T;

% Compute the integral 
N=1000;
[tal1,Q1,KSum1,Ad1,Bd1,Cd1,Dd1] = paramW(N,h,A1z,B1z);

% W2 integral
S2 = 4/(s-1);

% Model SS S2
dssP2 = ss(S2, 'min');
A2z = dssP2.A; %Aj=T^-1*A*T
B2z = dssP2.B; %Bj=T^-1*B;
C2z = dssP2.C; %Cj=C*T;

% Compute the integral 
[tal2,Q2,KSum2,Ad2,Bd2,Cd2,Dd2] = paramW(N,h,A2z,B2z);

sim('manitius_ex1_simu.slx');

%% Figures
figure(fig); fig = fig + 1;

subplot(2,1,1)
plot(t, y_open_loop, 'k', 'LineWidth', 2);
xlabel('Time (s)','interpreter','Latex');
ylabel('Output','interpreter','Latex');
title('Open loop response','interpreter','Latex');
grid
subplot(2,1,2)
plot(t,y_closed_loop_SS, 'k', 'LineWidth', 2);
hold on
plot(t,y_closed_loop_FT, '--r', 'LineWidth', 2);
xlabel('Time (s)','interpreter','Latex');
ylabel('Output','interpreter','Latex');
grid
leg1 = legend({'Closed loop state space response','Closed loop tranfer function response'},'FontName','Times New Roman','FontSize',12,'location','southeast');
set(leg1(1),'Interpreter','latex');
legend boxoff
