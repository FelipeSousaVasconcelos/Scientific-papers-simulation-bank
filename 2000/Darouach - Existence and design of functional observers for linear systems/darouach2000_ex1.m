% @article{darouach2000existence,
%   title={Existence and design of functional observers for linear systems},
%   author={Darouach, Mohamed},
%   journal={IEEE Transactions on Automatic Control},
%   volume={45},
%   number={5},
%   pages={940--943},
%   year={2000},
%   publisher={IEEE}
% }

% Felipe José de Sousa Vasconcelos
% First modification: 05/07/2024
% Last modification: 05/07/2024

%% Initialization
clear all
close all
clc

warning off

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot,'defaultTextInterpreter','latex');
rng('default');

%% System
A = [-2 1 0; 1 -3 1; 0 0 -1];
C = [1 0 0; 0 1 0];
L = [1 1 0];

Pi = [L; C*A; C];
Theta = L*A;

Abar = A*(eye(3)-pinv(L)*L);
Cbar = C*(eye(3)-pinv(L)*L);
Sig = [C*Abar; Cbar];

%% Checking conditions

% Conditions discussed in ref [12]
if rank([L*A; C; L]) == rank([C; L])
    disp('Condition in [12] is satisfied.')
else
    disp('Condition in [12] is NOT satisfied.')
end

% Conditions discussed of Theorem 1
if rank([L*A; C*A; C; L]) == rank([C*A; C; L])
    disp('Condition (10) of Theorem 1 is satisfied.')
else
    disp('Condition (10) of Theorem 1 is NOT satisfied.')
end

for sigma = 0:100
    for omega = -100:0.01:100
        if rank([(sigma+1i*omega)*L-L*A; C*A; C]) == rank([C*A; C; L])
            test = 1;
        else
            test = 0;
            break
        end
    end
end

if test == 1
    disp('Condition (20) of Theorem 1 is satisfied.')
else
    disp('Condition (20) of Theorem 1 is NOT satisfied.')
end


%% Design of the functional observer
F = L*A*pinv(L) - L*Abar*pinv(Sig)*[C*A; C]*pinv(L)
G = ( eye(4) - Sig*pinv(Sig)) * [C*A; C] *pinv(L)
p = -3;

Z4 = (F - p) / G(4);

Z = [0 0 0 Z4]

N = F - Z*G

EK = L*Abar*pinv(Sig) + Z*( eye(4) - Sig*pinv(Sig));
E = EK(1:2)
K = EK(3:4)

J = K+N*E

