% @article{darouach2019existence,
%   title={On the existence and design of functional observers},
%   author={Darouach, Mohamed and Fernando, Tyrone},
%   journal={IEEE Transactions on Automatic Control},
%   volume={65},
%   number={6},
%   pages={2751--2759},
%   year={2019},
%   publisher={IEEE}
% }

% Felipe José de Sousa Vasconcelos
% First modification: 04/07/2024
% Last modification: 04/07/2024

%% System
A = [-2 1 0; 1 -3 1; 0 0 -1];
C = [1 0 0; 0 1 0];
L = [1 1 0];

Pi = [L; C*A; C];
Theta = L*A;

Abar = A*(eye(3)-pinv(L)*L);
Cbar = C*(eye(3)-pinv(L)*L);
Sig = [C*Abar; Cbar];

%% Method in ref. [1]

F = L*A*pinv(L) - L*Abar*pinv(Sig)*[C*A; C]*pinv(L)
G = ( eye(4) - Sig*pinv(Sig)) * [C*A; C] *pinv(L)
p = -3;

Z4_ref1 = (F - p) / G(4);

Z_ref1 = [0 0 0 Z4_ref1];

N_ref1 = F - Z_ref1*G

EK_ref1 = L*Abar*pinv(Sig) + Z_ref1*( eye(4) - Sig*pinv(Sig));
E_ref1 = EK_ref1(1:2)
K_ref1 = EK_ref1(3:4)

J_ref1 = K_ref1+N_ref1*E_ref1

%% Method in the article

Piinv = [ ( eye(3)-( eye(3)-pinv(L)*L ) * pinv(Sig) * [C*A; C] )*pinv(L)  ( eye(3)-pinv(L)*L )*pinv(Sig)];

A1_ = Theta*Piinv * [1 0 0 0 0]';
A1 = L*A*pinv(L) - L*Abar*pinv(Sig)*[C*A; C]*pinv(L)

B1_ = (eye(5)-Pi*Piinv) * [1 0 0 0 0]';
B1 = [0; -(eye(4) - Sig * pinv(Sig)) * [C*A; C] *pinv(L)]

Z4 = (A1 - p) / B1(5);
Z = [0 0 0 Z4]

N = A1 - Z*B1(2:end)

EK = L*Abar*pinv(Sig) - Z*( eye(4) - Sig*pinv(Sig));
E = EK(1:2)
K = EK(3:4)

J = K+N*E
