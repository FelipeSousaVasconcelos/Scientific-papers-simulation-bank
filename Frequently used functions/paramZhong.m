function [tal,Q,KSoma,Ad,Bd,Cd,Dd] = paramZhong(N,L,A,B)
tal = L/N;

Q = cell(N,1);
for e = 1:N
    i = e-1;
    Q{e,1} = expm(((i/N)*L)*A)*B;
end

if (0==det(A))
    KSoma=0;
    for nn=1:5
        KSoma=KSoma+(L/N)*(A*(L/N))^(nn-1)/factorial(nn);
    end
else
    KSoma = (expm(A*(L/N))-eye(size(A)))/A;
end
% Control signal buffer
% for past d_up
Ad = zeros(N,N);
for e = 2:N
    Ad(e,e-1) = 1;
end
Bd = zeros(N,1); Bd(1,1) = 1;
Cd = eye(N,N); Cd = Ad;
Dd = zeros(N,1); Dd(1,1) = 1;
