
maxiter = 2000;
maxtime = 70;

n = 10;

rng(5,"twister")

%Real input
% A = randn(n);
% B = randn(n);

%Complex input
A = randn(n) + 1i*randn(n) ;
B = randn(n) + 1i*randn(n);

[S,T,distance,time_seconds,Q,infotable] = nearest_schur_stable(A, B,maxiter,maxtime);

S_tri = Q(:,:,1)*S*Q(:,:,2);
T_tri = Q(:,:,1)*T*Q(:,:,2);

% eig(S_tri,-T_tri)

%Compute the eigenvalues
-diag(S_tri)./diag(T_tri)


