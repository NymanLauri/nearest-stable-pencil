
maxiter = 2000;
maxtime = 100;

n = 10;

%rng(2,"twister") with real input, n=10 -> all real parts are 0
% rng(3,"twister") with real input, n=10 -> many eigenvalues are 0
rng(4,"twister")

%Real input
A = randn(n);
B = randn(n);

%Complex input
% A = randn(n) + 1i*randn(n) ;
% B = randn(n) + 1i*randn(n);


[S,T,distance,time_seconds,Q,infotable] = nearest_hurwitz_stable(A, B,maxiter,maxtime);

S_tri = Q(:,:,1)*S*Q(:,:,2);
T_tri = Q(:,:,1)*T*Q(:,:,2);

% eig(S_tri,-T_tri)

%Compute the eigenvalues
-diag(S_tri)./diag(T_tri)


