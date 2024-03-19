function [S,T,distance,time_seconds,Q,infotable] = nearest_hurwitz_stable(A, B, maxiter, timemax_seconds, x0)
% Computes a (locally) nearest pencil T*x + S to the 
% square pencil B*x + A with eigenvalues in the prescribed set
% 
% Input:
%   A, B (matrix of size (n,n))
%       the coefficients of the pencil B*x + A in consideration
%   maxiter (integer)
%        maximum amount of outer iterations, used in Manopt's
%        "trustregions" solver
%   timemax_seconds (integer)
%        maximum amount of time (in seconds) the algorithm can run before
%        the iteration is stopped
%   x0 (tensor of size (n,n,2))
%       contains the initial value of the optimization variables s.t. x0(:,:,1)
%       and x0(:,:,2) are unitary
%
% Output:
%   S, T (matrix of size (n,n))
%       the coefficients of the nearest pencil T*x + S with eigenvalues in
%       the prescribed set
%   distance (vector)
%       the distance to the constructed pencil after each iteration
%   time_seconds (vector)
%       elapsed time after each iteration 
%   Q (tensor of size (n,n,2))
%       contains the final value of the optimization variables s.t. Q(:,:,1)
%       and Q(:,:,2) are unitary
%   infotable (table)
%       contains various types of information for diagnostic purposes
%
% Requirement: Manopt needs to be imported

n = length(B);

if not(exist('maxiter', 'var'))
    maxiter = 1000;
end
if not(exist('timemax_seconds', 'var'))
    timemax_seconds = 1000;
end
if not(exist('x0', 'var'))
    x0 = [];
end

% Rescale the pencil to be of norm 100
P_norm = norm([A B], 'f')*1e-2;
A = A / P_norm;
B = B / P_norm;

problem.M = stiefelcomplexfactory(n, n, 2);
problem.cost = @cost;
% The code uses the Euclidean gradient. Projection to 
% the tangent space of U(n) is handled automatically (see 
% stiefelcomplexfactory documentation)
problem.egrad = @egrad;
% Euclidean Hessian. Projection is handled automatically.
% problem.ehess = @ehess;

options.tolgradnorm = 1e-10;
options.maxiter = maxiter;
options.maxtime = timemax_seconds;
options.verbosity = 2; % 2 = Default; 0 = No output; 

[Q, xcost, info, ~] = trustregions(problem, x0, options);

infotable = struct2table(info);
distance = sqrt(infotable.cost);
time_seconds = infotable.time;

% Construct the nearest singular pencil
T = Q(:,:,1)*B*Q(:,:,2);
S = Q(:,:,1)*A*Q(:,:,2);

alphas = (abs(diag(T)).^2 + abs(diag(S)).^2) ./...
    (2*(real(diag(T).*conj(diag(S))) ));

signs = real(diag(T).*conj(diag(S))) < 0;
lambdas = signs.*(alphas + sqrt(alphas.^2 - 1));

% The case lambda = -1 should be dealt with separately.
if any(abs(lambdas + 1) < 1e-6)
    keyboard
end

[S, T] = proj(S, T);

T = Q(:,:,1)'*T*Q(:,:,2)';
S = Q(:,:,1)'*S*Q(:,:,2)';

% Squared distance to the singular pencil. Should be equal to xcost
assert(abs(xcost -  (norm(B-T,'f')^2 + norm(A-S,'f')^2)) < 1e-10)

% Rescale back
T = T*P_norm;
S = S*P_norm;
distance = distance*P_norm;

% ---------

function f = cost(Q)

    T = Q(:,:,1)*B*Q(:,:,2);
    S = Q(:,:,1)*A*Q(:,:,2);

    alphas = (abs(diag(T)).^2 + abs(diag(S)).^2) ./...
    (2*(real(diag(T).*conj(diag(S)))));

    signs = real(diag(T).*conj(diag(S))) < 0;
    lambdas = signs.*(alphas + sqrt(alphas.^2 - 1));

    f = norm(tril(T,-1),'fro')^2 + norm(tril(S,-1),'fro')^2 + real(diag(T).*conj(diag(S))).'*lambdas;
end

function g = egrad(Q)
    
    Q1 = Q(:,:,1);
    Q2 = Q(:,:,2);
    
    M11 = B*Q2;
    M01 = A*Q2;
    
    M12 = Q1*B;
    M02 = Q1*A;
    
    T = Q1*B*Q2;
    S = Q1*A*Q2;
    
    [PS, PT] = proj(S, T);

    g = zeros(size(Q));

    g(:,:,1) = 2* (T - PT) * M11' + 2* (S - PS) * M01';
    g(:,:,2) = 2* M12' * (T - PT) + 2* M02' * (S - PS);

end

function H = ehess(Q, d)

    % not implemented yet

end

function [PS, PT] = proj(S, T)

    alphas = (abs(diag(T)).^2 + abs(diag(S)).^2) ./...
    (2*(real(diag(T).*conj(diag(S)))));

    signs = real(diag(T).*conj(diag(S))) < 0;
    lambdas = signs.*(alphas + sqrt(alphas.^2 - 1));

    PS = triu(S,1) + diag(1./(1-lambdas.^2).* (diag(S) - lambdas.*diag(T)));
    PT = triu(T,1) + diag(1./(1-lambdas.^2).* (diag(T) - lambdas.*diag(S)));
    
end

end
