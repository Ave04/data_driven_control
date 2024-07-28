clear all, close all, clc

%% Initialise system 
n = 2 ; m = 1;

%A = randn(n,n) ; %B = randn(n,m);

A = [0.9 0;0 1.1];
B = [0;1];
% K = [1 0.83577];
% CL = A - B*K;  % closed-loop

%% Generate Data

% Sample size
T = 10;

% initial value
x(:,1) = [10;10];

for k=1:T
    u(:,k) = rand(1);
    % u(:,k) = -K*x(:,k); % randomize -> may not meet rank condition
    x(:,k+1) = A*x(:,k) + B*u(:,k);
end

disp(x)

% plot(x')

%% DMDc

X = x(:,1:end-1); % X0,T
X2 = x(:,2:end); % X1,T

U = u;
Acorrect = (X2 - B*U)*pinv(X);

%% CVX solver for K (equation 15)

cvx_begin sdp
    variable Q(T,n)
    %main constraint
    [X*Q X2*Q;
    Q'*X2' X*Q] > 10^(-9)*eye(2*n)
cvx_end

%% Solve for K

K2 = U*Q*inv(X*Q) % only 1 of the stabilising K (can have multiple)

RHS = A + B*K2; 
Gk = pinv(X2)*RHS % Decision variable Gk

%% Right-Inverse (* might need correction => does not return same value of B)

B2 = X2*pinv(U);
A2 = X2*pinv(X);

%% LQR Q and R

Qx = diag([20,20]);  % same size as A (n by n)
R = 0.01*eye(m);  % size = (m by m)
Rroot = R^0.5;

% check dlqr K
[Klqr,~,~] = dlqr(A,B,Qx,R)

%% LQR Convex Optimisation (equation 27)
% W is the controllability grammian --> use 'gram' to compute --> n x n
% size
% Xdum here is R^(0.5)*K*W*K'*R^(0.5) --> R is specified

cvx_begin sdp
    variables Q2(T,n) Xdum(m,m)  
    minimize( trace(Qx*X*Q2) + trace(Xdum) ) 
    subject to
        [Xdum Rroot*U*Q2;Q2'*U'*Rroot X*Q2] >= 0
        [X*Q2-eye(n) X2*Q2;Q2'*X2' X*Q2] >= 0
cvx_end

%% Solve for K

K3 = U*Q2*inv(X*Q2)
