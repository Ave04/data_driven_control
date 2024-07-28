clear all, close all, clc


%% Initialise system 
n = 2 ; m = 1;

%A = randn(n,n) ; %B = randn(n,m);

A = [0.9 0;0 1.1];
B = [0;1];
K = [0 0.3];
CL = A - B*K;  % closed-loop

%% Generate Data

% Sample size
T = 10;

% initial value
x(:,1) = [10;10];

for k=1:T
    u(:,k) = -K*x(:,k); % randomize -> may not meet rank condition
    x(:,k+1) = A*x(:,k) + B*u(:,k);
end


disp(x)

%plot(x')

%% DMDc

X = x(:,1:end-1); % X0,T
X2 = x(:,2:end); % X1,T

U = u;
Acorrect = (X2 - B*U)*pinv(X);

%% CVX solver

cvx_begin sdp
    variable Q(T,n)
    %main constraint
    [X*Q X2*Q;Q'*X2' X*Q] > 10^(-9)*eye(2*n)
cvx_end

%% Solve for K

K2 = U*Q*inv(X*Q); % only 1 of the stabilising K (can have multiple)

RHS = A + B*K2; 
Gk = pinv(X2)*RHS % Decision variable Gk

%% Right-Inverse

B2 = X2*pinv(U);
A2 = X2*pinv(X);

%% LQR

% W is the controllability grammian --> use 'gram' to compute --> n x n
% size
% X here is R^(0.5)*K*W*K'*R^(0.5) --> R is specified
cvx_begin
    variables Q(n,n) W(n,n)  ;
    minimize( trace(Qx*W + R^(0.5)*K*W*K'*R^(0.5)) ) 
    subject to
        (A+B*K)*W*(A+B*K)'-W+eye(n) <= 0
        W >= eye(n)
        X-R^(0.5)*K*W*K'*R^(0.5) >= 0
cvx_end
