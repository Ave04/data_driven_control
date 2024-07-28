%% Housekeeping
clear, clc

%% Initialise matrix rows/columns
n = 4; m = 2;
A = randn(n,n); B = randn(n,m);

%% Set sequence size and cell array (optional)
T = 15;
xk = {}; % create cell for state space 

%% Populate state + input matrices
xk = randn(n,1);
uk = randn(m,1);

%% DT loop (matrix)
xk1 = [];
xk1(1) = xk(1);

for k = 1:T-1
     xk1(k+1) = A*xk(k) + B*uk(k);  
     disp(xk1(k));
end

%% Populate Cell arrays
for i = 1:T
    xk{i} = randn(n,1); % **wrong** - successive data points depend on the previous ones
end

uk = {};
for i = 1:T
    uk{i} = randn(m,1);
end

%% DT loop (cells)
xk1{1} = xk{1};

for k = 1:T-1
     xk1{k+1} = A*xk{k} + B*uk{k};  
     disp(xk1{k});
end

%% CVX soln

cvx_begin
    
cvx_end
