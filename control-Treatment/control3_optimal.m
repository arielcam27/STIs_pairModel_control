function [Tx, X, U] = control3_optimal(pars, tFinal, B, uMax)

% INPUT:
% -pars: parameter values
% -tFinal: final time
% -B: weight parameter
% OUTPUT: 
% -(Tx, X): state variables at times Tx
% -U: control variable at times Tx (interpolated)

%-- model parameters --%

%pars = [rho, s, phi, mu, N, h, gama];
%pars = [5, 2, 52, 1/9, 1e6, 0.115, 0.72];

%-- behavior parameters --%
rho = pars(1);
s   = pars(2);
phi = pars(3);
mu  = pars(4);
N   = pars(5); %- N = nu/mu -> nu = N*mu

%-- STI parameters --%

% HPV -- no curable STI
%h = 0.073;
%gamma = 0.2;

% trichomoniasis
%h = 0.115;
%gama = 0.72;

h    = pars(6);
gama = pars(7);

%-- steady-state --%
Xe = N*(s + 2*mu)/(s + 2*mu + rho);

% note: with those parameters we have
% - analitically:
% X* = 307692.307692
% P* = 346153.846154
% - numerically: 
% X1*  = 228824.758596
% P01* = 60295.2300822
% II*   = 804499.674566

% Initial conditions
% X1, P01, II
initx = [1e5 3e5 7e5]; 
initp = [0. 0. 0.];

options = odeset('RelTol', 1e-8, 'AbsTol',[1e-8 1e-8 1e-8]);

tf = tFinal;

%-- optimal control: FBSM --%

u_min = 0;
u_max = uMax;

n = tf*100+1;
u = u_max*ones(1,n);
Tu = linspace(0, tf, n);

max_iteration = 50;
TOL = 1e-2;

for i = 1:max_iteration
    % forward solve: stateEq (0, tf), x0=initx
    [Tx,X] = ode45(@(t,x) stateEq(t,x,u,Tu), [0 tf], initx, options);
    
    % save ONLY if u=zeros (states without control)
    if i==1
        Tx0 = Tx;
        X0 = X;
        [Tx,X] = ode45(@(t,x) stateEq(t,x,u,Tu), [0 tf], initx, options);
    end
   
    % backward solve: costateEq, (tf, 0), l0=initp
    x1 = X(:,1);
    x2 = X(:,2);
    x3 = X(:,3);

    [Tp,P] = ode45(@(t,p) costateEq(t,p,u,x1,x2,x3,Tu,Tx), [tf 0], initp, options);
    
    % control update
    p1 = P(:,1);
    p2 = P(:,2);
    p3 = P(:,3);

    x1 = interp1(Tx,x1,Tu);
    x2 = interp1(Tx,x2,Tu);
    x3 = interp1(Tx,x3,Tu);

    p1 = interp1(Tp,p1,Tu);
    p2 = interp1(Tp,p2,Tu);
    p3 = interp1(Tp,p3,Tu);
    
    u_old = u;

    % third cost functional
    u = 0.5*(x1.*p1 + p2.*(2.0*x2 - x3 + x1) + x3.*p3)/B;

    % applying bounds
    u = min(u_max, u);
    u = max(u_min, u);
    ERROR = norm(u_old-u);
    u = 0.5*u+0.5*u_old;

    % convergence criterion
    if ERROR < TOL
       disp(i)
       U = interp1(Tu,u,Tx);
       break;
    end
end

% check if FBSM converged
if i == max_iteration
    disp('Stopped before required residual is obtained.');
end


%--------------------
% auxiliary functions
%--------------------

%%-- Model equations with control u(t) --%%
function dx = stateEq(t,x,u,Tu)
dx = zeros(3,1);
X1 = x(1);
P01 = x(2);
II = x(3);

uu = interp1(Tu,u,t);

dx(1) = (s + mu)*II - (2*mu + rho + s)*X1 - (gama + uu)*X1;
dx(2) = rho*(1 - h)*X1*(1 - X1/Xe) - (s + phi*h + 2*mu)*P01 + (gama + uu)*(II - X1 - 2*P01);
dx(3) = rho*h*X1*(1 - X1/Xe) + phi*h*P01 - mu*II - (gama + uu)*II;
end

%%-- Costate equations --%%
function dp = costateEq(t,p,u,x1,x2,x3,Tu,Tx)
dp = zeros(3,1);

x1 = interp1(Tx,x1,t);
x2 = interp1(Tx,x2,t);
x3 = interp1(Tx,x3,t);

uu = interp1(Tu,u,t);

X1 = x1;
P01 = x2;
II = x3;

l1 = p(1);
l2 = p(2);
l3 = p(3);

% third cost functional
dp(1) = -l1*(-gama - 2.0*mu - rho - s - uu) - l2*(-X1*rho*(-h + 1.0)/Xe - gama + rho*(-h + 1.0)*(-X1/Xe + 1.0) - uu) - l3*(-X1*h*rho/Xe + h*rho*(-X1/Xe + 1.0));

dp(2) = -h*l3*phi - l2*(-2.0*gama - h*phi - 2.0*mu - s - 2.0*uu);

dp(3) = -l1*(mu + s) - l2*(gama + uu) - l3*(-gama - mu - uu) - 1;

end

end
