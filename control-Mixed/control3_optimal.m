function [Tx, X, U1, U2] = control3_optimal(pars, tFinal, B, u1Max, u2Max)

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

options = odeset('RelTol', 1e-7, 'AbsTol',[1e-7 1e-7 1e-7]);

tf = tFinal;

%-- optimal control: FBSM --%

u1_min = 0;
u1_max = u1Max;

u2_min = 0;
u2_max = u2Max;

n = tf*5+1;

u1 = u1_max*ones(1,n);
u2 = u2_max*ones(1,n);

Tu = linspace(0, tf, n);

max_iteration = 300;
TOL = 1e-2;

for i = 1:max_iteration
    % forward solve: stateEq (0, tf), x0=initx
    [Tx,X] = ode45(@(t,x) stateEq(t,x,u1,u2,Tu), [0 tf], initx, options);
    
    % save ONLY if u=zeros (states without control)
    if i==1
        Tx0 = Tx;
        X0 = X;
        [Tx,X] = ode45(@(t,x) stateEq(t,x,u1,u2,Tu), [0 tf], initx, options);
    end
   
    % backward solve: costateEq, (tf, 0), l0=initp
    x1 = X(:,1);
    x2 = X(:,2);
    x3 = X(:,3);

    [Tp,P] = ode45(@(t,p) costateEq(t,p,u1,u2,x1,x2,x3,Tu,Tx), [tf 0], initp, options);
    
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
    
    u1_old = u1;
    u2_old = u2;

    % third cost functional
    u1 = (x3.*p3 + x1.*p1 + p2.*(2.*x2 - x3 + x1))/(2*B);
    u2 = -(p2.*(x2.*h.*phi - x1.*h.*rho.*(x1/Xe - 1)) - p3.*(x2.*h*phi - x1.*h.*rho.*(x1/Xe - 1)))/(2*B);


    % applying bounds
    u1 = min(u1_max, u1);
    u1 = max(u1_min, u1);

    u2 = min(u2_max, u2);
    u2 = max(u2_min, u2);

    ERROR = max(norm(u2_old - u2), norm(u1_old - u1));

    u1 = 0.25*u1+0.75*u1_old;
    u2 = 0.25*u2+0.75*u2_old;
    
    disp([i, ERROR]);

    % convergence criterion
    if ERROR < TOL
       disp(i)
       U1 = interp1(Tu,u1,Tx);
       U2 = interp1(Tu,u2,Tx);
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
function dx = stateEq(t,x,u1,u2,Tu)
dx = zeros(3,1);

X1 = x(1);
P01 = x(2);
II = x(3);

uu1 = interp1(Tu,u1,t);
uu2 = interp1(Tu,u2,t);

dx(1) = (s + mu)*II - (2*mu + rho + s)*X1 - (gama + uu1)*X1;
dx(2) = rho*(1 - h*(1-uu2))*X1*(1 - X1/Xe) - (s + phi*h*(1-uu2) + 2*mu)*P01 + (gama + uu1)*(II - X1 - 2*P01);
dx(3) = rho*h*(1-uu2)*X1*(1 - X1/Xe) + phi*h*(1-uu2)*P01 - mu*II - (gama + uu1)*II;
end

%%-- Costate equations --%%
function dp = costateEq(t,p,u1,u2,x1,x2,x3,Tu,Tx)
dp = zeros(3,1);

x1 = interp1(Tx,x1,t);
x2 = interp1(Tx,x2,t);
x3 = interp1(Tx,x3,t);

uu1 = interp1(Tu,u1,t);
uu2 = interp1(Tu,u2,t);

X1 = x1;
P01 = x2;
II = x3;

l1 = p(1);
l2 = p(2);
l3 = p(3);

dp(1) = l2*(gama + uu1 + rho*(X1/Xe - 1)*(h*(uu2 - 1) + 1) + (X1*rho*(h*(uu2 - 1) + 1))/Xe) - l3*(h*rho*(X1/Xe - 1)*(uu2 - 1) + (X1*h*rho*(uu2 - 1))/Xe) + l1*(gama + 2*mu + rho + s + uu1);

dp(2) = l2*(2*gama + 2*mu + s + 2*uu1 - h*phi*(uu2 - 1)) + h*l3*phi*(uu2 - 1);

dp(3) = l3*(gama + mu + uu1) - l1*(mu + s) - l2*(gama + uu1) - 1;

end

end
