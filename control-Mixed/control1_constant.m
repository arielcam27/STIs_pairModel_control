function [Tx0, X0] = control1_constant(pars, tFinal, uConst)

% INPUT:
% -pars: parameter values
% -tFinal: final time
% -uConst: constant control
% OUTPUT: 
% -(Tx, X): state variables at times Tx

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

%-- integration parameters --%
options = odeset('RelTol', 1e-8, 'AbsTol',[1e-8 1e-8 1e-8]);

t0 = 0; tf = tFinal;
n = tf*100+1;
Tu = linspace(t0, tf, n);

% initial conditions
initx = [1e5 3e5 7e5];

% control initial value
u = uConst*ones(1,n);

% integration
[Tx0,X0] = ode45(@(t,x) stateEq(t,x,u,Tu), [t0 tf], initx, options);

%---end---%

%---------------------%
% auxiliary functions %
%---------------------%

%%-- model equations with control u(t) --%%
function dx = stateEq(t,x,u,Tu)
dx = zeros(3,1);
X1 = x(1);
P01 = x(2);
II = x(3);

uu = interp1(Tu,u,t);

dx(1) = (s + mu)*II - (2*mu + rho + s)*X1 - gama*X1;
dx(2) = rho*(1 - h*(1-uu))*X1*(1 - X1/Xe) - (s + phi*h*(1-uu) + 2*mu)*P01 + gama*(II - X1 - 2*P01);
dx(3) = rho*h*(1-uu)*X1*(1 - X1/Xe) + phi*h*(1-uu)*P01 - mu*II - gama*II;
end

end
