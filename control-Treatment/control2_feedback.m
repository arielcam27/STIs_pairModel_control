function [Tx, X, U, PhiShape] = control2_feedback(pars, tFinal, I0, M)

% INPUT:
% -pars: parameter values
% -tFinal: final time
% -I0: alert level parameter
% OUTPUT: 
% -(Tx, X): state variables at times Tx
% -U: control variable at times Tx
% -PhiShape: Phi([0,N]) discretized with steps 10

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

options = odeset('RelTol', 1e-8, 'AbsTol',[1e-8 1e-8 1e-8]);

tf = tFinal;

[Tx,X] = ode45(@(t,x) stateEq(t,x), [0 tf], initx, options);

U = Phi(X(:,3)); % actual control function values for simulation

IIdummie = 0:10:N;
PhiShape = Phi(IIdummie); % general shape of Phi

%---------------------%
% auxiliary functions %
%---------------------%

%-- feedback control --%
function res = Phi(II)
%M = 0.5;
k = 1e-4;

res = M./(1+exp(-k*(II-I0))) - M./(1+exp(k*I0));

% k = (1/I0)*log(10*M-1);
% M = 0.8;
% n = 4;
% k = 1./(1-I0);
% res = M*(II.^n./(1-II))./(k*I0.^n+(II.^n./(1-II)));
end

%%-- model equations with control u(t) --%%
function dx = stateEq(t,x)
dx  = zeros(3,1);

X1  = x(1);
P01 = x(2);
II  = x(3);
uu  = Phi(II);

dx(1) = (s + mu)*II - (2*mu + rho + s)*X1 - (gama + uu)*X1;
dx(2) = rho*(1 - h)*X1*(1 - X1/Xe) - (s + phi*h + 2*mu)*P01 + (gama + uu)*(II - X1 - 2*P01);
dx(3) = rho*h*X1*(1 - X1/Xe) + phi*h*P01 - mu*II - (gama + uu)*II;
end

end
