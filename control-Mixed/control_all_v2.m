%--- MASTER.m ---%
clear all;

% parameter names
%pars = [rho, s, phi, mu, N, h, gama];

%-- parameter sets --%

% HPV %
%pars = [5, 2, 52, 1/9, 1e6, 0.073, 0.5];

% trichonomiasis %
pars = [5, 2, 52, 1/9, 1e6, 0.115, 0.727];
% h = (19.2 + 3.86)/2 %
% 1 year = 52 weeks
% 1/gamma = (2.08 + 14.6 + 20.1 + 249)/4 weeks
% 1/gamma = 1.37 years
% gamma = 1.373942308

% gonorrhea %
%pars = [5, 2, 52, 1/9, 1e6, 0.348, 1.538];
% h = (45.9 + 23.7)/2 %
% 1 year = 52 weeks
% 1/gamma = (34.0 + 33.6)/2 weeks
% 1/gamma = 0.65 years
% gamma = 1.538461538

% chlamydia %
%pars = [5, 2, 52, 1/9, 1e6, 0.129, 0.855];
% h = (16.2 + 9.75)/2 %
% 1 year = 52 weeks
% 1/gamma = (15.0 + 106.6)/2 weeks
% 1/gamma = 1.169 years
% gamma = 0.855263158

%-- global parameters --%

totalYears = 20.0;
u1Max = 1.0;
u2Max = 0.75; 
cost = 5e5;

%------------------------%
%-- constant functions --%
%------------------------%

% "function [Tx0, X0] = control1_constant(pars, tFinal, uConst)"

% STEP 1: Computation

[Tx0, X0] = control1_constant(pars, totalYears, 0.0);
I0 = X0(:,3);

[Tx1Old, X1Old, U11Old, U12Old] = control3_optimal(pars, totalYears, 1e5, u1Max, u2Max);

Tx1 = linspace(Tx1Old(1), Tx1Old(end), 100);

X1 = interp1(Tx1Old, X1Old, Tx1);

U1  = interp1(Tx1Old, U11Old, Tx1);
U2  = interp1(Tx1Old, U12Old, Tx1);

I1 = X1(:,3);

% STEP 2: Plots

%figure('PaperPositionMode', 'auto');

figure;

% 2.1: control
subplot(1,2,1);
hold on;

plot(Tx1, U1, '-r', 'LineWidth',2);
plot(...
Tx1(1:floor(size(Tx1,2)/10):size(Tx1,2)),...
 U1(1:floor(size(Tx1,2)/10):size(Tx1,2)),...
'xr');

plot(Tx1, U2, '-m', 'LineWidth',2);
plot(...
Tx1(1:floor(size(Tx1,2)/10):size(Tx1,2)),...
 U2(1:floor(size(Tx1,2)/10):size(Tx1,2)),...
'dm');

ylabel('Control functions $u_T$, $u_C$','Interpreter','latex');
ylim([0.0, u1Max]);
yticks([0, u1Max/2, u1Max]);

xlabel('Time $t$','Interpreter','latex');
xlim([0, totalYears]);
xticks([0, totalYears/2, totalYears]);

set(gca,'TickLabelInterpreter','latex');
box on; 

h    = zeros(2, 1);
h(1) = plot(0,NaN,'xr');
h(2) = plot(0,NaN,'dm');

%axP2 = get(gca,'Position'); 

lh   = legend(h, ...
             {'$u_T$: Treatment', '$u_C:$ Condom'}, ...
              'Interpreter','latex', ...
              'Location','NorthEast', ...
              'Box','off');
          
%set(gcf,'units','points','position',[0,0,600,400]);

% 2.2: infected population
subplot(1,2,2);
hold on;
plot(Tx0,X0(:,3),':k','LineWidth',2);

plot(Tx1,I1,'-k','LineWidth',2);

ylabel('Infected individuals $I$','Interpreter','latex');
ylim([0., 1e6]);
yticks([0, 5e5, 1e6]);

xlabel('Time $t$','Interpreter','latex');
xlim([0, totalYears]);
xticks([0, totalYears/2, totalYears]);

set(gca,'TickLabelInterpreter','latex');
box on; 

h    = zeros(2, 1);
h(1) = plot(0,NaN,':k');
h(2) = plot(0,NaN,'-k');

%axP2 = get(gca,'Position'); 

lh   = legend(h, ...
             {'No control','Control'}, ...
              'Interpreter','latex', ...
              'Location','NorthEast', ...
              'Box','off');

set(gcf, 'Renderer', 'painters', 'Position', [0 0 750 200]);
%%
disp('Cost:');
disp(trapz(Tx1, cost*U1.^2 + cost*U2.^2));
disp('Averted people:');
disp(X0(end,3) - I1(end));
