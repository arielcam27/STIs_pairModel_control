%--- MASTER.m ---%
clear all;

% parameter names
%pars = [rho, s, phi, mu, N, h, gama];

%-- parameter sets --%

% HPV %
%pars = [5, 2, 52, 1/9, 1e6, 0.073, 0.5];

% trichonomiasis %
%pars = [5, 2, 52, 1/9, 1e6, 0.115, 0.727];
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
pars = [5, 2, 52, 1/9, 1e6, 0.129, 0.855];
% h = (16.2 + 9.75)/2 %
% 1 year = 52 weeks
% 1/gamma = (15.0 + 106.6)/2 weeks
% 1/gamma = 1.169 years
% gamma = 0.855263158

%-- global parameters --%

totalYears = 20.0;
uMax = 1.0;
alert = 1e5;
cost = 5e5;

%------------------------%
%-- constant functions --%
%------------------------%

% "function [Tx0, X0] = control1_constant(pars, tFinal, uConst)"

% STEP 1: Computation

[Tx0, X0] = control1_constant(pars, totalYears, 0.0);

[TxOld, XOld] = control1_constant(pars, totalYears, uMax);

Tx = linspace(TxOld(1), TxOld(end), 100);
X  = interp1(TxOld, XOld, Tx);
U  = uMax*ones(1,size(Tx,2));

II = X(:,3)'; 

% STEP 2: Plots

%figure('PaperPositionMode', 'auto');

figure;

% 2.1: control
subplot(1,2,1);
hold on;

plot(Tx, U, '-k', 'LineWidth',2);
plot(...
Tx(1:floor(size(Tx,2)/10):size(Tx,2)),...
 U(1:floor(size(Tx,2)/10):size(Tx,2)),...
'ok');

ylabel('Control function $u_T(t)$','Interpreter','latex');
ylim([0.0, uMax]);
yticks([0, uMax/2, uMax]);

xlabel('Time $t$','Interpreter','latex');
xlim([0, totalYears]);
xticks([0, totalYears/2, totalYears]);

set(gca,'TickLabelInterpreter','latex');
box on; 

% 2.2: infected population
subplot(1,2,2);
hold on;
plot(Tx0,X0(:,3),':k','LineWidth',2);

plot(Tx,II,'-k','LineWidth',2);
plot(...
Tx(1:floor(size(Tx,2)/10):size(Tx,2)),...
 II(1:floor(size(Tx,2)/10):size(Tx,2)),...
'ok');

ylabel('Infected individuals $I(t)$','Interpreter','latex');
ylim([0., 1e6]);
yticks([0, 5e5, 1e6]);

xlabel('Time $t$','Interpreter','latex');
xlim([0, totalYears]);
xticks([0, totalYears/2, totalYears]);

set(gca,'TickLabelInterpreter','latex');
box on; 

% STEP 3: "naive" cost functional

disp('Cost of constant control:');
disp(trapz(Tx, cost*U.^2));
disp('Averted people:');
disp(X0(end,3) - II(end));


%-- density-dependent function --%

[TxPhIIOld, XPhIIOld, UPhIIOld, PhiShape3] = control2_feedback(pars, totalYears, alert, uMax);

TxPhII = linspace(TxPhIIOld(1), TxPhIIOld(end), 100);
XPhII  = interp1(TxPhIIOld, XPhIIOld, TxPhII);
UPhII  = interp1(TxPhIIOld, UPhIIOld, TxPhII);

II = XPhII(:,3)';

% STEP 2: plots

% 2.1: control Phi

subplot(1,2,1);
hold on;

plot(TxPhII, UPhII,'-g','LineWidth',2);
plot(...
TxPhII(1:floor(size(TxPhII,2)/10):size(TxPhII,2)),...
 UPhII(1:floor(size(TxPhII,2)/10):size(TxPhII,2)),...
'sg');

% 2.2: infected population plots

subplot(1,2,2);
hold on;

plot(TxPhII,XPhII(:,3),'-g','LineWidth',2);

plot(...
TxPhII(1:floor(size(TxPhII,2)/10):size(TxPhII,2)),...
    II(1:floor(size(TxPhII,2)/10):size(TxPhII,2)),...
'sg');



% STEP 3: "naive" cost functional

disp('Cost of feedback control:');
disp(trapz(TxPhII, cost*UPhII.^2));
disp('Averted people:');
disp(X0(end,3) - II(end));

%-- optimal control --%

[TxOld, XOld, u_opt3Old] = control3_optimal(pars, totalYears, cost, uMax);

Tx = linspace(TxOld(1), TxOld(end), 100);
X  = interp1(TxOld, XOld, Tx);

u_opt3 = interp1(TxOld, u_opt3Old, Tx);
II = X(:,3)';

% recovery functions plots
subplot(1,2,1);
hold on;

plot(Tx, u_opt3,'-b','LineWidth',2);
plot(Tx(1:floor(size(Tx,2)/10):size(Tx,2)),...
u_opt3(1:floor(size(Tx,2)/10):size(Tx,2)),...
'pb');

% infected population plots
subplot(1,2,2);
hold on;

plot(Tx,X(:,3),'-b','LineWidth',2);
plot(Tx(1:floor(size(Tx,2)/10):size(Tx,2)),...
II(1:floor(size(Tx,2)/10):size(Tx,2)),...
'pb');

% cost

disp('Cost of optimal control:');
disp(trapz(Tx, cost*u_opt3.^2));
disp('Averted people:');
disp(X0(end,3) - II(end));

h    = zeros(4, 1);
h(1) = plot(0,NaN,':k');
h(2) = plot(0,NaN,'o-k');
h(3) = plot(0,NaN,'s-g');
h(4) = plot(0,NaN,'p-b');

%axP2 = get(gca,'Position'); 

lh   = legend(h, ...
             {'No control','Constant','Feedback','Optimal'}, ...
              'Interpreter','latex', ...
              'Location','NorthEast', ...
              'Box','off');
          
%set(gcf,'units','points','position',[0,0,600,400]);

set(gcf, 'Renderer', 'painters', 'Position', [0 0 750 200]);
%set(gca, 'Position', axP2);

% %-- Phi Shape --%
% figure;
% hold on; 
% N=1e6;
% plot(0:100000:N,PhiShape1(1:10000:end),'ok');
% plot(0:100000:N,PhiShape2(1:10000:end),'sk');
% plot(0:100000:N,PhiShape3(1:10000:end),'pk');
% 
% plot(0:10:N,PhiShape1,'k','LineWidth',2);
% plot(0:10:N,PhiShape2,'k','LineWidth',2);
% plot(0:10:N,PhiShape3,'k','LineWidth',2);
% 
% plot([0 1],[0.1 0.1],'--','color',[0 0 0] + 0.5);
% plot([0 1],[0.4 0.4],'--','color',[0 0 0] + 0.5);
% plot([0 1],[0.8 0.8],'--','color',[0 0 0] + 0.5);
% 
% set(gca,'TickLabelInterpreter','latex');
% legend({ ...
%     '$I_{alert} = 8e5$',  ...
%     '$I_{alert} = 5e5$', ...
%     '$I_{alert} = 4e5$'}, 'Interpreter', 'latex', 'Location', 'NorthEastOutside');
% xlabel('Infection level $I$','Interpreter','latex');
% ylabel('Control function $\Phi(I)$','Interpreter','latex');
% % ylim([0,1]);
% box on;
% legend boxoff;
% 
% set(gcf,'units','points','position',[0,0,300,150]);
% set(gcf, 'PaperSize', [10 10]);
