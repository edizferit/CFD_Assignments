% FLUID MECHANICS COMPUTATIONAL HOMEWORK QUESTION-2
clear all, close all, clc
format long

% Assumptions
% -Steady
% -Incompressible
% -Fully-developed
% -Non-laminar
% -Pipe flow

% PIPE PARAMETERS
epsilon = 0.00015; % surface roughness (m)
D = 0.4; % diameter (m)
L = 30; % length (m)

gravity = 9.81; % gravity (m/s^2)

% FLUID PARAMETERS 1
d = 995.7; % density (kg/m^3)
nu = 0.801*10^(-3); % dynamic viscosity (N*s/m)

% FLUID PARAMETERS 2
mflow = 15:1:30; % mass flow rate (kg/s)
Q = mflow/d; % volumetric flow rate (m^3/s)
v = Q/(pi/4*D^2); % velocity (m/s)
Re = d*v*D/nu; % Reynold's number 
f = zeros(size(mflow)); % friction factor
N = length(f);


%%
%========= COLEBROOK FORMULA =========%
cb = @(f,Re) f.^(-1/2) + 2*log10(epsilon/D/3.7 + 2.51./Re.*f^(-1/2));

% Secant method to find the roots 
delta_abs = 1e-8;
delta_rel = 1e-1;
maxI = 10000;
for j=1:N
    err = 100;
    relerr = 100;
    p0 = 0.1;
    p1 = 0.2;
    for i=1:maxI
        p2 = p1 - cb(p1,Re(j))*(p1-p0)/(cb(p1,Re(j))-cb(p0,Re(j)))/100;
        p0 = p1;
        p1 = p2;
        err = abs(p1 - p0);
        relerr = 2*abs(p1-p0)/(abs(p1)+abs(p0));
        if err < delta_abs && relerr < delta_rel
            break;
        end
    end
    f(j) = p2;
end

% Print the f values calculated from Colebrook formula
f

% Major head loss for Colebrook: Darcy-Weisbach equation
headLoss_cb = f.*(L/D).*v.^2/2/gravity;

% PLOT
figure
plot(mflow,headLoss_cb,'b','linewidth',1.5)
grid on
%ylim([2 10]*1e-4)
xlabel('Mass Flow Rate (kg/s)')
ylabel('Major Head Loss (m)')
title('Major Head Loss vs Mass Flow Rate','Colebrook Formula')

%%
%========== HAALAND ==========%
f_ha = (-1.8*log10((epsilon/D/3.7)^1.11 + 6.9./Re)).^(-2)

% Major head loss for Haaland: Darcy-Weisbach equation
headLoss_ha = f_ha.*(L/D).*v.^2/2/gravity;

% PLOT
figure
plot(mflow,headLoss_ha,'color',[0.9 0.4 0.17],'linewidth',1.5)
grid on
%ylim([2 10]*1e-4)
xlabel('Mass Flow Rate (kg/s)')
ylabel('Major Head Loss (m)')
title('Major Head Loss vs Mass Flow Rate','Haaland')

%%
%=========== PERCENT ERROR ===========%
error = abs(headLoss_cb-headLoss_ha);
percentError = error./headLoss_cb*100;

% PLOT
figure
plot(mflow,percentError,'k','linewidth',1.5)
grid on
xlabel('Mass Flow Rate (kg/s)')
ylabel('Percent Error (%)')
title('Mass Flow Rate vs Percent Error','wrt. Colebrook')

