% FLUID MECHANICS COMPUTATIONAL HOMEWORK QUESTION-3
clear all, close all, clc

% Laminar boundary layer
%
% U - free stream velocity in x direction
% v - kinematic viscosity
% d - density
% delta - boundary layer thickness
% deltaStar - boundary layer displacement thickness
% teta - boundary layer momentum thickness

% Symbolic Values
syms Re;
syms x;
syms U;
syms density;
syms b;         % width
syms nu;        % dynamic viscosity

% Table
blasius_table = [0 0 0 0.3321; 0.5 0.0415 0.1659 0.3309;...
1 0.1656 0.3298 0.323; 1.5 0.3701 0.4868 0.3026;...
2 0.65 0.6298 0.2668; 2.5 0.9963 0.7513 0.2174;...
3 1.3968 0.846 0.1614; 3.5 1.8377 0.913 0.1078;...
4 2.3057 0.9555 0.0642; 4.5 2.7901 0.9795 0.034;...
5 3.2833 0.9915 0.0159; 5.5 3.7806 0.9969 0.0066;...
6 4.2796 0.999 0.0024; 6.5 4.7793 0.9997 0.0008;...
7 5.2792 0.9999 0.0002; 7.5 5.7792 1 0.0001;...
8 6.2792 1 0];

% From table
eta = blasius_table(:,1);
f = blasius_table(:,2);
fprime = blasius_table(:,3);
fprime2 = blasius_table(:,4);
deta = 0.5; % step size
N = length(eta);


% Solve Differential Equation
tspan = [0 8];
y0 = [0 0 0.3321];
[t,y] = ode45(@(t,y) odecfn(y),tspan,y0);

% Plot Blasius solution
figure
plot(y(:,1),t,'linewidth',2)
hold on, grid on
plot(y(:,2),t,'linewidth',2)
plot(y(:,3),t,'linewidth',2)
xlabel('f, u/U, f"','Fontweight','bold')
ylabel('eta','Fontweight','bold')
title('Blasius Solution')
xlim([0 1])
legend('f','fprime','fprime2')


%%
%======= BOUNDARY LAYER THICKNESS =======%
% y = delta when u/U = 0.99
delta = 5/sqrt(Re)*x 


%%
%======= DISPLACEMENT THICKNESS =======%
% The displacement thickness is defined by:
% deltaStar = integral of (1-u/U)dy from 0 to infinity.
% deltaStar = integral of (1-u/U)dy from 0 to delta.
% Changing variables:
% f'=u/U, eta=y*sqrt(U/v/x), y=eta*sqrt(v*x/U), dy=deta*sqrt(v*x/U)
% deltaStar = x*sqrt(1/Re)*(integral of (1-f')deta from 0 to 5)

% Numeric Integral with Trapezoidal Rule for Displacement Thickness
T_d = 0;
for i=1:N-1
    T_d = T_d + deta/2*(1-fprime(i+1)+1-fprime(i));
end

T_d     % print the result of numeric integration

deltaStar = T_d*x/sqrt(Re)

%========== ERROR ANALYSIS ==========%
% Apply numeric differentiation to find f'''. 
% Maximum value of f''' will be used in error calculation.

fprime3 = zeros(size(fprime2));
for i=1:N-1
    fprime3(i+1) = (fprime2(i+1)-fprime2(i))/(eta(i+1)-eta(i));
end

% Estimated integration error is:
error_T_d = -(eta(end)-eta(1)).*max(abs(fprime3)).*deta.^2/12


%%
%======= MOMENTUM THICKNESS =======%
% The momentum thickness is defined by:
% teta = integral of u/U(1-u/U)dy from 0 to infinity.
% teta = integral of u/U(1-u/U)dy from 0 to delta.
% Changing variables:
% f'=u/U, eta=y*sqrt(U/v/x), y=eta*sqrt(v*x/U), dy=deta*sqrt(v*x/U)
% deltaStar = x*sqrt(1/Re)*(integral of f'(1-f')deta from 0 to 5)

% Numeric Integral with Trapezoidal Rule for Momentum Thickness
T_m = 0;
for i=1:N-1
    T_m = T_m + deta/2*(fprime(i+1)*(1-fprime(i+1))+...
        fprime(i)*(1-fprime(i)));
end

T_m     % print the result of numeric integration

teta = T_m*x/sqrt(Re)

%========== ERROR ANALYSIS ==========%
% Apply numeric differentiation to find second derivative of f'(1-f'). 
% Maximum value of this derivative will be used in error calculation.

f_error = zeros(size(fprime2));
for i=1:N-1
    f_error(i+1) = (fprime2(i+1)-2*fprime(i+1)*fprime2(i+1) - ...
        (fprime2(i)-2*fprime(i)*fprime2(i)))/(eta(i+1)-eta(i));
end

% Estimated integration error is:
error_T_m = -(eta(end)-eta(1)).*max(abs(f_error)).*deta.^2/12


%=============== SHEAR STRESS ================%
shear = nu*U*sqrt(Re)/x

%=============== DRAG FORCE ================%
dragForce_momentum = density*b*U^2*teta      % from momentum
dragForce_shear = 0.664*b*density*U^2*x/sqrt(Re)

%=============== LOCAL FRICTION COEFFICIENT ===============%
localFrictionCoefficient_shear = shear/(density*U^2/2)
localFrictionCoefficient_momentum = 0.6572*2/sqrt(Re)


%========== VELOCITY PROFILE and SHEAR PROFILE ========%
figure
plot(fprime,eta,'linewidth',2)
grid on, hold on
plot(fprime,eta*2,'linewidth',2)
plot(fprime,eta*3,'linewidth',2)
plot(fprime,eta*4,'linewidth',2)
xlabel('u/U ~ u')
ylabel('eta ~ y')
title('Velocity Profile in Boundary Layer')
legend('x1=x','x2=4x','x3=9x','x4=16x')

figure
plot(fprime2,eta,'linewidth',2)
grid on, hold on
plot(fprime2./2,eta,'linewidth',2)
plot(fprime2./3,eta,'linewidth',2)
plot(fprime2./4,eta,'linewidth',2)
xlabel('fprime2 ~ shear stress')
ylabel('eta ~ y')
title('Shear Stress Profile in Boundary Layer')
legend('x1=x','x2=4x','x3=9x','x4=16x')


%========== DIFFERENTIAL EQUATIONS ==========%
% This function, stores the 3 first order differential equations which will
% be used in ode45. These equations come from 2f''' + ff'' = 0.
function dydt = odecfn(y)
dydt = zeros(3,1);
dydt(1) = y(2);
dydt(2) = y(3);
dydt(3) = -y(1)*y(3)/2;
end
