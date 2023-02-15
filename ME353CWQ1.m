% FLUID MECHANICS COMPUTATIONAL HOMEWORK QUESTION-1
clear all, close all, clc

%=========== PARAMETERS ===========%
% VELOCITY (Column vector)
v = [7.07 10.00 12.25 14.14 15.81 17.32 18.71...
    20.00 21.21 22.36 23.45]'; % m/s
S = size(v); % size of data
N = length(v);
% MASS (Column vector)
m = [50 100 150 200 250 300 350 400 450 500 550]';
% DYNAMIC VISCOSITY
nu = 1.8*10^(-5); % Pa*s
% DENSITY
d = 1.225; % kg/m^3
% AREA
area = 25; % m^2

%=========== PI TERMS ===========%
pi_1 = m.*v/(nu*area); % 11x1
pi_2 = area^(3/2)*d./m; % 11x1


%=========== SECOND ORDER FIT ===========%
% y=ax^2+bx+c, y=pi_1, x=pi_2
A_sec = [pi_2.^2 pi_2 ones(size(pi_2))];
y_sec = pi_1;
x_sec = (A_sec'*A_sec)^(-1)*A_sec'*y_sec;
a_sec = x_sec(1);
b_sec = x_sec(2);
c_sec = x_sec(3);

% THE FUNCTION f(pi_2) = pi_1
pi_1_sec = @(pi_2) a_sec*pi_2.^2 + b_sec*pi_2 +c_sec;
t_sec = min(pi_2):0.01:max(pi_2);

% Pi-2 vs Pi-1 graph
figure
plot(pi_2,pi_1,'*r')                % As we see from the graph, second 
hold on, grid on                    % degree polynomial curve fitting  
plot(t_sec,pi_1_sec(t_sec'),'b')    % doesn't work well in our problem. 
xlabel('PI-2')                      % Later, the rms values also will
ylabel('PI-1')                      % lead us to the same conclusion.
title('Second Order Fit')           

rms_sec = sqrt(sum((pi_1_sec(pi_2) - pi_1).^2)/N) % RMS VALUE FOR SECOND D.


%%
%=========== POWER SERIES FIT ===========%
% y=a*x^b, y=pi_1, x=pi_2
A_pow = [log(pi_2) ones(size(pi_2))]; % lineerization done
y_pow = log(pi_1);
x_pow = (A_pow'*A_pow)^(-1)*A_pow'*y_pow;
b_pow = x_pow(1);
a_pow = exp(x_pow(2));

% THE FUNCTION f(pi_2) = pi_1
pi_1_pow = @(pi_2) a_pow*pi_2.^b_pow;
t_pow = min(pi_2):0.01:max(pi_2);

% Pi-2 vs Pi-1 graph
figure
plot(pi_2,pi_1,'*r')
hold on, grid on
plot(t_pow,pi_1_pow(t_pow'),'b')
xlabel('PI-2')
ylabel('PI-1')
title('Power Series Fit')

rms_pow = sqrt(sum((pi_1_pow(pi_2) - pi_1).^2)/N) % RMS VALUE FOR POWER S.

%%
%============== EXPONENTIAL FIT ==============%
% y=a*e^(b*x), y=pi_1, x=pi_2
A_exp = [pi_2 ones(size(pi_2))]; % lineerization done
y_exp = log(pi_1);
x_exp = (A_exp'*A_exp)^(-1)*A_exp'*y_exp;
b_exp = x_exp(1);
a_exp = exp(x_exp(2));

% THE FUNCTION f(pi_2) = pi_1
pi_1_exp = @(pi_2) a_exp*exp(b_exp*pi_2);
t_exp = min(pi_2):0.01:max(pi_2);

% Pi-2 vs Pi-1 graph
figure
plot(pi_2,pi_1,'*r')
hold on, grid on
plot(t_exp,pi_1_exp(t_exp'),'b')
xlabel('PI-2')
ylabel('PI-1')
title('Exponential Fit')

rms_exp = sqrt(sum((pi_1_exp(pi_2) - pi_1).^2)/N) % RMS VALUE FOR POWER S.



% COMMENT ON RMS VALUES:
% rms_sec = 4.2106e+06   Here, it can be seen that rms value of the
% rms_pow = 1.1315e+03   power series fit is much less than second degree 
% rms_exp = 4.7491e+06             polynomial fit and the exponential fit.
%                        
%                        So, we will continue with the function that we 
%                        found with power series fit for the rest of the 
%                        problem, since it gives more accurate results.


%%
% VELOCITY VS MASS graph
v_m = @(m) nu.*area./m.*a_pow.*(area.^(3/2).*d./m).^b_pow;
tm=min(m):0.01:max(m);

figure
plot(m,v,'*r')
hold on
plot(tm,v_m(tm),'b')
grid on
xlabel('Mass (kg)')
ylabel('Velocity (m/s)')
title('Velocity vs Mass')

% VELOCITY VS DENSITY graph
v_d = @(d) nu.*area./m.*a_pow.*(area.^(3/2).*d./m).^b_pow;
td=d-1:0.01:d+1;

figure
plot(td,v_d(td),'b')
grid on
xlabel('Density (kg/m^3)')
ylabel('Velocity (m/s)')
title('Velocity vs Density')
legend('Closest to origin: m=50 kg','Furthest from origin: m=550 kg')

% VELOCITY VS AREA graph
v_area = @(area) nu.*area./m.*a_pow.*(area.^(3/2).*d./m).^b_pow;
tarea=area-10:0.01:area+10;

figure
plot(tarea,v_area(tarea),'b')
grid on
xlabel('Area (m^2)')
ylabel('Velocity (m/s)')
title('Velocity vs Area')
legend('Closest to origin: m=50 kg','Furthest from origin: m=550 kg')


%%
% PART-B OF THE PROBLEM:
m_b = 10000; % (kg)
v_b = 30; % (m/s)

% Leave area term alone in the function that relates pi_1 to pi_2 which we 
% found by power series fitting technique.

area_b = (m_b^(b_pow+1)*v_b/nu/a_pow/d^b_pow)^(1/((3*b_pow+2)/2)); 
                                                    
% The result is area_b = 65.4986 m^2 which is a reasonable result if we 
% look from an conceptual perspective. On the other hand, we should note 
% that an extrapolation far beyond the data points we have is needed for 
% the given mass and velocity of the prototype. So, we can't be sure the 
% area calculated from the relationship would be enough to carry the
% prototype.


