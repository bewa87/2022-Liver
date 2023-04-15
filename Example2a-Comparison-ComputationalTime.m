% Step 0: Close All Windows and Clear Everything

clear all;
close all;

% Step 1: Definition of System Parameters and Initial Conditions

k01 = 1.0;
k02 = 1.0;
k03 = 3.0;
k04 = 3.0;
k05 = 1.0;
k06 = 1.0;
k07 = 1.0;
k08 = 1.0;
k09 = 1.0;
k10 = 1.0;
k11 = 1.0;
k12 = 1.0;

x_start = 0.5; 
y_start = 3.0; 
z_start = 3.0;

% Step 2: Definition of Time Intervals and Solution Vectors

T    = 10;
h_v  = [2 1 0.5 0.1 0.0001];

% Step 2.1: Different Time Vectors for Different Time Steps

t_5  = 0:h_v(5):T;
t_5  = (t_5)';

% Step 2.2: Solution Vectors for Different Time Steps

% Vectors For h(5)

x_5_nsfd    = zeros(length(t_5),1);
y_5_nsfd    = zeros(length(t_5),1);
z_5_nsfd    = zeros(length(t_5),1);
x_5_nsfd(1) = x_start;
y_5_nsfd(1) = y_start;
z_5_nsfd(1) = z_start;

x_5_ee      = zeros(length(t_5),1);
y_5_ee      = zeros(length(t_5),1);
z_5_ee      = zeros(length(t_5),1);
x_5_ee(1)   = x_start;
y_5_ee(1)   = y_start;
z_5_ee(1)   = z_start;

x_5_rk2     = zeros(length(t_5),1);
y_5_rk2     = zeros(length(t_5),1);
z_5_rk2     = zeros(length(t_5),1);
x_5_rk2(1)  = x_start;
y_5_rk2(1)  = y_start;
z_5_rk2(1)  = z_start;

% Step 3: Loops Over Times For All Three Algorithms

% Step 3.5: Loop For h = h_v(4)

h = h_v(5);

tic();
for j = 1:1:length(t_5)-1
  
  % Method 1: Loop For NSFD-Method
  
  A             = 1 + h*(k03*y_5_nsfd(j)/(x_5_nsfd(j)+k04) ...
                  + k05*z_5_nsfd(j)/(x_5_nsfd(j)+k06));
  B             = 1 + h*k08/(y_5_nsfd(j)+k09);
  C             = 1 + h*k11/(z_5_nsfd(j)+k12);
  x_5_nsfd(j+1) = x_5_nsfd(j)*(1+h*k02)/A + (h*k01)/A;
  y_5_nsfd(j+1) = y_5_nsfd(j)/B + h*k07*x_5_nsfd(j+1)/B;
  z_5_nsfd(j+1) = z_5_nsfd(j)/C + h*k10*x_5_nsfd(j+1)/C;

endfor
time_nsfd = toc()

tic();
for j = 1:1:length(t_5)-1
 
  % Method 2: Loop For Explicit Eulerian Method
  
  x_5_ee(j+1)   = x_5_ee(j) + h*(k01+k02*x_5_ee(j)-k03*(x_5_ee(j)*y_5_ee(j))/(x_5_ee(j)+k04)-k05*(x_5_ee(j)*z_5_ee(j))/(x_5_ee(j)+k06));
  y_5_ee(j+1)   = y_5_ee(j) + h*(k07*x_5_ee(j)-k08*(y_5_ee(j))/(y_5_ee(j)+k09));
  z_5_ee(j+1)   = z_5_ee(j) + h*(k10*x_5_ee(j)-k11*(z_5_ee(j))/(z_5_ee(j)+k12));

endfor
time_ee = toc()

tic();
for j = 1:1:length(t_5)-1
  
  % Method 3: Loop For Runge-Kutta Method Of Second Order
  
  t_h           = t_5(j) + 0.5*h;
  t_f           = t_5(j) + 1.0*h;
  x_h           = x_5_rk2(j) + 0.5*h*(k01+k02*x_5_rk2(j)-k03*(x_5_rk2(j)*y_5_rk2(j))/(x_5_rk2(j)+k04)-k05*(x_5_rk2(j)*z_5_rk2(j))/(x_5_rk2(j)+k06));
  y_h           = y_5_rk2(j) + 0.5*h*(k07*x_5_rk2(j)-k08*(y_5_rk2(j))/(y_5_rk2(j)+k09));
  z_h           = z_5_rk2(j) + 0.5*h*(k10*x_5_rk2(j)-k11*(z_5_rk2(j))/(z_5_rk2(j)+k12));
  x_5_rk2(j+1)  = x_5_rk2(j) + h*(k01+k02*x_h-k03*(x_h*y_h)/(x_h+k04)-k05*(x_h*z_h)/(x_h+k06));
  y_5_rk2(j+1)  = y_5_rk2(j) + h*(k07*x_h-k08*(y_h)/(y_h+k09));
  z_5_rk2(j+1)  = z_5_rk2(j) + h*(k10*x_h-k11*(z_h)/(z_h+k12));
  
endfor
time_rk2 = toc()