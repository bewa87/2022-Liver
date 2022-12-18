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

T    = 4;
h_v  = [2 1 0.5 0.1 0.01];

% Step 2.1: Different Time Vectors for Different Time Steps

t_1  = 0:h_v(1):T;
t_1  = (t_1)';
t_2  = 0:h_v(2):T;
t_2  = (t_2)';
t_3  = 0:h_v(3):T;
t_3  = (t_3)';
t_4  = 0:h_v(4):T;
t_4  = (t_4)';
t_5  = 0:h_v(5):T;
t_5  = (t_5)';

% Step 2.2: Solution Vectors for Different Time Steps

% Vectors For h(1)

x_1_nsfd    = zeros(length(t_1),1);
y_1_nsfd    = zeros(length(t_1),1);
z_1_nsfd    = zeros(length(t_1),1);
x_1_nsfd(1) = x_start;
y_1_nsfd(1) = y_start;
z_1_nsfd(1) = z_start;

x_1_ee      = zeros(length(t_1),1);
y_1_ee      = zeros(length(t_1),1);
z_1_ee      = zeros(length(t_1),1);
x_1_ee(1)   = x_start;
y_1_ee(1)   = y_start;
z_1_ee(1)   = z_start;

x_1_rk2     = zeros(length(t_1),1);
y_1_rk2     = zeros(length(t_1),1);
z_1_rk2     = zeros(length(t_1),1);
x_1_rk2(1)   = x_start;
y_1_rk2(1)   = y_start;
z_1_rk2(1)   = z_start;

% Vectors For h(2)

x_2_nsfd    = zeros(length(t_2),1);
y_2_nsfd    = zeros(length(t_2),1);
z_2_nsfd    = zeros(length(t_2),1);
x_2_nsfd(1) = x_start;
y_2_nsfd(1) = y_start;
z_2_nsfd(1) = z_start;

x_2_ee      = zeros(length(t_2),1);
y_2_ee      = zeros(length(t_2),1);
z_2_ee      = zeros(length(t_2),1);
x_2_ee(1)   = x_start;
y_2_ee(1)   = y_start;
z_2_ee(1)   = z_start;

x_2_rk2     = zeros(length(t_2),1);
y_2_rk2     = zeros(length(t_2),1);
z_2_rk2     = zeros(length(t_2),1);
x_2_rk2(1)   = x_start;
y_2_rk2(1)   = y_start;
z_2_rk2(1)   = z_start;

% Vectors For h(3)

x_3_nsfd    = zeros(length(t_3),1);
y_3_nsfd    = zeros(length(t_3),1);
z_3_nsfd    = zeros(length(t_3),1);
x_3_nsfd(1) = x_start;
y_3_nsfd(1) = y_start;
z_3_nsfd(1) = z_start;

x_3_ee      = zeros(length(t_3),1);
y_3_ee      = zeros(length(t_3),1);
z_3_ee      = zeros(length(t_3),1);
x_3_ee(1)   = x_start;
y_3_ee(1)   = y_start;
z_3_ee(1)   = z_start;

x_3_rk2     = zeros(length(t_3),1);
y_3_rk2     = zeros(length(t_3),1);
z_3_rk2     = zeros(length(t_3),1);
x_3_rk2(1)   = x_start;
y_3_rk2(1)   = y_start;
z_3_rk2(1)   = z_start;

% Vectors For h(4)

x_4_nsfd    = zeros(length(t_4),1);
y_4_nsfd    = zeros(length(t_4),1);
z_4_nsfd    = zeros(length(t_4),1);
x_4_nsfd(1) = x_start;
y_4_nsfd(1) = y_start;
z_4_nsfd(1) = z_start;

x_4_ee      = zeros(length(t_4),1);
y_4_ee      = zeros(length(t_4),1);
z_4_ee      = zeros(length(t_4),1);
x_4_ee(1)   = x_start;
y_4_ee(1)   = y_start;
z_4_ee(1)   = z_start;

x_4_rk2     = zeros(length(t_4),1);
y_4_rk2     = zeros(length(t_4),1);
z_4_rk2     = zeros(length(t_4),1);
x_4_rk2(1)   = x_start;
y_4_rk2(1)   = y_start;
z_4_rk2(1)   = z_start;

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
x_5_rk2(1)   = x_start;
y_5_rk2(1)   = y_start;
z_5_rk2(1)   = z_start;

% Step 3: Loops Over Times For All Three Algorithms

% Step 3.1: Loop For h = h_v(1)

h = h_v(1);

for j = 1:1:length(t_1)-1
  
  % Method 1: Loop For NSFD-Method
  
  A             = 1 + h*(k03*y_1_nsfd(j)/(x_1_nsfd(j)+k04) ...
                  + k05*z_1_nsfd(j)/(x_1_nsfd(j)+k06));
  B             = 1 + h*k08/(y_1_nsfd(j)+k09);
  C             = 1 + h*k11/(z_1_nsfd(j)+k12);
  x_1_nsfd(j+1) = x_1_nsfd(j)*(1+h*k02)/A + (h*k01)/A;
  y_1_nsfd(j+1) = y_1_nsfd(j)/B + h*k07*x_1_nsfd(j+1)/B;
  z_1_nsfd(j+1) = z_1_nsfd(j)/C + h*k10*x_1_nsfd(j+1)/C;
  
  % Method 2: Loop For Explicit Eulerian Method
  
  x_1_ee(j+1)   = x_1_ee(j) + h*(k01+k02*x_1_ee(j)-k03*(x_1_ee(j)*y_1_ee(j))/(x_1_ee(j)+k04)-k05*(x_1_ee(j)*z_1_ee(j))/(x_1_ee(j)+k06));
  y_1_ee(j+1)   = y_1_ee(j) + h*(k07*x_1_ee(j)-k08*(y_1_ee(j))/(y_1_ee(j)+k09));
  z_1_ee(j+1)   = z_1_ee(j) + h*(k10*x_1_ee(j)-k11*(z_1_ee(j))/(z_1_ee(j)+k12));
  
  % Method 3: Loop For Runge-Kutta Method Of Second Order
  
  t_h           = t_1(j) + 0.5*h;
  t_f           = t_1(j) + 1.0*h;
  x_h           = x_1_rk2(j) + 0.5*h*(k01+k02*x_1_rk2(j)-k03*(x_1_rk2(j)*y_1_rk2(j))/(x_1_rk2(j)+k04)-k05*(x_1_rk2(j)*z_1_rk2(j))/(x_1_rk2(j)+k06));
  y_h           = y_1_rk2(j) + 0.5*h*(k07*x_1_rk2(j)-k08*(y_1_rk2(j))/(y_1_rk2(j)+k09));
  z_h           = z_1_rk2(j) + 0.5*h*(k10*x_1_rk2(j)-k11*(z_1_rk2(j))/(z_1_rk2(j)+k12));
  x_1_rk2(j+1)  = x_1_rk2(j) + h*(k01+k02*x_h-k03*(x_h*y_h)/(x_h+k04)-k05*(x_h*z_h)/(x_h+k06));
  y_1_rk2(j+1)  = y_1_rk2(j) + h*(k07*x_h-k08*(y_h)/(y_h+k09));
  z_1_rk2(j+1)  = z_1_rk2(j) + h*(k10*x_h-k11*(z_h)/(z_h+k12));
  
endfor

% Step 3.2: Loop For h = h_v(2)

h = h_v(2);

for j = 1:1:length(t_2)-1
  
  % Method 1: Loop For NSFD-Method
  
  A             = 1 + h*(k03*y_2_nsfd(j)/(x_2_nsfd(j)+k04) ...
                  + k05*z_2_nsfd(j)/(x_2_nsfd(j)+k06));
  B             = 1 + h*k08/(y_2_nsfd(j)+k09);
  C             = 1 + h*k11/(z_2_nsfd(j)+k12);
  x_2_nsfd(j+1) = x_2_nsfd(j)*(1+h*k02)/A + (h*k01)/A;
  y_2_nsfd(j+1) = y_2_nsfd(j)/B + h*k07*x_2_nsfd(j+1)/B;
  z_2_nsfd(j+1) = z_2_nsfd(j)/C + h*k10*x_2_nsfd(j+1)/C;
  
  % Method 2: Loop For Explicit Eulerian Method
  
  x_2_ee(j+1)   = x_2_ee(j) + h*(k01+k02*x_2_ee(j)-k03*(x_2_ee(j)*y_2_ee(j))/(x_2_ee(j)+k04)-k05*(x_2_ee(j)*z_2_ee(j))/(x_2_ee(j)+k06));
  y_2_ee(j+1)   = y_2_ee(j) + h*(k07*x_2_ee(j)-k08*(y_2_ee(j))/(y_2_ee(j)+k09));
  z_2_ee(j+1)   = z_2_ee(j) + h*(k10*x_2_ee(j)-k11*(z_2_ee(j))/(z_2_ee(j)+k12));
  
  % Method 3: Loop For Runge-Kutta Method Of Second Order
  
  t_h           = t_2(j) + 0.5*h;
  t_f           = t_2(j) + 1.0*h;
  x_h           = x_2_rk2(j) + 0.5*h*(k01+k02*x_2_rk2(j)-k03*(x_2_rk2(j)*y_2_rk2(j))/(x_2_rk2(j)+k04)-k05*(x_2_rk2(j)*z_2_rk2(j))/(x_2_rk2(j)+k06));
  y_h           = y_2_rk2(j) + 0.5*h*(k07*x_2_rk2(j)-k08*(y_2_rk2(j))/(y_2_rk2(j)+k09));
  z_h           = z_2_rk2(j) + 0.5*h*(k10*x_2_rk2(j)-k11*(z_2_rk2(j))/(z_2_rk2(j)+k12));
  x_2_rk2(j+1)  = x_2_rk2(j) + h*(k01+k02*x_h-k03*(x_h*y_h)/(x_h+k04)-k05*(x_h*z_h)/(x_h+k06));
  y_2_rk2(j+1)  = y_2_rk2(j) + h*(k07*x_h-k08*(y_h)/(y_h+k09));
  z_2_rk2(j+1)  = z_2_rk2(j) + h*(k10*x_h-k11*(z_h)/(z_h+k12));
  
endfor

% Step 3.3: Loop For h = h_v(3)

h = h_v(3);

for j = 1:1:length(t_3)-1
  
  % Method 1: Loop For NSFD-Method
  
  A             = 1 + h*(k03*y_3_nsfd(j)/(x_3_nsfd(j)+k04) ...
                  + k05*z_3_nsfd(j)/(x_3_nsfd(j)+k06));
  B             = 1 + h*k08/(y_3_nsfd(j)+k09);
  C             = 1 + h*k11/(z_3_nsfd(j)+k12);
  x_3_nsfd(j+1) = x_3_nsfd(j)*(1+h*k02)/A + (h*k01)/A;
  y_3_nsfd(j+1) = y_3_nsfd(j)/B + h*k07*x_3_nsfd(j+1)/B;
  z_3_nsfd(j+1) = z_3_nsfd(j)/C + h*k10*x_3_nsfd(j+1)/C;
  
  % Method 2: Loop For Explicit Eulerian Method
  
  x_3_ee(j+1)   = x_3_ee(j) + h*(k01+k02*x_3_ee(j)-k03*(x_3_ee(j)*y_3_ee(j))/(x_3_ee(j)+k04)-k05*(x_3_ee(j)*z_3_ee(j))/(x_3_ee(j)+k06));
  y_3_ee(j+1)   = y_3_ee(j) + h*(k07*x_3_ee(j)-k08*(y_3_ee(j))/(y_3_ee(j)+k09));
  z_3_ee(j+1)   = z_3_ee(j) + h*(k10*x_3_ee(j)-k11*(z_3_ee(j))/(z_3_ee(j)+k12));
  
  % Method 3: Loop For Runge-Kutta Method Of Second Order
  
  t_h           = t_3(j) + 0.5*h;
  t_f           = t_3(j) + 1.0*h;
  x_h           = x_3_rk2(j) + 0.5*h*(k01+k02*x_3_rk2(j)-k03*(x_3_rk2(j)*y_3_rk2(j))/(x_3_rk2(j)+k04)-k05*(x_3_rk2(j)*z_3_rk2(j))/(x_3_rk2(j)+k06));
  y_h           = y_3_rk2(j) + 0.5*h*(k07*x_3_rk2(j)-k08*(y_3_rk2(j))/(y_3_rk2(j)+k09));
  z_h           = z_3_rk2(j) + 0.5*h*(k10*x_3_rk2(j)-k11*(z_3_rk2(j))/(z_3_rk2(j)+k12));
  x_3_rk2(j+1)  = x_3_rk2(j) + h*(k01+k02*x_h-k03*(x_h*y_h)/(x_h+k04)-k05*(x_h*z_h)/(x_h+k06));
  y_3_rk2(j+1)  = y_3_rk2(j) + h*(k07*x_h-k08*(y_h)/(y_h+k09));
  z_3_rk2(j+1)  = z_3_rk2(j) + h*(k10*x_h-k11*(z_h)/(z_h+k12));
  
endfor

% Step 3.4: Loop For h = h_v(4)

h = h_v(4);

for j = 1:1:length(t_4)-1
  
  % Method 1: Loop For NSFD-Method
  
  A             = 1 + h*(k03*y_4_nsfd(j)/(x_4_nsfd(j)+k04) ...
                  + k05*z_4_nsfd(j)/(x_4_nsfd(j)+k06));
  B             = 1 + h*k08/(y_4_nsfd(j)+k09);
  C             = 1 + h*k11/(z_4_nsfd(j)+k12);
  x_4_nsfd(j+1) = x_4_nsfd(j)*(1+h*k02)/A + (h*k01)/A;
  y_4_nsfd(j+1) = y_4_nsfd(j)/B + h*k07*x_4_nsfd(j+1)/B;
  z_4_nsfd(j+1) = z_4_nsfd(j)/C + h*k10*x_4_nsfd(j+1)/C;
  
  % Method 2: Loop For Explicit Eulerian Method
  
  x_4_ee(j+1)   = x_4_ee(j) + h*(k01+k02*x_4_ee(j)-k03*(x_4_ee(j)*y_4_ee(j))/(x_4_ee(j)+k04)-k05*(x_4_ee(j)*z_4_ee(j))/(x_4_ee(j)+k06));
  y_4_ee(j+1)   = y_4_ee(j) + h*(k07*x_4_ee(j)-k08*(y_4_ee(j))/(y_4_ee(j)+k09));
  z_4_ee(j+1)   = z_4_ee(j) + h*(k10*x_4_ee(j)-k11*(z_4_ee(j))/(z_4_ee(j)+k12));
  
  % Method 3: Loop For Runge-Kutta Method Of Second Order
  
  t_h           = t_4(j) + 0.5*h;
  t_f           = t_4(j) + 1.0*h;
  x_h           = x_4_rk2(j) + 0.5*h*(k01+k02*x_4_rk2(j)-k03*(x_4_rk2(j)*y_4_rk2(j))/(x_4_rk2(j)+k04)-k05*(x_4_rk2(j)*z_4_rk2(j))/(x_4_rk2(j)+k06));
  y_h           = y_4_rk2(j) + 0.5*h*(k07*x_4_rk2(j)-k08*(y_4_rk2(j))/(y_4_rk2(j)+k09));
  z_h           = z_4_rk2(j) + 0.5*h*(k10*x_4_rk2(j)-k11*(z_4_rk2(j))/(z_4_rk2(j)+k12));
  x_4_rk2(j+1)  = x_4_rk2(j) + h*(k01+k02*x_h-k03*(x_h*y_h)/(x_h+k04)-k05*(x_h*z_h)/(x_h+k06));
  y_4_rk2(j+1)  = y_4_rk2(j) + h*(k07*x_h-k08*(y_h)/(y_h+k09));
  z_4_rk2(j+1)  = z_4_rk2(j) + h*(k10*x_h-k11*(z_h)/(z_h+k12));
  
endfor

% Step 3.5: Loop For h = h_v(4)

h = h_v(5);

for j = 1:1:length(t_5)-1
  
  % Method 1: Loop For NSFD-Method
  
  A             = 1 + h*(k03*y_5_nsfd(j)/(x_5_nsfd(j)+k04) ...
                  + k05*z_5_nsfd(j)/(x_5_nsfd(j)+k06));
  B             = 1 + h*k08/(y_5_nsfd(j)+k09);
  C             = 1 + h*k11/(z_5_nsfd(j)+k12);
  x_5_nsfd(j+1) = x_5_nsfd(j)*(1+h*k02)/A + (h*k01)/A;
  y_5_nsfd(j+1) = y_5_nsfd(j)/B + h*k07*x_5_nsfd(j+1)/B;
  z_5_nsfd(j+1) = z_5_nsfd(j)/C + h*k10*x_5_nsfd(j+1)/C;
  
  % Method 2: Loop For Explicit Eulerian Method
  
  x_5_ee(j+1)   = x_5_ee(j) + h*(k01+k02*x_5_ee(j)-k03*(x_5_ee(j)*y_5_ee(j))/(x_5_ee(j)+k04)-k05*(x_5_ee(j)*z_5_ee(j))/(x_5_ee(j)+k06));
  y_5_ee(j+1)   = y_5_ee(j) + h*(k07*x_5_ee(j)-k08*(y_5_ee(j))/(y_5_ee(j)+k09));
  z_5_ee(j+1)   = z_5_ee(j) + h*(k10*x_5_ee(j)-k11*(z_5_ee(j))/(z_5_ee(j)+k12));
  
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

% Step 4: Plots For Comparison

##figure(1)
##plot(t_1,x_1_nsfd)
##hold on
##plot(t_1,x_1_ee)
##hold on
##plot(t_1,x_1_rk2)
##title("Example 2: Comparison Of Three Numerical Algorithms (h = 1)")
##xlabel("t")
##ylabel("x(t)")
##legend("(M1)", "(M2)", "(M3)")
##
##figure(2)
##plot(t_2,x_2_nsfd)
##hold on
##plot(t_2,x_2_ee)
##hold on
##plot(t_2,x_2_rk2)
##title("Example 2: Comparison Of Three Numerical Algorithms (h = 0.75)")
##xlabel("t")
##ylabel("x(t)")
##legend("(M1)", "(M2)", "(M3)")
##
##figure(3)
##plot(t_3,x_3_nsfd)
##hold on
##plot(t_3,x_3_ee)
##hold on
##plot(t_3,x_3_rk2)
##title("Example 2: Comparison Of Three Numerical Algorithms (h = 0.5)")
##xlabel("t")
##ylabel("x(t)")
##legend("(M1)", "(M2)", "(M3)")
##
##figure(4)
##plot(t_4,x_4_nsfd)
##hold on
##plot(t_4,x_4_ee)
##hold on
##plot(t_4,x_4_rk2)
##title("Example 2: Comparison Of Three Numerical Algorithms (h = 0.1)")
##xlabel("t")
##ylabel("x(t)")
##legend("(M1)", "(M2)", "(M3)")

figure(1)
subplot(2,2,1);
plot(t_1,x_1_nsfd)
hold on
plot(t_1,x_1_ee)
hold on
plot(t_1,x_1_rk2)
title("Example 2: Time Step h = 2")
xlabel("t")
ylabel("x(t)")
%legend({"(M1)", "(M2)", "(M3)"},"location","northwest")
hold off
subplot(2,2,2);
plot(t_2,x_2_nsfd)
hold on
plot(t_2,x_2_ee)
hold on
plot(t_2,x_2_rk2)
title("Example 2: Time Step h = 1")
xlabel("t")
ylabel("x(t)")
%legend({"(M1)", "(M2)", "(M3)"},"location","northwest")
hold off
subplot(2,2,3);
plot(t_3,x_3_nsfd)
hold on
plot(t_3,x_3_ee)
hold on
plot(t_3,x_3_rk2)
title("Example 2: Time Step h = 0.5")
xlabel("t")
ylabel("x(t)")
%legend({"(M1)", "(M2)", "(M3)"},"location","northwest")
hold off
subplot(2,2,4);
plot(t_4,x_4_nsfd)
hold on
plot(t_4,x_4_ee)
hold on
plot(t_4,x_4_rk2)
title("Example 2: Time Step h = 0.1")
xlabel("t")
ylabel("x(t)")
%legend({"(M1)", "(M2)", "(M3)"},"location","northwest")

figure(2)
plot(t_5,x_5_nsfd)
hold on
plot(t_5,x_5_ee)
hold on
plot(t_5,x_5_rk2)
title("Example 2: Comparison Of Three Numerical Algorithms (h = 0.01)")
xlabel("t")
ylabel("x(t)")
legend("(M1)", "(M2)", "(M3)")