% Step 0: Close All Windows and Clear Everything

%clear all;
%close all;

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

% Known Equilibrium

x_s = 0.641065863794874;
y_s = 1.786026457590849;
z_s = 1.786026457590849;

% Step 2: Definition of Time Intervals and Solution Vectors

T    = 25;
h_v  = [2 1 0.5 0.25 0.125 0.0625];

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
t_6  = 0:h_v(6):T;
t_6  = (t_6)';

% Step 2.2: Solution Vectors for Different Time Steps

% Vectors For h(1)

x_1_nsfd    = zeros(length(t_1),1);
y_1_nsfd    = zeros(length(t_1),1);
z_1_nsfd    = zeros(length(t_1),1);
x_1_nsfd(1) = x_start;
y_1_nsfd(1) = y_start;
z_1_nsfd(1) = z_start;

% Vectors For h(2)

x_2_nsfd    = zeros(length(t_2),1);
y_2_nsfd    = zeros(length(t_2),1);
z_2_nsfd    = zeros(length(t_2),1);
x_2_nsfd(1) = x_start;
y_2_nsfd(1) = y_start;
z_2_nsfd(1) = z_start;

% Vectors For h(3)

x_3_nsfd    = zeros(length(t_3),1);
y_3_nsfd    = zeros(length(t_3),1);
z_3_nsfd    = zeros(length(t_3),1);
x_3_nsfd(1) = x_start;
y_3_nsfd(1) = y_start;
z_3_nsfd(1) = z_start;

% Vectors For h(4)

x_4_nsfd    = zeros(length(t_4),1);
y_4_nsfd    = zeros(length(t_4),1);
z_4_nsfd    = zeros(length(t_4),1);
x_4_nsfd(1) = x_start;
y_4_nsfd(1) = y_start;
z_4_nsfd(1) = z_start;

% Vectors For h(5)

x_5_nsfd    = zeros(length(t_5),1);
y_5_nsfd    = zeros(length(t_5),1);
z_5_nsfd    = zeros(length(t_5),1);
x_5_nsfd(1) = x_start;
y_5_nsfd(1) = y_start;
z_5_nsfd(1) = z_start;

% Vectors For h(6)

x_6_nsfd    = zeros(length(t_6),1);
y_6_nsfd    = zeros(length(t_6),1);
z_6_nsfd    = zeros(length(t_6),1);
x_6_nsfd(1) = x_start;
y_6_nsfd(1) = y_start;
z_6_nsfd(1) = z_start;

% Step 3: Loops Over Times For All Three Algorithms

% Step 3.1: Loop For h = h_v(1)

h = h_v(1);

for j = 1:1:length(t_1)-1
  
  A             = 1 + h*(k03*y_1_nsfd(j)/(x_1_nsfd(j)+k04) ...
                  + k05*z_1_nsfd(j)/(x_1_nsfd(j)+k06));
  B             = 1 + h*k08/(y_1_nsfd(j)+k09);
  C             = 1 + h*k11/(z_1_nsfd(j)+k12);
  x_1_nsfd(j+1) = x_1_nsfd(j)*(1+h*k02)/A + (h*k01)/A;
  y_1_nsfd(j+1) = y_1_nsfd(j)/B + h*k07*x_1_nsfd(j+1)/B;
  z_1_nsfd(j+1) = z_1_nsfd(j)/C + h*k10*x_1_nsfd(j+1)/C;
  
endfor

% Step 3.2: Loop For h = h_v(2)

h = h_v(2);

for j = 1:1:length(t_2)-1
  
  A             = 1 + h*(k03*y_2_nsfd(j)/(x_2_nsfd(j)+k04) ...
                  + k05*z_2_nsfd(j)/(x_2_nsfd(j)+k06));
  B             = 1 + h*k08/(y_2_nsfd(j)+k09);
  C             = 1 + h*k11/(z_2_nsfd(j)+k12);
  x_2_nsfd(j+1) = x_2_nsfd(j)*(1+h*k02)/A + (h*k01)/A;
  y_2_nsfd(j+1) = y_2_nsfd(j)/B + h*k07*x_2_nsfd(j+1)/B;
  z_2_nsfd(j+1) = z_2_nsfd(j)/C + h*k10*x_2_nsfd(j+1)/C;
  
endfor

% Step 3.3: Loop For h = h_v(3)

h = h_v(3);

for j = 1:1:length(t_3)-1
  
  A             = 1 + h*(k03*y_3_nsfd(j)/(x_3_nsfd(j)+k04) ...
                  + k05*z_3_nsfd(j)/(x_3_nsfd(j)+k06));
  B             = 1 + h*k08/(y_3_nsfd(j)+k09);
  C             = 1 + h*k11/(z_3_nsfd(j)+k12);
  x_3_nsfd(j+1) = x_3_nsfd(j)*(1+h*k02)/A + (h*k01)/A;
  y_3_nsfd(j+1) = y_3_nsfd(j)/B + h*k07*x_3_nsfd(j+1)/B;
  z_3_nsfd(j+1) = z_3_nsfd(j)/C + h*k10*x_3_nsfd(j+1)/C;
  
endfor

% Step 3.4: Loop For h = h_v(4)

h = h_v(4);

for j = 1:1:length(t_4)-1
  
  A             = 1 + h*(k03*y_4_nsfd(j)/(x_4_nsfd(j)+k04) ...
                  + k05*z_4_nsfd(j)/(x_4_nsfd(j)+k06));
  B             = 1 + h*k08/(y_4_nsfd(j)+k09);
  C             = 1 + h*k11/(z_4_nsfd(j)+k12);
  x_4_nsfd(j+1) = x_4_nsfd(j)*(1+h*k02)/A + (h*k01)/A;
  y_4_nsfd(j+1) = y_4_nsfd(j)/B + h*k07*x_4_nsfd(j+1)/B;
  z_4_nsfd(j+1) = z_4_nsfd(j)/C + h*k10*x_4_nsfd(j+1)/C;
  
endfor

% Step 3.5: Loop For h = h_v(4)

h = h_v(5);

for j = 1:1:length(t_5)-1
  
  A             = 1 + h*(k03*y_5_nsfd(j)/(x_5_nsfd(j)+k04) ...
                  + k05*z_5_nsfd(j)/(x_5_nsfd(j)+k06));
  B             = 1 + h*k08/(y_5_nsfd(j)+k09);
  C             = 1 + h*k11/(z_5_nsfd(j)+k12);
  x_5_nsfd(j+1) = x_5_nsfd(j)*(1+h*k02)/A + (h*k01)/A;
  y_5_nsfd(j+1) = y_5_nsfd(j)/B + h*k07*x_5_nsfd(j+1)/B;
  z_5_nsfd(j+1) = z_5_nsfd(j)/C + h*k10*x_5_nsfd(j+1)/C;
  
endfor

% Step 3.6: Loop For h = h_v(6)

h = h_v(6);

for j = 1:1:length(t_6)-1
  
  A             = 1 + h*(k03*y_6_nsfd(j)/(x_6_nsfd(j)+k04) ...
                  + k05*z_6_nsfd(j)/(x_6_nsfd(j)+k06));
  B             = 1 + h*k08/(y_6_nsfd(j)+k09);
  C             = 1 + h*k11/(z_6_nsfd(j)+k12);
  x_6_nsfd(j+1) = x_6_nsfd(j)*(1+h*k02)/A + (h*k01)/A;
  y_6_nsfd(j+1) = y_6_nsfd(j)/B + h*k07*x_6_nsfd(j+1)/B;
  z_6_nsfd(j+1) = z_6_nsfd(j)/C + h*k10*x_6_nsfd(j+1)/C;
  
endfor

% Step 4: Solution Plots

figure(1)
plot(t_6,x_6_nsfd,"linestyle","--","linewidth",2);
hold on
plot(t_6,y_6_nsfd,"linestyle","-.","linewidth",2);
hold on
plot(t_6,z_6_nsfd,"linestyle",":","linewidth",2);
legend("x(t)","y(t)","z(t)")
title("Example 3: Convergence To Equilibrium (h = 0.0625)")
xlabel("t")
ylabel("Solution Components")

% Step 5: Convergence Plot

x_end    = [abs(x_1_nsfd(end)-x_s) abs(x_2_nsfd(end)-x_s) abs(x_3_nsfd(end)-x_s) ...
            abs(x_4_nsfd(end)-x_s) abs(x_5_nsfd(end)-x_s) abs(x_6_nsfd(end)-x_s)];
y_end    = [abs(y_1_nsfd(end)-y_s) abs(y_2_nsfd(end)-y_s) abs(y_3_nsfd(end)-y_s) ...
            abs(y_4_nsfd(end)-y_s) abs(y_5_nsfd(end)-y_s) abs(y_6_nsfd(end)-y_s)];
z_end    = [abs(z_1_nsfd(end)-z_s) abs(z_2_nsfd(end)-z_s) abs(z_3_nsfd(end)-z_s) ...
            abs(z_4_nsfd(end)-z_s) abs(z_5_nsfd(end)-z_s) abs(z_6_nsfd(end)-z_s)];
comp_end = [x_end; y_end; z_end];

figure(2)
plot(log10(h_v),log10(x_end),"linestyle","--","linewidth",2);
hold on
plot(log10(h_v),log10(y_end),"linestyle","-.","linewidth",2);
hold on
plot(log10(h_v),log10(z_end),"linestyle",":","linewidth",2);
hold on
plot(log10(h_v),1*log10(h_v),"linestyle","-","linewidth",2);
hold off
legend("x","y","z","Theoretical Order 1","location","southeast")
title("Example 3: Double-Logarithmic Convergence Plot")
xlabel("log(h)")
ylabel("Logarithmic Distance Of Sol.-Components And Equilibrium")
