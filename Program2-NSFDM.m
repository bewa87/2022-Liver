% Step 0: Close All Windows and Clear Everything

clear all;
close all;

% Step 1: Definition of System Parameters and Initial Conditions

k01 = 1.0;
k02 = 1.5;
k03 = 1.0;
k04 = 1.0;
k05 = 1.0;
k06 = 1.0;
k07 = 1.0;
k08 = 1.0;
k09 = 1.0;
k10 = 1.0;
k11 = 1.0;
k12 = 1.0;

x_start = 0.1; %0.1;
y_start = 0.0; %0.0;
z_start = 0.0; %0.0;

% Step 2: Definition of Time Intervals and Solution Vectors

T    = 1000;
h    = 0.001;
t    = 0:h:T;
t    = t';
x    = zeros(length(t),1);
y    = zeros(length(t),1);
z    = zeros(length(t),1);
x(1) = x_start;
y(1) = y_start;
z(1) = z_start;

% Step 3: Non-Standard Finite Difference Numerical Solution Approximation

for j = 1:1:length(t)-1
  A      = 1 + h*(k03*y(j)/(x(j)+k04) + k05*z(j)/(x(j)+k06));
  B      = 1 + h*k08/(y(j)+k09);
  C      = 1 + h*k11/(z(j)+k12);
  x(j+1) = x(j)*(1+h*k02)/A + (h*k01)/A;
  y(j+1) = y(j)/B + h*k07*x(j+1)/B;
  z(j+1) = z(j)/C + h*k10*x(j+1)/C;
endfor

figure(1)
plot(t(1:1:end),x(1:1:end))

figure(2)
plot(t(1:1:end),y(1:1:end))

figure(3)
plot(t(1:1:end),z(1:1:end))

figure(4)
plot(t(1:100:end),x(1:100:end))

figure(5)
plot(t(1:100:end),y(1:100:end))

figure(6)
plot(t(1:100:end),z(1:100:end))

figure(7)