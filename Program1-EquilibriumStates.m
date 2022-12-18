% Program 1: Calculating Equilibrium States And Decision Whether It Is
% Locally Asymptotically Stable

clear all;
close all;

% Definition Of Problem Parameters

k01 = 1.00;
k02 = 1.00;
k03 = 1.00;
k04 = 1.00;
k05 = 1.00;
k06 = 1.00;
k07 = 1.00;
k08 = 1.00;
k09 = 1.00;
k10 = 1.00;
k11 = 1.00;
k12 = 1.00;
E   = min([k08/k07 k11/k10]);

% Step 1: Finding x_s By Newton-Type Approach

% Definition Of Funktion w(x) Whose Unique Zero We Want to Find

w   = @(x)   x*((k03*k07*k09*x)/((x+k04)*(k08-k07*x)) ...
             + (k05*k10*k12*x)/((x+k06)*(k11-k10*x))) ...
             - k02*x - k01;
dw  = @(x)   ((k03*k07*k09*x)/((x+k04)*(k08-k07*x)) ...
             + (k05*k10*k12*x)/((x+k06)*(k11-k10*x))) ...
             + x*((k03*k07*k09*(k07*x*x+k04*k08))/((x+k04)^2*(k08-k07*x)^2) ...
             + (k05*k10*k12*(k10*x*x+k06*k11))/((x+k06)^2*(k11-k10*x)^2)) ...
             - k02;
           
% Determination Of Starting Values

x0  = 0.5*E;
A   = w(x0);

while(A<0)
  x0 = x0+0.5*(E-x0);
  A  = w(x0);
endwhile

% Newton-Type Approach For w(x) = 0

N_max  = 10;
z      = zeros(3,N_max+1);
z(1,1) = x0;
z(2,1) = w(x0);

for k = 1:1:N_max
  z(1,k+1) = z(1,k) - w(z(1,k))/dw(z(1,k));
  z(2,k+1) = w(z(1,k+1));
  z(3,k+1) = abs(z(1,k+1)-z(1,k));
endfor

x_s = z(1,end);
y_s = k07*k09*x_s/(k08-k07*x_s);
z_s = k10*k12*x_s/(k11-k10*x_s);

% Step 2: Calculating All Three Eigenvalues Of The Jacobian Of Our
%         Dynamical System

% Definition Of Jacobian Of Our Dynamical System

J11 = @(x,y,z) k02 - k03*k04*y/((x+k04)^2) - k05*k06*z/((x+k06)^2);
J12 = @(x)     -k03*x/(x+k04);
J13 = @(x)     -k05*x/(x+k06);
J21 = k07;
J22 = @(y)     -k08*k09/((y+k09)^2);
J23 = 0;
J31 = k10;
J32 = 0;
J33 = @(z)     -k11*k12/((z+k12)^2);
J   = [J11(x_s,y_s,z_s) J12(x_s) J13(x_s); ...
       J21 J22(y_s) J23; ...
       J31 J32 J33(z_s)];
       
% Compute Real Parts Of Eigenvalues

Eig_RP = real(eig(J));

if(Eig_RP(1)<0 && Eig_RP(2)<0 && Eig_RP(3)<0)
  disp("All eigenvalues of the Jacobian have negative real parts. Hence, the dynamical system is asymptotically stable.")
else
  disp("At least one eigenvalue of the Jacobian is non-negative. Hence, the dynamical system is not asymptotically stable.")
endif