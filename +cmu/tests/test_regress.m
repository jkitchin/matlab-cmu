clear all; close all; clc

u = cmu.units('m','kg','s','kmol','coul','cd','K');
% 
% x = [0 1 2 3 4]'*u.s;
% y = 1*u.kmol/u.m^3 + 4*u.kmol/u.m^3/u.s*x + rand(size(x));
% 
% F = [x.^0 x];
% 
% [b bint r rint] = regress(y,F)

%%
%F\y

%%
clc
C = [0.2 0.02 0.01 0.005 0.002]'*u.kmol/u.m^3;
r = -[1.08 0.55 0.38 0.2 0.09]'*u.kmol/u.m^3/u.s;

z1 = C./-r

x = 1./C;
y = 1./-r;

plot(x,y)

z2 = y./x

z3 = [x.^0 x]\y

F = [x.^0 x]

F\y

b = linsolve(F,y)
intercept = b(1)
vmax = 1/intercept

slope = b(2)
Km = vmax*slope