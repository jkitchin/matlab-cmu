clc; clear all; close all;

u = cmu.units;
Cao = 1*u.kmol/u.m^3;

%% raw data provided
t = [0 0.5 1.0 1.5 2.0 3.0 4.0 6.0 10.0]*u.s;
Cc = 2*u.mol/u.m^3/u.s*t;


%% getting the slope
P = polyfit(t,Cc,2)

polyder(P)

a = polyval(P,1*u.s)