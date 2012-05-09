%% temperature function tests
% integration test
clear all
clc

u = cmu.unit.units;

%% 1
M = 10*u.mol;
Cp = 2.9e4*u.J/u.mol/u.K;

T1 = u.degC(60);
T2 = u.degC(30);

dH = M*Cp*(T2-T1);
as(dH, u.kJ,'%1.2g')

%% 2
T = linspace(T1,T2);
C = Cp*ones(size(T));

dH2 = M*trapz(T,C)
as(dH2,u.J)

%% 3
f = @(T) Cp*ones(size(T));
dH3 = M*quad(f,T1,T2)
as(dH3,u.J)
