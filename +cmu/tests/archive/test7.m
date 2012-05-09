function test7
clear all
close all

u = cmu.unit.units;

Fao = 2.2*u.mol/u.min;
Vspan = [0 20]*u.L;

[V,Fa] = ode45(@pfr,Vspan,Fao);

plot(V/u.L,Fa/(u.mol/u.min))
xlabel('Volume (L)')
ylabel('Fa (mol/min)')

X = (Fao - Fa(end))/Fao

function dFadV = pfr(V,Fa)
%% plug flow reactor mole balance
u = cmu.unit.units;
k = 1.2/u.s;
vo = 12*u.L/u.s;
Ca = Fa/vo;
ra = -k*Ca;

dFadV = ra;