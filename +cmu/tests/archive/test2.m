%% ode integration

clear all
close all
clc
u = cmu.unit.units;

tspan = [0*u.s 0.15*u.s 0.3*u.s];
init = 5*u.kg;

ode = @(t,M) 0.5/M;

[t,M] = ode45(ode,tspan,init);

plot(t,M)

%% stiff test
tspan = [0 0.3]*u.s;
init = 5*u.kg;

ode = @(t,M) 0.5*M;

[t,M] = ode15s(ode,tspan,init);

plot(t,M)
