function test_suite = test_ode_solvers
clc
initTestSuite;

function test1
u = cmu.units;
k = 0.1/u.s;

f = @(t,Ca) -k*Ca;

tspan = [0 0.5 1]*u.s;
Cao = 1*u.mol/u.L;

[t,Ca] = ode45(f,tspan,Cao);
assertElementsAlmostEqual(Ca(end),Cao*exp(-k*1*u.s),'relative',1e-3) %ode solver tolerance is low


function test1a
u = cmu.units;
k = 0.1/u.s;

f = @(t,Ca) -k*Ca;

tspan = [0 0.5 1]*u.s;
Cao = 1*u.mol/u.L;

sol = ode45(f,tspan,Cao);
t = sol.x;
Ca = sol.y;

assertElementsAlmostEqual(t(1),0*u.s,'relative',1e-6)
assertElementsAlmostEqual(Ca(end),Cao*exp(-k*1*u.s),'relative',1e-3) %ode solver tolerance is low

ca = deval(sol,1*u.s);
assertElementsAlmostEqual(Ca(end),ca,'relative',1e-3)

function test2
u = cmu.units;
k = 0.1/u.s;

f = @(t,Ca) -k*Ca;

tspan = [0 0.5 1]*u.s;
Cao = 1*u.mol/u.m^3;

[t,Ca] = ode15s(f,tspan,Cao);
assertElementsAlmostEqual(Ca(end),Cao*exp(-k*1*u.s),'relative',1e-3)

function test2a
u = cmu.units;
k = 0.1/u.s;

f = @(t,Ca) -k*Ca;

tspan = [0 0.5 1]*u.s;
Cao = 1*u.mol/u.L;

sol = ode15s(f,tspan,Cao);
t = sol.x;
Ca = sol.y;

assertElementsAlmostEqual(t(1),0*u.s,'relative',1e-6)
assertElementsAlmostEqual(Ca(end),Cao*exp(-k*1*u.s),'relative',1e-3) %ode solver tolerance is low

ca = deval(sol,1*u.s);

assertElementsAlmostEqual(Ca(end),ca,'relative',1e-3)

