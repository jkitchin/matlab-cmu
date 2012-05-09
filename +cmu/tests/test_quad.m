function test_suite = test_quad

initTestSuite;


function test1
% integration with units on function and limits
u = cmu.units;
v = 1*u.m/u.s;

f = @(t) v*ones(size(t));

t1 = 0*u.s; t2 = 10*u.s;
d = quad(f,t1,t2);
assertAlmostEqual(d,10*u.m,1e-6)

function test2
clear all
u = cmu.units;

X1 = 0*u.dimensionless; X2 = 0.5*u.dimensionless;

f = @(X) 1*u.m^3*ones(size(X));
V = quad(f,X1,X2);

assertAlmostEqual(V,0.5*u.m^3,1e-6)

