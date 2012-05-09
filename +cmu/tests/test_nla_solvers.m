function test_suite = test_nla_solvers
clc

initTestSuite;


function test1

u = cmu.units;

k = 1*u.N/u.m;
E = @(x) 0.0 - 0.5*k*x^3;
p = fzero(E,0.1*u.m);
% i don't know why I cant say two zeros are nearly equal so I use the
% asserts here.
assert(abs(double(p)) <= 1e-6)


function test2

u = cmu.units;
k = 1*u.N/u.m;
E = @(x) 0.0 - 0.5*k*x^3;
options = optimset('Display','off','TolFun',1e-12);
p = fsolve(E,0.1*u.m,options);

assert(abs(double(p)) <= 0.05)