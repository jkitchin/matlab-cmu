function test_suite = test_linsolve
clear all; clc
initTestSuite;

function test1
u = cmu.units;
a = [[1 0]; [0 1]]*u.m;
b = [3; 4]*u.m/u.s;
x = a\b;
assertEqual(x, [3; 4]/u.s)

function test2
u = cmu.units;
a = [[1 0]; [0 1]]*u.m;
b = [3; 4]*u.m/u.s;
x2 = linsolve(a,b);
assertEqual(x2, [3; 4]/u.s)