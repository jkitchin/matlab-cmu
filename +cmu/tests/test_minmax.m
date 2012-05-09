function test_suite = test_minmax
clear all; clc
initTestSuite;

function test1
u = cmu.units;

a = [1 -1 2]*u.s;

amin = min(a);
assertEqual(amin,-1*u.s)

[amin,aind] = min(a);
assertEqual(aind,2)

function test2
u = cmu.units;

a = [1 -1 2]*u.s;

amax = max(a);
assertEqual(amax,2*u.s)

[amax,aind] = max(a);
assertEqual(aind,3)