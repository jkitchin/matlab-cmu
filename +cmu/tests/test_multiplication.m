function test_suite = test_operators
clc
initTestSuite;

function test1
u = cmu.units;
assertEqual(5*u.s + 5*u.s, 10*u.s)

function test1a
u = cmu.units;
a = [1 1]*u.s;
assertEqual(a + a, [2 2]*u.s)

function test2
u = cmu.units;
assertEqual(5*u.s - 5*u.s, 0*u.s)

function test2a
u = cmu.units;
a = [1 1]*u.s;
assertEqual(a - a, [0 0]*u.s)

function test3
u = cmu.units;
assertEqual(5*u.s * 5*u.s, 25*u.s^2)

function test3a
u = cmu.units;
a = [1 1]*u.m;
b = [2; 2]/u.s;
assertEqual(a*b,4*u.m/u.s)

function test3b
u = cmu.units;
a = [1 1]*u.m;
b = [2 2]/u.s;
assertEqual(a.*b, [2 2]*u.m/u.s)

function test4
u = cmu.units;
assertEqual(5*u.s / (5*u.s), 1)

function test4a
%mrdivide
u = cmu.units;
a = [[1 0]; [0 1]]*u.m;
b = [3 4]*u.m/u.s;
x = b/a;
assertEqual(x, [3 4]/u.s)

function test4b
% ./
u = cmu.units;
a = [1 1]*u.m;
b = [2 2]*u.s;
assertEqual(a./b, [0.5 0.5]*u.m/u.s)