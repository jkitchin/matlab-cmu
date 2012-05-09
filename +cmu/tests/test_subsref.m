function test_suite = test_subsref
clc; clear all
initTestSuite;


function test0

u = cmu.units;

a = 1*u.m;
assertEqual(a(1),1*u.m)
assertEqual(a(1,1),1*u.m)

function test1
u = cmu.units;

a = [1 2]*u.m;
assert(a(1)==1*u.m)
assert(a(2)==2*u.m)

assert(a(1,1)==1*u.m)
assert(a(1,2)==2*u.m)

function test2
u = cmu.units;

a = [1; 2]*u.m;
assert(a(1)==1*u.m)
assert(a(2)==2*u.m)

assert(a(1,1)==1*u.m)
assert(a(2,1)==2*u.m)

function test3
u = cmu.units;
a = [1*u.m 2*u.s];

assert(a(1)==1*u.m)
assert(a(2)==2*u.s)
assert(a(1,1)==1*u.m)
assert(a(1,2)==2*u.s)

function test4
u = cmu.units;
a = [1*u.m; 2*u.s];

assert(a(1)==1*u.m)
assert(a(2)==2*u.s)
assert(a(1,1)==1*u.m)
assert(a(2,1)==2*u.s)

function test5
u = cmu.units;
a = [[1*u.m 2*u.s];[u.kg u.coul]];

assert(a(1)==1*u.m)
assert(a(2)==u.kg)
assert(a(3)==2*u.s)
assert(a(4)==u.coul)

assert(a(2,2) == u.coul)

function test6
u = cmu.units;
x = [0 1 2 3 4]'*u.s;
X = [x.^0 x];

assert(X(1)==1)
assert(X(1,2)==0*u.s)
assert(X(10) == 4*u.s)

