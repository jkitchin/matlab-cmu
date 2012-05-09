function test_suite = test_divide

initTestSuite;

function test0
u = cmu.units;

a = 5*u.m;

assert(a/5 == 1*u.m)
assert(5/a == 1/u.m)

function test1
u = cmu.units;

a = [10 5]*u.m;

b1 = a/5;
assert(b1(1) == 2*u.m)
assert(b1(2) == 1*u.m)

function test2
u = cmu.units;

a = [10 5]*u.m;
b = [5 5]*u.s;

c = a./b;
assert(c(1) == 2*u.m/u.s)
assert(c(2) == 1*u.m/u.s)

function test3
u = cmu.units;

a = [10*u.m 5*u.m];
b = [5*u.s 5*u.s];

c = a./b;
assert(c(1) == 2*u.m/u.s)
assert(c(2) == 1*u.m/u.s)

c = a.\b;
assert(abs(c(1)-0.5*u.s/u.m) < 1e-3*u.s/u.m)
assert(c(2) == 1*u.s/u.m)

function test4
u = cmu.units;

a = [10*u.m 5/u.s];
b = [5*u.s 5/u.m];

c = a./b;
assert(c(1) == 2*u.m/u.s)
assert(c(2) == 1*u.m/u.s)
