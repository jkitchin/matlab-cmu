function test_suite = test_times

initTestSuite;

function test0

u = cmu.units;
a = [1 2 3]*u.m;
b = [1; 2; 3]/u.s;

a1 = a*b;
assertEqual(a1, 14*u.m/u.s)

function test1
assertExceptionThrown(@t1,'units:plus')

function t1
u = cmu.units;
a = [u.m u.s u.kg];
a*a'; % this should not be allowed because it implies m^2 + s^2 + kg^2

function test2
u = cmu.units;

a = [u.m u.s u.kg];
b = [1/u.s; u.m/u.s^2; u.m/u.s/u.kg];

assert(a*b == 3*u.m/u.s)



