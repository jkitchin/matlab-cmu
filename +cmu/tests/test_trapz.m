function test_suite = test_trapz

initTestSuite;

function test1
u = cmu.units;
y = [1 1 1 1]*u.s;
assertAlmostEqual(trapz(y),3*u.s,1e-6)

function test2
u = cmu.units;
t = [0 1 2 3]*u.s;
v = [1 1 1 1]*u.m/u.s;

d = trapz(t,v);
assertAlmostEqual(trapz(t,v),3*u.m,1e-6)

function test3
u = cmu.units;

Y = [0 1 2
    3 4 5]*u.s;

assertAlmostEqual(trapz(Y,1), [1.5000    2.5000    3.5000]*u.s, 1e-6)
assertAlmostEqual(trapz(Y,2), [2;8]*u.s, 1e-6)