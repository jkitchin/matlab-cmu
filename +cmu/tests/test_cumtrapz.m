function test_suite = test_trapz
clear all
initTestSuite;

function test1
u = cmu.units;
y = [1 1 1 1]*u.s;
assertAlmostEqual(cumtrapz(y),[0 1 2 3]*u.s,1e-6)

function test2
u = cmu.units;
x = [0 1 2 3]*u.s;
y = [1 1 1 1]*u.m/u.s;
assertAlmostEqual(cumtrapz(x,y),[0 1 2 3]*u.m,1e-6)

function test3
u = cmu.units;
Y = [[0 1 2];[3 4 5]] * u.kg;
assertAlmostEqual(cumtrapz(Y,1),[[0 0 0];[1.5 2.5 3.5]]*u.kg,1e-6)
assertAlmostEqual(cumtrapz(Y,2),[[0 0.5 2];[0 3.5 8]]*u.kg,1e-6)
