function test_suite = test_operators

initTestSuite;


function test0
u = cmu.units;

a = 1*u.m;
% unary operators
assertEqual(-a, -1*u.m)
assertEqual(+a, 1*u.m)

function test1
% transpose
u = cmu.units;
a = [1 1]*u.m;
assertEqual(a',[1;1]*u.m)
assertEqual(transpose(a),[1;1]*u.m)
assertEqual(ctranspose(a),[1;1]*u.m)

function test2
u = cmu.units;
a = [1 2 3]*u.coul;
assertEqual(a.^2, [1 4 9]*u.coul^2)

b = [1 4 9]*u.coul^2;
assertEqual(sqrt(b), [1 2 3]*u.coul)