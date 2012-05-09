function test_suite = test_creation

initTestSuite;

function test0
u = cmu.units;

% row vector
a = [u.m u.cm u.mm];
assertEqual(a,[1*u.m 0.01*u.m 0.001*u.m])

% column vector
a = [u.m; u.cm; u.mm];
assertEqual(a,[1*u.m; 0.01*u.m; 0.001*u.m])

% array
a=[[u.m u.cm];[u.nm u.um]];
assertEqual(a,[[1*u.m 0.01*u.m];[1e-9*u.m 1e-6*u.m]])

