function test_suite = test_indexing

initTestSuite;

function test1
u = cmu.units;

a = [1 2 3 4 5]*u.s;

assertEqual(a(1),1*u.s)

assertEqual(a(2:4),[2 3 4]*u.s)

assertEqual(a(1:3:5), [1 4]*u.s)