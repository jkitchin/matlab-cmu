function test_suite = test_comparisons
clc
initTestSuite;

function test1
u = cmu.units;

assertTrue(1*u.s == 1*u.s)
assertFalse(1*u.s == 1*u.kg)
assertTrue(1*u.s ~= 1*u.kg)

function test2
u = cmu.units;

assertTrue(1*u.hr > 1*u.s)
assertFalse(1*u.hr < 1*u.s)

function test3
u = cmu.units;

assertTrue(1*u.hr >= 1*u.s)
assertTrue(1*u.s >= 1*u.s)
assertTrue(1*u.s <= 1*u.s)

assertFalse(1*u.hr <= 1*u.s)