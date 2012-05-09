function test_as
u=cmu.units;
a = 60*u.s;
assertEqual(a.as(u.min),'1.000*min')