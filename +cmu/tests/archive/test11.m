clear all
u = cmu.unit.units;
a1 = [1 2 3] * u.m;
a2 = [1; 2; 3] * u.m;


a3 = [1*u.m 2*u.m 3*u.m]
a4 = [1*u.m; 2*u.m; 3*u.m]

b1 = a3*a4
b2 = a3.*a3
b3 = a3./a3
b4 = a3 / a3