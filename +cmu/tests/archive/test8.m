clear all
u = cmu.unit.units;
a = 5*u.m;
b = 5*u.s;

%assert(a ~= b)

sprintf('%f %1.2f',1,a)