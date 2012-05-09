% true = 1
% false = 0
clear all
u = cmu.unit.units;
assert(~(u.m == u.s))

assert(u.m == u.m)

[1 3]*u.m == [10 300]*u.dm

[u.cm u.s] ~= [0.01*u.m u.s]

assert(all(([u.cm u.s] ~= [0.01*u.m u.min]) - [0 1])==0)

assert(u.m ~= u.cm)

assert(1 == u.m)