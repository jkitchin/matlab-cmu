clear all;
u = cmu.unit.units;

a = 5*u.m
b = 6*u.s

try
    a + b
catch
    disp('error caught')
end
a*b

a > b

a < b

gt(a,b)