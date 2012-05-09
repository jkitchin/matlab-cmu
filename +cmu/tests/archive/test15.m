clear all
u = cmu.units;

a = [ 1 3]*u.J;
b = [3 3]*u.s;

c = a./b;

c(1).as(u.J./u.s)

a = [ 1*u.J 3/u.s];
b = [3*u.s 3/u.J];

c = a./b;

c(1).as(u.J./u.s)