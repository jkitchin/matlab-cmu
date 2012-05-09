clear all; clc

u = cmu.units;
x = [0 1 2 3 4]'*u.s;
X = [x.^0 x];

X(10)