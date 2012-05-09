clear all

u = cmu.units;
Wa = [0 1 2 3]*u.m
Ca = [0 1 2 3]*u.s

Y = [Wa' Ca']

a = Y(:,1)