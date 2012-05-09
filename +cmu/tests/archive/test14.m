clear all
u = cmu.units;

%%GT
5*u.kg > 4*u.kg
5*u.kg > 4*u.ton
%5*u.kg > 4*u.m

[1 2 3] > [0 1 4]
[1 2 3]*u.m > [0 1 4]*u.m
%[1 2 3]*u.m > [0 1 4]*u.kg

%% LT
5*u.kg < 4*u.kg
5*u.kg < 4*u.ton

%% GE
5*u.kg >= 4*u.kg
5*u.kg >= 5*u.kg
%5*u.kg > 4*u.m

%% LE
5*u.kg <= 4*u.kg
5*u.kg <= 5*u.kg
%5*u.kg > 4*u.m