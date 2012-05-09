clear all; clc; close all
% Given this data, we want to fit a line to the data to extract the slope. 
u = cmu.units;
x = [0 0.5 1 1.5 2 3 4 6 10]*u.s;  % time
y = [0 -0.157 -0.315 -0.472 -0.629 -0.942 -1.255 -1.884 -3.147]*u.dimensionless; % dimensionless


X = [(x.^0)' x'];

M = X'*X;

assertAlmostEqual(M(2,2),168.50*u.s^2,1e-3)

z = X'*y';

b = (X'*X)\(X'*y')

b2 = M\z


X\y'