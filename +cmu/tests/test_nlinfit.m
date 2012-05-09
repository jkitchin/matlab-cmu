clear all; close all; clc

u = cmu.units;
x = [0 1 2 3 4]'*u.s;
y = 1*u.mol/u.L + 4*u.mol/u.L/u.s*x + rand(size(x));

f = @(pars,x) pars(1) + pars(2)*x;

[pars residuals J] = nlinfit(x,y,f,[1.4*u.mol/u.L  3.4*u.mol/u.L/u.s]);

alpha = 0.05;
ci = nlparci(pars,residuals,'jacobian',J,'alpha',alpha)

[yp delta] = nlpredci(f,0.5*u.s,pars,residuals,'jacobian',J)