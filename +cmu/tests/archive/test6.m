%% CSTR solution using fsolve and fzero

function test5
clear all
u = cmu.unit.units;

Ca_guess = 0.09*u.mol/u.L;
Cb_guess = 0.2*u.mol/u.L;
C_exit=fsolve(@cstr,[Ca_guess, Cb_guess])


function z = cstr(C)
Ca = C(1);
Cb = C(2);
u = cmu.unit.units;
V = 10*u.gal;
vo = 0.25*u.gal/u.s;
Cao = 2*u.mol/u.L;
Cbo = 2*u.mol/u.L;
k = 2.2*u.L/u.mol/u.min;

ra = -k*Ca*Cb;
rb = ra;
%V = vo*(Cao - Ca)/(k*Ca)
z1 = V - vo*(Cao - Ca)/(-ra);
z2 = V - vo*(Cbo - Cb)/(-rb);
z = [z1; z2];