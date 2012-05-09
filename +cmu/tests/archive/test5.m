%% CSTR solution using fsolve and fzero

function test5
clear all
u = cmu.unit.units;

Ca_guess = 0.09*u.mol/u.L;
cstr(Ca_guess)

C_exit=fsolve(@cstr,Ca_guess)
as(C_exit,u.mol/u.L)

fzero(@cstr,Ca_guess)


function z = cstr(Ca)
u = cmu.unit.units;
V = 10*u.gal;
vo = 0.25*u.gal/u.s;
Cao = 2*u.mol/u.L;
k = 2.2/u.min;

ra = -k*Ca;
%V = vo*(Cao - Ca)/(k*Ca)
z = V - vo*(Cao - Ca)/(-ra);