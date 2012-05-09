u = cmu.unit.units;
a = 45*u.BTU/u.lb;

as(a,u.cal/u.kg)
as(a,u.J/u.kg)
as(a,u.kWh/u.kg)
as(a,u.ft*u.lbf/u.lb)
a.as(u.ft*u.lbf/u.lb)

%% example 22.5
W = 31.3*u.lb/u.s*32.2*u.ft/u.s^2*25*u.ft
as(W,u.kW)
as(W,u.hp)