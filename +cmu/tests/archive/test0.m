%% basic tests, unary operators

clear all
% unit.base_units
u = cmu.unit.units('SI');
% unit.base_units
u = cmu.unit.units('cm','s','gm','K','mol','coul');
cmu.unit.base_units
u = cmu.unit.units('s','cm','gm','K','mol','coul');
cmu.unit.base_units
u = cmu.unit.units('SI');

%%
u.m
+u.mm
-u.mm
a = u.cm
a.as(u.m)

%%
[u.m u.cm u.mm]
[u.m; u.cm; u.mm]
a=[[u.m u.cm];[u.nm u.um]]

a.as(u.cm)

%%
1*u.cm
u.m*0.5

%%
1/u.s
u.coul/5

%%
5*u.in/u.s

%%
a = 12*u.in
as(a,u.ft)
a.as(u.ft)

%%
as(6*u.in + 0.5*u.ft, u.ft)
as(6*u.in - 0.5*u.ft, u.ft)


