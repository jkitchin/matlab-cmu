clc; clear all
u = cmu.units;

Cc =2 *(u.mol/u.L);

Cc.display

Cc.display('%1.10f')


%%
s = sprintf('%0.2f' , Cc)
assert(strcmp(s,'2000.00/m^3*mol'))

s2 = sprintf('%s', Cc.as(u.mol/u.L));
 
assert(strcmp(s2,'2.000*mol/L'))

s3 = sprintf('%s', Cc.as(u.mol/u.L,'%0.2f'));
assert(strcmp(s3,'2.00*mol/L'))

%%
k = 1/u.s
s4 = k.as(1/u.min)
assert(strcmp(s4,'60.000*min^-1'))