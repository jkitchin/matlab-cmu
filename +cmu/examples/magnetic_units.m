clear all; close all; clc
u = cmu.units;

%% define new units
% these will eventually be defined as built in derived units

u.H = u.m^2*u.kg/(u.s^2*u.A^2);
u.H.displaystring = 'H';

u.ohm = u.V/u.A;
u.ohm.displaystring = 'ohm';

%% unit operations
a = 1*u.H;
a.as(u.H)

%% convert to some other units containing V and A
% if we divide the H by a V, we can see what is left over
a/u.V
% we could express a H as V*s^2/coul
a.as(u.V*u.s^2/u.coul)
%%
% That it is not too helpful. we
% recognize s/coul as  1/A. so we can see what is left over after
% that.
a/u.V/(1/u.A)
% now all that is left is a second
%%
% So, a Henry should be V*s/A
a.as(u.V*u.s/u.A)

%% another example. 
% Let's put an H in terms of an ohm. An ohm is V/A, so based on the
% example above, an H should be an ohm*s
a/u.ohm

a.as(u.ohm*u.s)

