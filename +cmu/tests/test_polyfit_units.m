function test_suite = test_polyfit_units
clc; clear all; close all;
initTestSuite;


function test1
clear all; close all;
TOLERANCE = 1e-6;
u = cmu.units;
%% raw data provided
t = [0 0.5 1.0 1.5 2.0 3.0 4.0 6.0 10.0]*u.s;
Cc = 2*u.mol/u.m^3/u.s*t;

%% getting the slope
Pd = polyfit(double(t),double(Cc),1);

P = polyfit(t,Cc,1);
assert(abs(double(P(1) - 2*u.mol/u.m^3/u.s)) <= TOLERANCE);
assert(abs(double(P(2) - 0*u.mol/u.m^3)) <= TOLERANCE);

function test2
clear all; close all;
TOLERANCE = 1e-6;
u = cmu.units;

%% raw data provided
t = [0 0.5 1.0 1.5 2.0 3.0 4.0 6.0 10.0]; % no dimension
Cc = 2*u.mol/u.m^3*t;

%% getting the slope
Pd = polyfit(double(t),double(Cc),1);

P = polyfit(t,Cc,1);
assert(abs(double(P(1) - 2*u.mol/u.m^3)) <= TOLERANCE);
assert(abs(double(P(2) - 0*u.mol/u.m^3)) <= TOLERANCE);

function test3
clear all; close all;
TOLERANCE = 1e-6;
u = cmu.units;

%% raw data provided
t = [0 0.5 1.0 1.5 2.0 3.0 4.0 6.0 10.0]*u.dimensionless; % no dimension
Cc = 2*u.mol/u.m^3*t;

%% getting the slope
Pd = polyfit(double(t),double(Cc),1);

P = polyfit(t,Cc,1);
assert(abs(double(P(1) - 2*u.mol/u.m^3)) <= TOLERANCE);
assert(abs(double(P(2) - 0*u.mol/u.m^3)) <= TOLERANCE);

function test4
clear all; close all;
TOLERANCE = 1e-6;
u = cmu.units;

%% raw data provided
t = [0 0.5 1.0 1.5 2.0 3.0 4.0 6.0 10.0]; 
Cc = 2*t;

%% getting the slope
P = polyfit(t*u.s,Cc,1);
assert(abs(double(P(1) - 2/u.s)) <= TOLERANCE);
assert(abs(double(P(2) - 0)) <= TOLERANCE);