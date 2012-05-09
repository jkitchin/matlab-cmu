function u = conversion_factors(varargin)
% load the cmu.units package. Returns a struct containing the units.
%
% typical usage:
% >> u = cmu.units;
% >> velocity = 5*u.m/u.s % velocity in m/s
% 5*m/s
% >> t = 5*u.min % time in minutes
% 300*s
% >> d = velocity*t  %distance traveled
% 1500*m
% >> d.as(u.ft)
%
% ans =
%
% 4921.259843 ft

if ~(datenum(version('-date')) >= datenum('March 18, 2011'))
    error('You need at least Matlab R2011a to use cmu.units')
end

if nargin > 0
    u = cmu.unit.simple_units(varargin{:});
else
    u = cmu.unit.simple_units;
end