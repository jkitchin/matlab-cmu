function c = constants(base)
% load the cmu.constants package for Matlab
%
% typical usage:
% >> c = cmu.constants;
% 
if ~(datenum(version('-date')) >= datenum('March 18, 2011'))
    error('You need at least Matlab R2011a to use cmu.units')
end

if nargin == 1
    c = cmu.unit.constants(base)
else
    c = cmu.unit.constants;
end