function v = version
v = '2.3';
% 2/15/2012 Added error check to force initial conditions to be column
% vectors for ode45 and ode15s because with mixed units there is a bug
% that prevents it from working.

% 11/19/2011 Fixed bug in odesolvers when there are mixed units in the
% dependent variables. added a few other odesolver placeholders. some minor
% cleanup of code.

% 10/25/2011 I fixed some bugs in subsref that were causing problems with
% multiplication of mixed row/column unit matrices.
% 10/07/2011 added support for regress, nlinfit, nlparci and nlpredci
% 10/06/2011 added some electromagnetic units, and the luminosity base
% unit.
% 9/27/2011 fixed arguments to the odesolvers to be consistent with
% matlab, and not use varargin. 
%
% 9/25/2011 fixed linsolve to work with mixed row/column units.
%
% 9/23/2011 added polyder, polyval and polyint
%
% 9/22/2011 bug fixes for unusual polyfit cases, e.g. with
% dimensionless numbers.
% 9/22/2011 added limited support for polyfit with units, fixed some
% sprintf bugs, modified display function.
%
% 9/13/2011 added the colors function
% 9/11/2011 added linspace, logspace, minor bug fixes to mtimes
%
% 9/10/2011 modified display function to better print arrays of units.
% Added some Matlab version checking to units. cmu.units does not work with
% Matlab 2009 because we use ~ for arguments to ignore. Changed order of
% arguments for derc, so it is derc(x,y).
%
% 9/9/2011 rewrote times function to allow mixed unit vectors.
% modified transpose and ctranspose to make sure the unit exponents
% are transposed. rewrote the subsref function to get mixed unit
% vectors correct.
%
% 9/6/2011 added trapz, deval functions. Added support for solution
% structures to ode45 and ode15s
%
% 9/5/2011 bug fix in cumtrapz functions. added tests directory.
%
% 9/4/2011 bug fix in fzero and fsolve when no output arguments are
% specified.
%
% 9/4/2011 added cumtrapz function

