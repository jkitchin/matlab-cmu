matlab-cmu
==========

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.


Repository for the +cmu matlab package. this package contains a units class that enforces unit algebra for engineering calculations in Matlab. Some features are:

1. Automatic unit conversion to a set of user-defined units (default is SI).
2. Automatic unit algebra, e.g. kg*kg = kg^2, kg + m is not allowed.
3. Many overloaded Matlab functions so that units work in ode solvers, linear algebra, nonlinear algebra, etc...
4. Defined functions where units are not allowed, e.g. trigonemetric functions, log, exp...
5. Small library of commonly used units

INSTALLATION
============

Get the code with one of these options.

1. with git:

        git clone git@github.com:jkitchin/matlab-cmu.git

2. as a zip file:

        wget -O matlab-cmu.zip https://github.com/jkitchin/matlab-cmu/zipball/master
        unzip matlab-cmu.zip

Now, make sure to add the directory formed in 1 or 2 to your Matlab path.

Example usage
=============

[See this m-file](https://github.com/jkitchin/matlab-cmu/blob/master/%2Bcmu/examples/unit_tutorials.m)

Some gotchas
============

Not every function in Matlab is overloaded to use units. That means some functions may silently drop the units, and leave you with a plain number. I do not know a solution to that problem. Some functions will not work with units, e.g. det, eig, and others. I haven't implemented those yet.

Mixed algebra with dimensionless numbers and units. There are times where it is essential to tell matlab a number is dimensionless, e.g. 0.5*u.dimensionless. This is necessary so that the appropriate overloaded function is selected.

units can make Matlab very slow! In the simple cases, using units is nearly as fast as not using them. However, if you have vectors or matrices with mixed units, any linear algebra will be slow because it is coded in the units package through loops.

units are just not defined for some operations, e.g. an LU decomposition. you will have to figure out what you mean and how to use units in cases like that.