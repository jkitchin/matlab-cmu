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
4. Small library of commonly used units

INSTALLATION
============

Get the code with one of these options.

1. with git:

        git clone git@github.com:johnkitchin/matlab-cmu.git

2. as a zip file:

        wget -O matlab-cmu.zip https://github.com/johnkitchin/matlab-cmu/zipball/master
        unzip matlab-cmu.zip

Now, make sure to add the directory formed in 1 or 2 to your Matlab path.

Example usage
=============

[See this html](+cmu/examples/html/unit_tutorials.html)