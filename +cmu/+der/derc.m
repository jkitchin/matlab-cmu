function dydx = derc(x,y)
% dydx = derc(y)
% dydx = derc(x,y) % numerical derivative by centered difference formula.
% compute the numerical derivative of the function y(x) using a centered
% difference formula. The end points are computed with backward
% differences. The function requires that data are in increasing order, and
% non repeating. It is best if the x-data is even spaced, but it is not
% required.

dydx = zeros(size(y)); %preallocate the array

% check that data are sorted in ascending order
if any(sort(x) ~= x)
    error('x data is not sorted in ascending order')
end

% check that there are no repeats
if any(diff(x) == 0)
    error('Repeated x-values detected')
end

% compute the derivative using centered difference
if double(abs((x(2) - x(1)) - (x(3) - x(2)))) < 1e-5 % even spacing
    dydx(1) = (-3*y(1) + 4*y(2) - y(3))/(x(3) - x(1)); 
else
    dydx(1) = (y(2) - y(1))/(x(2) - x(1));
end

dydx(2:end-1) = (y(3:end) - y(1:end-2))./(x(3:end) - x(1:end-2));

if double(abs((x(end) - x(end-1)) - (x(end-2) - x(end-1)))) < 1e-5 % even spacing
    dydx(end) = (y(end-2) - 4*y(end-1) + 3*y(end))/(x(end-2) - x(end)); 
else
    dydx(end) = (y(end) - y(end-1))/(x(end) - x(end-1));
end

