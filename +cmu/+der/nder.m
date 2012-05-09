% Numerical derivative of a function handle by finite differentiation

function dfdx = nder(functionhandle,x)
% Compute the derivative of the function handle at x
% so far, only a one-dimensional function is handled
k = [-2 -1 0 1 2];
h = 0.1;
xx = h*k + x;
f = arrayfun(functionhandle,xx);
dfdx = 1/(12*h)*(f(1)-8*f(2) + 8*f(4) - f(5));

% now we iterate by shrinking h and stoping when the difference in the
% new dfdx and old dfdx is sufficiently small.
count = 0;
while count <= 10
    count = count + 1;
    h = h/2;
    xx = h*k + x;
    f = arrayfun(functionhandle,xx);
    dfdx_new = 1/(12*h)*(f(1)-8*f(2) + 8*f(4) - f(5));
    if abs(dfdx_new - dfdx) < 1e-6
        return
    else
        dfdx = dfdx_new;
    end
end
 
error('10 iterations exceeded without converging in derivative')
    
end