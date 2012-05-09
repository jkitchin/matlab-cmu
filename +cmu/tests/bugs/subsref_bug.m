%% bug in subsref
%{
when indexing a matrix of mixed row/columns, the units are not
properly indexed out.

Y = [unit1; unit2]

Y(:,1) has both units attached to it.
%}
