%% class to implement units for Matlab
% see help cmu.unit.units
% or  help cmu.units
%
% Usually used as (both of these do the same thing):
% u = cmu.units;
% u = cmu.unit.units; 

%{
John Kitchin (jkitchin@andrew.cmu.edu)
Started December 1, 2010

Contributing authors: Bruno Calfa

KNOWN BUGS:
    Many functions in Matlab strip the units by treating the unit object as
    a double and returning a double. 
%}

classdef (InferiorClasses = {?double}) unit < double
    
    properties (Constant, Hidden)
        TOLERANCE = 1e-6 %used to determine zero and integers
    end
    
    properties 
        exponents  %exponents on each of the base_units that defines the unit.
        displaystring=''    % a string for display
    end
    
    methods
        function U = unit(data,exponents,displaystring)
            % constructor function - users do not usually use this.
            % data is the numeric data you are adding units too
            % exponents is either:
            % 1. a string of a unit, e.g. 'm', 'kg'. Not 'kg/s'
            % 2. a single cell containing the exponent vector, in which
            % case all data has the same units
            % 3. a multidimensional cell with an exponent vector for each
            % data point.
            %            
            % displaystring is used for display only, particularly when you
            % want to display in non-base units. This string may not be
            % useful for products of units.
            %
            % some examples:
            % U = cmu.unit(1,'m','m') defines a base unit of a meter
            % U = cmu.unit(0.5,[1 -1 0 0 0 0],'m/s') defines 0.5*m/s
            
            U = U@double(data); %initialize superclass with the data
                        
            if nargin > 1
                if isempty(cmu.unit.base_units)
                    base_units = cmu.unit.base_units('SI');
                else
                    base_units = cmu.unit.base_units;
                end
                
                if ischar(exponents)
                    U.exponents = {strcmp(exponents,base_units)};
                elseif iscell(exponents)
                    U.exponents = exponents;
                elseif isvector(exponents)
                    U.exponents = {exponents};
                elseif isempty(exponents)
                    U.exponents = [];
                else
                    error('that kind of exponent is not supported yet')
                end
                
                if nargin == 3
                    U.displaystring = displaystring;
                end
            end         
        end
                
        function s = display(U,format)
            % pretty print the units
            if nargin == 1
                format = '%g';
            end
            
            n = length(U);
            if n == 1
                % print a single unit
                s = sprintf(format,double(U));
                
                e1 = U.exponents;
                
                if ~isempty(e1)
                    e = e1{1};
                else
                    e = [0 0 0 0 0 0 0];
                end
                
                for i=1:numel(e) %loop over each exponent
                    %check for exponent greater than zero
                    if abs(e(i)) > cmu.unit.TOLERANCE
                        %check if it is an integer
                        if abs(mod(e(i),1)) < cmu.unit.TOLERANCE
                            % it is an integer. check if it is one
                            if abs(e(i)-1) < cmu.unit.TOLERANCE
                                s = strcat(s,'*%s'); %+1
                                s = sprintf(s,cmu.unit.base_units{i});
                            elseif abs(e(i)+1) < cmu.unit.TOLERANCE
                                s = strcat(s,'/%s'); %-1
                                s = sprintf(s,cmu.unit.base_units{i});
                            elseif e(i) > 1
                                s = strcat(s,'*%s^%i');
                                s = sprintf(s,cmu.unit.base_units{i},e(i));
                            elseif e(i) < -1
                                s = strcat(s,'/%s^%i');
                                s = sprintf(s,cmu.unit.base_units{i},abs(e(i)));
                            end
                        else
                            % exponent is a float
                            s = strcat(s,'*[%s^%g]');
                            s = sprintf(s,cmu.unit.base_units{i},e(i));
                        end
                    end
                end
                if (nargout == 0)
                    disp(s)
                end
            else %this means a vector or matrix is displayed
                if all(size(U.exponents) == [1 1])
                    % this means only one unit for all data
                    if (nargout == 0)
                        s1 = num2str(double(U));
                        s = '';
                        sz = size(s1);
                        for i = 1:sz(1)
                            s = strcat(s,'|',s1(i,:),'|','\n');
                        end
                        s = strcat(s(1:end-2),'*',U.displaystring,'\n');
                        fprintf(s)
                    end
                elseif all(size(U.exponents) == size(U))
                    % every element has its own unit
                    sz = size(U);
                    s = '|';
                    for i = 1:sz(1)
                        for j = 1:sz(2)
                            s = strcat(s,sprintf('%1.2f\t',U(i,j)),'\t');
                        end
                        s = strcat(s(1:end-2),'|\n');
                    end
                    fprintf(s)  
                    %disp(U)
                else
                    % this means there are different units in a
                    % vector/matrix, but all rows or columns share the
                    % units.
                    disp('matrix with mixed row or column units')
                    double(U)
                    % U.exponents{:}
                end
            end     
        end
        
        function varargout = as(U1, U2, format)
            % function to display U1 as if it had units of U2. returns a
            % string. you will get a funny result for a matrix or vector.
            e1 = U1.exponents;
            e2 = U2.exponents;
            [m,n] = size(e1);
            for i = 1:m
                for j = 1:n
                    if any(abs(e1{i,j} - e2{1})) > cmu.unit.TOLERANCE
                        error('units do not match')
                    end
                end
            end
            
            if nargin < 3
                format = '%1.3f';
            end
            s = sprintf('%s*',format); %gets format string into the string
            s = sprintf(strcat(s,'  ','%s'),double(U1)./double(U2),U2.displaystring);
            if nargout == 0
                disp(s);
            else
                varargout{1} = s;
            end
        end
        
        function U = plus(U1,U2)
            % U1 + U2
            v = double(U1) + double(U2);
            
            if isa(U1,'cmu.unit') && isa(U2,'cmu.unit')
                e1 = U1.exponents;
                e2 = U2.exponents;
                
                % Get dimensions of highest-dimensional argument
                if (numel(e1) > numel(e2))
                    [m,n] = size(e1);
                else
                    [m,n] = size(e2);
                end
                
                % Check if one of the arguments is scalar
                if (numel(e1) == 1)
                    for i=1:m
                        for j=1:n
                            if any(abs(e1{1,1} - e2{i,j})) > cmu.unit.TOLERANCE
                                error('units:plus','incompatible units for addition')
                            end
                        end
                    end
                elseif (numel(e2) == 1)
                    for i=1:m
                        for j=1:n
                            if any(abs(e1{i,j} - e2{1,1})) > cmu.unit.TOLERANCE
                                error('units:plus','incompatible units for addition')
                            end
                        end
                    end
                else
                    for i=1:m
                        for j=1:n
                            if any(abs(e1{i,j} - e2{i,j})) > cmu.unit.TOLERANCE
                                error('units:plus','incompatible units for subtraction')
                            end
                        end
                    end
                end
                
                e = e1;
                ds = U1.displaystring;
            elseif ~isa(U1,'cmu.unit') && isa(U2,'cmu.unit')
                e = U2.exponents;
                ds = U2.displaystring;
            elseif isa(U1,'cmu.unit') && ~isa(U2,'cmu.unit')
                e = U1.exponents;
                ds = U1.displaystring;
            end
            
            U = cmu.unit(v,e,ds);
        end
        
        function U = uplus(U1)
            % +U1
            U = cmu.unit(+double(U1),U1.exponents);
        end

        function U = minus(U1,U2)
            % U1 - U2
            v = double(U1) - double(U2);
            if isa(U1,'cmu.unit') && isa(U2,'cmu.unit')
                e1 = U1.exponents;
                e2 = U2.exponents;  
                
                % Get dimensions of highest-dimensional argument
                if (numel(e1) > numel(e2))
                    [m,n] = size(e1);
                else
                    [m,n] = size(e2);
                end
                
                % Check if one of the arguments is scalar
                if (numel(e1) == 1)
                    for i=1:m
                        for j=1:n
                            if any(abs(e1{1,1} - e2{i,j})) > cmu.unit.TOLERANCE
                                error('units:minus','incompatible units for subtraction')
                            end
                        end
                    end
                elseif (numel(e2) == 1)
                    for i=1:m
                        for j=1:n
                            if any(abs(e1{i,j} - e2{1,1})) > cmu.unit.TOLERANCE
                                error('units:minus','incompatible units for subtraction')
                            end
                        end
                    end
                else
                    for i=1:m
                        for j=1:n
                            if any(abs(e1{i,j} - e2{i,j})) > cmu.unit.TOLERANCE
                                error('units:minus','incompatible units for subtraction')
                            end
                        end
                    end
                end
                
                e = e1;
                ds = U1.displaystring;
            elseif ~isa(U1,'cmu.unit') && isa(U2,'cmu.unit')
                e = U2.exponents;  
                ds = U2.displaystring;
            elseif isa(U1,'cmu.unit') && ~isa(U2,'cmu.unit')
                e = U1.exponents;  
                ds = U1.displaystring;
            end
            U = cmu.unit(v,e,ds);
        end
        
        function U = uminus(U1)
            % -U1
            U = cmu.unit(-double(U1),U1.exponents);
        end
        
        function U = mtimes(U1, U2)
            % U1*U2
            if ~isa(U1,'cmu.unit')
                % constant * U2
                %'case a'
                U = cmu.unit(U1*double(U2),U2.exponents,U2.displaystring);
                return
            elseif ~isa(U2,'cmu.unit')
                % U1*constant
                %'case b'
                U = cmu.unit(double(U1)*U2,U1.exponents,U1.displaystring);
                return
            end
            
            % these cases are for units times units
            % we need to check the size of the units exponents to decide
            % what to do.
            % case 1: single unit * single unit
            %        for this, we can use fast matrix algebra, and handle
            %        the units after
            if isa(U1,'cmu.unit') && isa(U2,'cmu.unit')
                e1 = U1.exponents;
                e2 = U2.exponents;
                if all(size(e1) == [1 1]) && all(size(e2) == [1 1])
                    % 'case 1'
                    v = double(U1)*double(U2);
                    e = e1{1} + e2{1};
                    % we need parentheses to make sure units are
                    % grouped correctly. we check if a displaystring
                    % should be grouped with parentheses. any display
                    % string with a *,/,^ or ( in it should be grouped
                    % with parentheses
                    t = regexp(U2.displaystring,'\^|*|/|(', 'once');
                    if isempty(t)
                        ds = strcat(U1.displaystring,'*',U2.displaystring);
                    else
                        ds = strcat(U1.displaystring,'*(',U2.displaystring,')');
                    end
                    U = cmu.unit(v,e,ds);
                    if cmu.unit.isDimensionless(U)
                       U = double(U);
                    end
                    return
                else
                    % 'case 2'
                    % at least one of the units has a unit defined for each
                    % element. Here, we use the definition of matrix
                    % multiplication to perform the multiplication.
                    % 'case2'
                    [m1 p1] = size(e1);
                    [p2 n2] = size(e2);
                    
                    % The first unit is a scalar
                    if isscalar(double(U1))
                        %'scalar'
                        e = cell(p2,n2);
                        for i = 1:p2
                            for j = 1:n2
                                e{i,j} =  e1{1,1} + e2{1,j};
                            end
                        end
                        
                        v = double(U1)*double(U2);
                        t = regexp(U2.displaystring,'\^|*|/|(', 'once');
                        if isempty(t)
                            ds = strcat(U1.displaystring,'*',U2.displaystring);
                        else
                            ds = strcat(U1.displaystring,'*(',U2.displaystring,')');
                        end
                        U = cmu.unit(v,e,ds);
                        if cmu.unit.isDimensionless(U)
                           U = double(U);
                        end
                    elseif isscalar(double(U2))
                        %'u2 scalar'
                        e = cell(m1,p1);
                        for i = 1:m1
                            for j = 1:p1
                                e{i,j} =  e1{i,1} + e2{1,1};
                            end
                        end
                        
                        v = double(U1)*double(U2);
                        t = regexp(U2.displaystring,'\^|*|/|(', 'once');
                        if isempty(t)
                            ds = strcat(U1.displaystring,'*',U2.displaystring);
                        else
                            ds = strcat(U1.displaystring,'*(',U2.displaystring,')');
                        end
                        U = cmu.unit(v,e,ds);
                        if cmu.unit.isDimensionless(U)
                           U = double(U);
                        end
                    else
                        %'u1*u2 matrix with mixed column/row units'
                        [m1 p1] = size(U1);
                        [p2 n2] = size(U2);
                        
                        if p1 ~= p2
                            error('units: improper matrix dimensions for multiplication')
                        end
                        % m1 x p1 matrix times p2 x n2 matrix gives m1 x n2
                        % matrix.
                        U = []*cmu.unit([],{[0 0 0 0 0 0 0]},'');
                        for i = 1:m1
                            for j = 1:n2
                                Uij = 0;
                                for k=1:p1
                                    s = struct; s.type='()'; s.subs={i,k};
                                    u1_ik = subsref(U1,s);
                                    s = struct; s.type='()'; s.subs={k,j};
                                    u2_kj = subsref(U2,s);
                                    
                                    Uij = Uij + u1_ik*u2_kj;
                                end
                                s = struct; s.type = '()'; s.subs={i,j};
                                if isa(Uij,'cmu.unit')
                                    U = subsasgn(U,s,Uij);
                                else
                                    U = subsasgn(U,s,cmu.unit(Uij,{[0 0 0 0 0 0 0]},''));
                                end
                            end
                        end
                        
                        if cmu.unit.isDimensionless(U)
                            % 'dimensionless'
                            U = double(U);
                        end
                        return
                    end            
                end       
            end
        end
        
        function U = times(U1,U2)
            % U1.*U2
            v = double(U1).*double(U2);
            
            if ~isa(U1,'cmu.unit') && isa(U2,'cmu.unit')
                e = U2.exponents;
                ds = U2.displaystring;
            elseif isa(U1,'cmu.unit') && ~isa(U2,'cmu.unit')
                e = U1.exponents;
                ds = U1.displaystring;
            elseif isa(U1,'cmu.unit') && isa(U2,'cmu.unit')
                e1 = U1.exponents;
                e2 = U2.exponents;
                
                [m1,p1] = size(e1);
                [p2,n2] = size(e2);
                
                if ((m1 == p2) && (p1 == n2))
                    %same size arrays of exponents
                    e = cell(m1,n2);
                    for i = 1:numel(e1)
                        e{i} = e1{i} + e2{i};
                    end
  
                elseif ((m1 == 1) && (p1 == 1))
                    % First argument is scalar with units
                    e = cell(p2,n2);
                    for i = 1:p2
                        for j = 1:n2
                            e{i,j} =  e1{1,1} + e2{1,j};
                        end
                    end
                elseif ((p2 == 1) && (n2 == 1))
                    % Second argument is scalar with units
                    
                    e = cell(m1,p1);
                    for i = 1:m1
                        for j = 1:p1
                            e{i,j} =  e1{i,1} + e2{1,1};
                        end
                    end
                else
                    error('exponents are different sizes')
                end
                
                ds = strcat(U1.displaystring,'*(',U2.displaystring,')');
                
                t = regexp(U2.displaystring,'\^|*|/|(','once');
                if isempty(t)
                    ds = strcat(U1.displaystring,'*',U2.displaystring);
                else
                    ds = strcat(U1.displaystring,'*(',U2.displaystring,')');
                end
            end

            U = cmu.unit(v,e,ds);

            if cmu.unit.isDimensionless(U)
               U = double(U);
            end
        end

        function U = mrdivide(U1, U2)
            % U1 / U2
        
            if length(size(U1)) > 2
                error('dimensions greater than 2 not supported yet')
            end

            v = double(U1)/double(U2);            
            % we subtract all the exponents for several cases
            % case 1
            if ~isa(U1,'cmu.unit') && isa(U2,'cmu.unit')
                % number/unit
                e = U2.exponents{:};
                ds = strcat(U2.displaystring,'^-1');
                U = cmu.unit(v,{-e},ds);
            elseif isa(U1,'cmu.unit') && ~isa(U2,'cmu.unit')
                % unit/number
                U = cmu.unit(v,U1.exponents,U1.displaystring);
            elseif isa(U1,'cmu.unit') && isa(U2,'cmu.unit')
                % unit/unit
                e1 = U1.exponents;
                e2 = U2.exponents;
                
                % Get dimensions of highest-dimensional argument
                if (numel(e1) > numel(e2))
                    [m,n] = size(e1);
                else
                    [m,n] = size(e2);
                end
                
                e = cell(m,n);
                
                % Check if one of the arguments is scalar
                if (numel(e1) == 1)
                    for i=1:m
                        for j=1:n
                            e{i,j} = e1{1,1} - e2{i,j};
                        end
                    end
                elseif (numel(e2) == 1)
                    for i=1:m
                        for j=1:n
                            e{i,j} = e1{i,j} - e2{1,1};
                        end
                    end
                else
                    for i=1:m
                        for j=1:n
                            e{i,j} = e1{i,j} - e2{i,j};
                        end
                    end
                end
                
            t = regexp(U2.displaystring,'\^|*|/|(', 'once');
            if isempty(t)
                ds = strcat(U1.displaystring,'/',U2.displaystring);
            else
                ds = strcat(U1.displaystring,'/(',U2.displaystring,')');
            end
            U = cmu.unit(v,e,ds);

            else
                error('unsupported mrdivide operation')
            end
            if cmu.unit.isDimensionless(U)
               U = double(U);
            end
        end

        function U = rdivide(U1, U2)
            % U1 ./ U2
            if length(size(U1)) > 2
                error('dimensions greater than 2 not supported yet')
            end
            
            v = double(U1)./double(U2);
            % we subtract all the exponents for several cases
            % case 1
            if ~isa(U1,'cmu.unit') && isa(U2,'cmu.unit')
                % number/unit
                e = U2.exponents{:};
                ds = strcat('*',U2.displaystring,'^-1');
                U = cmu.unit(v,{-e},ds);
            elseif isa(U1,'cmu.unit') && ~isa(U2,'cmu.unit')
                % unit/number
                U = cmu.unit(v,U1.exponents,U1.displaystring);
            elseif isa(U1,'cmu.unit') && isa(U2,'cmu.unit')
                e1 = U1.exponents;
                e2 = U2.exponents;
                
                % Get dimensions of highest-dimensional argument
                if (numel(e1) > numel(e2))
                    [m,n] = size(e1);
                else
                    [m,n] = size(e2);
                end
                
                e = cell(m,n);
                
                % Check if one of the arguments is scalar
                if (numel(e1) == 1)
                    for i=1:m
                        for j=1:n
                            e{i,j} = e1{1,1} - e2{i,j};
                        end
                    end
                elseif (numel(e2) == 1)
                    for i=1:m
                        for j=1:n
                            e{i,j} = e1{i,j} - e2{1,1};
                        end
                    end
                else
                    for i=1:numel(e)
                        e{i} = e1{i} - e2{i};
                    end
                end
                
                t = regexp(U2.displaystring,'\^|*|/|(','once');
                if isempty(t)
                    ds = strcat(U1.displaystring,'/',U2.displaystring);
                else
                    ds = strcat(U1.displaystring,'/(',U2.displaystring,')');
                end
                
                U = cmu.unit(v,e,ds);
                
            else
                error('unsupported rdivide operation')
            end
            if cmu.unit.isDimensionless(U)
               U = double(U);
            end
        end
        
        function varargout = mldivide(U1,U2)
            % U1 \ U2
            % Check if U1 is scalar
            if isscalar(double(U1))
                % Perform U2/U1
                out = U2/U1;
            else
                % Call linsolve with no options
                out = linsolve(U1,U2);
            end
            
            if (nargout == 0)
                varargout = cell(1);
            end
            varargout{1} = out;
        end
        
        function varargout = ldivide(U1,U2)
            % U1 .\ U2
            out = U2./U1;
            
            if (nargout == 0)
                varargout = cell(1);
            end
            varargout{1} = out;            
        end
        
        function U = mpower(U1,k)
            % U^k
            v = double(U1)^k;
            e = U1.exponents;
            [m,n] = size(e);
            for i = 1:m
                for j = 1:n
                    e{i,j} = e{i,j}.*k;
                end
            end
            ds = strcat(U1.displaystring,sprintf('^%d',k));
            U = cmu.unit(v,e,ds);
        end
        
        function U = power(U1,k)
            % U.^k
            v = double(U1).^k;
            e = U1.exponents;
            [m,n] = size(e);
            for i = 1:m
                for j = 1:n
                    e{i,j} = e{i,j}.*k;
                end
            end
            ds = strcat(U1.displaystring,sprintf('^%d',k));
            U = cmu.unit(v,e,ds);
        end
        
        function U = abs(U1)
            v = abs(double(U1));
            U = cmu.unit(v,U1.exponents,U1.displaystring);
        end
        
        function U = horzcat(varargin)
            % [U1 U2 ...]
            % horizontal concatenation - defined recursively
            if length(varargin) == 1
                % single element so we are done
                U = varargin{1};
                return
            end
            
            % we concatenate two elements at a time
            u1 = varargin{1};
            u2 = varargin{2};
            
            values = [double(u1) double(u2)];
            if isa(u1,'cmu.unit') && isa(u2,'cmu.unit')
                exps = cat(2,u1.exponents, u2.exponents);
            elseif isempty(u1) && isa(u2,'cmu.unit')
                exps = [u2.exponents];
            elseif isempty(u2) && isa(u1,'cmu.unit')
                exps = [u1.exponents];
            else
                exps = {[0 0 0 0 0 0 0]};
                %error('horzcat is not defined for catenation of units and non-units')
            end
            
            U = cmu.unit(values,exps);
            
            if length(varargin) == 2
                return
            else
                % now we get the next ones
                U = horzcat(U,varargin{3:end});
            end
        end
        
        function U = vertcat(varargin)
            % [U1; U2; ...]
            % vertical concatenation - defined recursively
            if length(varargin) == 1
                % single element so we are done
                U = varargin{1};
                return
            end
            
            u1 = varargin{1};
            u2 = varargin{2};
            
            values = [double(u1); double(u2)];
            if isa(u1,'cmu.unit') && isa(u2,'cmu.unit')
                exps = [u1.exponents; u2.exponents];
            else
                error('vertcat for those exponents is not supported')
            end
            
            U = cmu.unit(values,exps);
            
            if length(varargin) == 2
                return
            else
                U = vertcat(U,varargin{3:end});
            end
        end

        function out = subsref(U,s)
           % U(:)
           % U(:,:)
           % U.as()
            switch s(1).type
                case '.'
                    if (strcmp(s(1).subs,'as'))
                        U2 = s(2).subs{1};
                        out = as(U,U2);
                        return
                    else
                        out = PrivateGet(U,s(1).subs);
                    end
                case '()'
                    data = double(U);
                    exps = U.exponents;
                    
                    % this is the indexed data no matter what
                    % what follows here is how we handle the indexing of
                    % the exponents
                    sdata = subsref(data,s(1));
                    
                    % we need to handle the exponents for two cases
                    % case 1, there is a single exponent for all the
                    % data.
                    % case 2: each data point has its own units or there
                    % are mixed exponents.
                    
                    [em,en] = size(exps);
                    if all([em en] == [1 1])
                        % 'one unit for all data.'
                        sexps = exps;
                        out = cmu.unit(sdata,sexps,U.displaystring);
                    elseif all([em en] == size(U))
                        % 'every number has its own units
                        sexps = subsref(exps,s);
                        out = cmu.unit(sdata,sexps,U.displaystring);  
                    else
                        %'third case'
                        % there are vectors (rows or columns) with
                        % different units. we have to modify the indexing on
                        % the exponents to avoid exceeding those
                        % dimensions. The dimensions of the exponent array
                        % should be 1xn if each column has units, and nx1
                        % if each row has units
                        % find index that equals 1
            
                        sz = size(exps);
                        if sz(1) == 1
                            %'column units'
                            N = 1;
                        else
                            %'row units'
                            N = 2;
                        end
                        
                        % there is a problem if you index by element with a
                        % single number. here we convert to i,j notation.
                        % and then set the appropriate element to 1. 
                        %% I think there is a bug in this code. I am not sure what the 10 is
                        if all(size(s.subs) == [1,1])
                            % 'single number indexing'
                            ind = ind2sub(10,size(data));
                            s.subs = {ind(1), ind(2)};
                        end    
                        
                        s.subs{N} = 1; % set right element to 1
                        % now we set the Nth element to 1, since in that
                        % dimension there is only one dimension
                        sexps = subsref(exps,s);
                        out = cmu.unit(sdata,sexps,U.displaystring);  
                    end
            end
            
            % recursive indexing 
            if (length(s) > 1)
                out = subsref(out,s(2:end));
            end

        end
           
        function V = PrivateGet(U,FullPropName)
            switch FullPropName
                case 'exponents'
                    V = U.exponents;
                case 'displaystring'
                    V = U.displaystring;
                otherwise
                    error([FullPropName ': Unknown Property.']);
            end
        end

        function U = subsasgn(U1,S,RHS)
            switch S(1).type
                case '.'
                    if (length(S) == 1)
                        U = PrivateSet(U1,S(1).subs,RHS);
                    else    
                        tmp = PrivateGet(U1,S(1));
                        tmp = subsasgn(tmp,S(2:end),RHS);
                        U = PrivateSet(U1,S(1).subs,tmp);
                    end
                    
                case '()'   
                    if (length(S) == 1)
                        sind = S.subs;
                        
                        if ~isempty(U1)
                            v = double(U1);
                            v(sind{:}) = double(RHS);
                            
                            if isa(U1,'cmu.unit')
                                e = U1.exponents;
                                e(sind{:}) = RHS.exponents(1);
                            else
                                e(sind{:}) = {[0 0 0 0 0 0 0]};
                            end
                        else
                            v = double(RHS);
                            
                            e = [];
                            if isa(RHS,'cmu.unit')
                                e = RHS.exponents;
                            else
                                e = {[0 0 0 0 0 0 0]};
                            end
                        end
                        
                        U = cmu.unit(v,e);
                    else
                        tmp = builtin('subsref',U1,S(1));
                        tmp = subsasgn(tmp,S(2:end),RHS);
                        U = builtin('subsasgn',U1,S(1),tmp);
                    end
            end
            
        end
        
        function U = PrivateSet(U,FullPropName,Value)
            switch FullPropName
                case 'exponents'
                      U.exponents = Value;
                case 'displaystring'
                      U.displaystring = Value;
                otherwise
                    error([FullPropName ': Unknown Property']);
            end
        end
        
        function Z = trapz(X,Y,DIM)
            if nargin == 1
                Y = X;
                Z = cmu.unit(trapz(double(Y)),Y.exponents, Y.displaystring);
            elseif nargin == 2
                if isscalar(Y)
                    DIM = Y; Y = X;
                    Zd = trapz(double(Y),DIM);
                    Z = cmu.unit(Zd,Y.exponents,Y.displaystring);
                else
                    Zd = trapz(double(X),double(Y));
                    U = X*Y';
                    Z = cmu.unit(Zd,U.exponents,U.displaystring);  
                end
            elseif nargin == 3
                % cumtrapz(x,y,dim)
                Zd = trapz(double(X),double(Y),DIM);
                U = X*Y';
                Z = cmu.unit(Zd,U.exponents,U.displaystring);
            end  
        end
                
        function Z = cumtrapz(X,Y,DIM)
            if nargin == 1
                % cumtrapz(Y)
                Y = X;
                Z = cmu.unit(cumtrapz(double(Y)),Y.exponents,Y.displaystring);
                
            elseif nargin == 2
                if isscalar(Y)
                    % cumtrapz(y,dim)
                    DIM = Y; Y = X; 
                    Zd = cumtrapz(double(Y),DIM);
                    Z = cmu.unit(Zd,Y.exponents,Y.displaystring);
                else
                    % cumtrapz(X,Y)
                    Zd = cumtrapz(double(X),double(Y));
                    U = X*Y';
                    Z = cmu.unit(Zd,U.exponents,U.displaystring);    
                end
            elseif nargin == 3
                % cumtrapz(x,y,dim)
                Zd = cumtrapz(double(X),double(Y),DIM);
                U = X*Y';
                Z = cmu.unit(Zd,U.exponents,U.displaystring);
            end   
        end
        
        function yiu = interp1(varargin)
            % Wrapper for interp1
            
            xu = varargin{1};
            Yu = varargin{2};
            xiu = varargin{3};
            
            x = double(xu);
            Y = double(Yu);
            xi = double(xiu);
            
            yi = interp1(x,Y,xi,varargin{4:end});
            
            % Put units back to output
            yiu = cmu.unit(yi,Yu.exponents,Yu.displaystring);
        end
        
        function U = diff(U1,y,dim)
            %overloaded diff to keep units
            data = double(U1);
            if nargin == 1
                v = diff(data);
            elseif nargin == 2
                v = diff(data,y);
            else
                v = diff(data,y,dim);
            end
            U = cmu.unit(v,U1.exponents,U1.displaystring);     
        end
        
        function U = sum(U1,DIM)
            % overloaded sum to keep units. does not work for mixed units,
            % but does not warn you of the error.
            data = double(U1);
            
            if nargin == 1
                v = builtin('sum',data);
            else
                v = builtin('sum',data,DIM);
            end
            
            U = cmu.unit(v,U1.exponents);
   
        end
        
        function U = floor(U1)
            U = cmu.unit(floor(double(U1)),U1.exponents,U1.displaystring);   
        end
        
        function U = ceil(U1)
            U = cmu.unit(ceil(double(U1)),U1.exponents,U1.displaystring);            
        end
        
        function U = round(U1)
            U = cmu.unit(round(double(U1)),U1.exponents,U1.displaystring);   
        end
        
        function U = fix(U1)
            U = cmu.unit(fix(double(U1)),U1.exponents,U1.displaystring);   
        end
        
        function [varargout] = min(varargin)
            % find min of a vector of units

            if nargin == 1
                U = varargin{1};
                [Y, I] = min(double(U));
            
                if nargout <= 1
                    varargout{1} = cmu.unit(Y,U.exponents,U.displaystring);
                elseif nargout == 2
                    varargout{1} = cmu.unit(Y,U.exponents,U.displaystring);
                    varargout{2} = I;
                end
                
            elseif nargin == 2
                U1 = varargin{1};
                U2 = varargin{2};
                Y = min(double(U1),double(U2));
                
                varargout{1} = cmu.unit(Y,U1.exponents,U1.displaystring);
                
            elseif nargin == 3
                U = varargin{1};
                dim = varargin{3};
                
                [Y, I] = min(double(U),[],dim);
            
                if nargout <= 1
                    varargout{1} = cmu.unit(Y,U.exponents,U.displaystring);
                elseif nargout == 2
                    varargout{1} = cmu.unit(Y,U.exponents,U.displaystring);
                    varargout{2} = I;
                end
            end 
        end
        
        function [varargout] = max(varargin)
            % find max of a vector of units

            if nargin == 1
                U = varargin{1};
                [Y, I] = max(double(U));
            
                if nargout <= 1
                    varargout{1} = cmu.unit(Y,U.exponents,U.displaystring);
                elseif nargout == 2
                    varargout{1} = cmu.unit(Y,U.exponents,U.displaystring);
                    varargout{2} = I;
                end
                
            elseif nargin == 2
                U1 = varargin{1};
                U2 = varargin{2};
                Y = max(double(U1),double(U2));
                
                if isa(U1,'cmu.unit')
                    varargout{1} = cmu.unit(Y,U1.exponents,U1.displaystring);
                else
                    varargout{1} = cmu.unit(Y,U2.exponents,U2.displaystring);
                end
                
            elseif nargin == 3
                U = varargin{1};
                dim = varargin{3};
                
                [Y, I] = max(double(U),[],dim);
            
                if nargout <= 1
                    varargout{1} = cmu.unit(Y,U.exponents,U.displaystring);
                elseif nargout == 2
                    varargout{1} = cmu.unit(Y,U.exponents,U.displaystring);
                    varargout{2} = I;
                end
            end 
        end

        function U = sqrt(U1)
            U = U1.^0.5;
        end
        
        function C = ne(U1,U2)
            % A ~= B
            C = ~eq(U1,U2);
        end
        
        function C = eq(U1,U2)
            % A == B
            % units are equal if their values are equal and units are equal
            C = double(U1) == double(U2);
            
            if isa(U1,'cmu.unit') && isa(U2,'cmu.unit')
                % we make sure the units are the same
                e1 = U1.exponents;
                e2 = U2.exponents;
                
                [m1,n1] = size(e1);
                [m2,n2] = size(e2);
                if m1 == m2 && n1 == n2
                    %same size exponents
                    CE = zeros(size(e1));
                    for i=1:m1
                        for j=1:n1
                            CE(i,j) = all(abs(e1{i,j} - e2{i,j})< cmu.unit.TOLERANCE);
                        end
                    end
                    C = C & CE;
                end
                
                
            elseif isa(U1,'cmu.unit') && ~isa(U2,'cmu.unit')
                % comparison with a number. we assume the units are the
                % same
                return
                
            elseif ~isa(U1,'cmu.unit') && isa(U2,'cmu.unit')
                % comparison with a number. we assume the units are the
                % same
                return
            end
        end
        
        function  [P,S,MU] = polyfit(X,Y,N)
            % overloaded polyfit for units
            % units for now are only on the parameters         
            dx = double(X);
            dy = double(Y);
            if nargout <= 1
                dP = polyfit(dx, dy, N);
            elseif nargout == 2
                [dP, S] = polyfit(dx, dy, N);
            elseif nargout == 3
                [dP, S, MU] = polyfit(dx, dy, N);
            end
            % units on parameters           
            P = [1 1]*cmu.unit(1,{[0 0 0 0 0 0 0]},'');
            for i=1:N+1
                % the ith element is P(i)*X^(N - (i - 1)) and the
                % units on that are the units of Y
                x1 = X(1);
                y1 = Y(1);
                if isa(X,'cmu.unit') && isa(Y,'cmu.unit')
                    % make unit units to compute units on parameter
                    ux = cmu.unit(1,x1.exponents,x1.displaystring);
                    uy = cmu.unit(1,y1.exponents,y1.displaystring);
                    U = uy/(ux^(N-(i-1)));
                elseif isa(X,'cmu.unit')
                    %'only units on x'
                    ux = cmu.unit(1,x1.exponents,x1.displaystring);
                    U = 1/(ux^(N-(i-1)));
                    if ~isa(U,'cmu.unit')
                        %'dimensionless'
                        U = cmu.unit(1,{[0 0 0 0 0 0 0]},'');
                    end
                elseif isa(Y,'cmu.unit')
                    %'only y has units'
                    uy = cmu.unit(1,y1.exponents,y1.displaystring);
                    U = uy;
                end
                uP = dP(i)*U;
                %P = [P uP];
                P(i) = uP;            
            end
        end
        
        function pder = polyder(P)
            if nargin > 1
                error('cmu.unit','polyder only supports one polynomial with units')
            end
            
            pd = polyder(double(P));
            
            % now construct the units. the output of polyder has one less
            % dimension than P, and the units are shifted
            
            for i = 1:length(pd)
                unit = P(i+1)/double(P(i+1));
                pder(i) = pd(i)*unit;
            end       
        end
        
        function pint = polyint(P,K)
            % I am not sure we can do the integration because we have to
            % figure out a new units. It may be possible to find a trend in
            % the exponents of each unit, but I havent figured it out yet.
           error('cmu.unit','polyint not supported for units') 
        end
        
        function pval = polyval(P,X)
            if nargout > 1
                error('cmu.unit','delta output is not supported yet')
            end
            
            N = length(P)-1;
            pval = 0;
            for i=1:N
                pval = pval + P(i)*X.^(N-(i-1));
            end      
            pval = pval + P(N+1);
        end
        
        function out = reshape(U,varargin)
            % Overloads reshape for unit objects
            data = double(U);

            data = reshape(data,varargin{:});
            out = cmu.unit(data,U.exponents,U.displaystring);
        end
        
        function str = sprintf(format,varargin)
            % Overloads sprintf for unit objects
            % this is a recursive function that calls itself until all
            % varargin that are units are no longer units.
            
            match = regexp(format,'%\d*\.?\d*\w{1,2}','match');
            
            % the number of matches should equal length of varargin
            if numel(match) ~= numel(varargin)
                error('cmu.unit','inputs do not match number of formats')
            end
            
            strings = {};
            for i = 1:length(varargin)
                u = varargin{i};
                if isa(varargin{i},'cmu.unit')
                    % convert unit to a string
                    strings{i} = u.display(match{i});
                else
                    strings{i} = sprintf(match{i},u);
                end
            end
            
            % now we have to fix up the format string to replace all
            % formats with %s.
            for i = 1:length(match)
                format = regexprep(format,match{i},'%s','once');
            end
            str = sprintf(format,strings{:});
        end
                    
        function U = transpose(U1)
            data = double(U1);
            data = builtin('transpose',data);
            U = cmu.unit(data,U1.exponents',U1.displaystring);
        end
        
        function U = ctranspose(U1)
            data = double(U1);
            data = builtin('ctranspose',data);
            U = cmu.unit(data,U1.exponents',U1.displaystring);
        end
        
        function varargout = ode45(fh,tspan,init,options)
            %wrapper for ode45
            if nargin == 3
                options = [];
            end
            % make sure everything is doubles
            function odeout = odefun(fh,t,y)
                % make sure there are units on t,y
                t1 = tspan(1);
                t = cmu.unit(t,t1.exponents);
                y = cmu.unit(y,init.exponents);
                odeout = double(feval(fh,t,y));
            end
            
            func = @(t,y) odefun(fh,t,y);
            ts = double(tspan);
            y0 = double(init);
            
            if ~iscolumn(y0)
                error('initial conditions must be a column vector!')
            end
            
            [varargout{1:nargout}] = ode45(func, ts, y0, options);
            if nargin == 3
                options = [];
            end
            
            %now we need to put units back onto the output
            if nargout == 1
                % this is a solution structure output
                sol = varargout{1};
                t1 = tspan(1);
                x = cmu.unit(sol.x,t1.exponents);
                y = sol.y;
                % this code allows there to be different units.
                for i=1:length(init)
                    u = cmu.unit(y(i,:),init.exponents{i});  
                    Y(i,:) = u;
                end
                                
                sol.x = x;
                sol.y = Y;
                varargout{1} = sol;
            elseif (nargout >= 2)
                t = varargout{1};
                y = varargout{2};
                t1 = tspan(1);
                varargout{1} = cmu.unit(t,t1.exponents);
                
                % this code allows there to be different units.
                for i=1:length(init)
                    u = cmu.unit(y(:,i),init.exponents{i});
                    Y(:,i) = u;
                end
                varargout{2} = Y;
                
            elseif (nargout == 5)
                TE = varargout{3};
                YE = varargout{4};
                t1 = tspan(1);
                varargout{3} = cmu.unit(TE,t1.exponents);
                warning('cmu.unit','check units on YE. they may not be correct')
                varargout{4} = cmu.unit(YE,init.exponents);
            end
        end
        
        function varargout = ode15s(fh,tspan,init,options)
            %wrapper for ode15s
            %fh = varargin{1};
            %tspan = varargin{2};
            %init  = varargin{3};
            if nargin == 3
                options = [];
            end
            
            % make sure everything is doubles for fh
            function odeout = odefun(fh,t,y)
                % make sure there are units on t,y
                t1 = tspan(1);
                t = cmu.unit(t,t1.exponents);
                y = cmu.unit(y,init.exponents);
                odeout = double(feval(fh,t,y));
            end
            
            func = @(t,y) odefun(fh,t,y);
            ts = double(tspan);
            y0 = double(init);
            
            if ~iscolumn(y0)
                error('initial conditions must be a column vector!')
            end
            [varargout{1:nargout}] = ode15s(func, ts, y0, options);
            
            %now we need to put units back onto the output
            if nargout == 1
                sol = varargout{1};
                               
                x = cmu.unit(sol.x,tspan(1).exponents);
                y = sol.y;
                % this code allows there to be different units.
                for i=1:length(init)
                    u = cmu.unit(y(i,:),init.exponents{i});  
                    Y(i,:) = u;
                end
                                
                sol.x = x;
                sol.y = Y;
                varargout{1} = sol;
                
            elseif (nargout >= 2)
                t = varargout{1};
                y = varargout{2};
                t1 = tspan(1);
                varargout{1} = cmu.unit(t,t1.exponents);
                % this code allows there to be different units.
                for i=1:length(init)
                    u = cmu.unit(y(:,i),init.exponents{i});
                    Y(:,i) = u;
                end
                varargout{2} = Y;
            
            elseif (nargout == 5)
                TE = varargout{3};
                YE = varargout{4};
                t1 = tspan(1);
                varargout{3} = cmu.unit(TE,t1.exponents);
                varargout{4} = cmu.unit(YE,init.exponents);
            end
        end
        
        function [uxint,upxint] = deval(sol,xint,idx)
            % overloaded deval
            y = sol.y;
            
            if nargin == 2
                sold = sol;
                sold.x = double(sol.x);
                sold.y = double(sol.y);
                if nargout == 1
                    uxint = deval(sold,double(xint));
                    uxint = cmu.unit(uxint,y.exponents,y.displaystring);
                else
                    [uxint,upxint] = deval(sold,double(xint));
                    uxint = cmu.unit(uxint,y.exponents,y.displaystring);
                    dydx = (sol.y)./(sol.x);
                    upxint = cmu.unit(upxint,dydx.exponents,dydx.displaystring);
                end
            elseif nargin == 3
                sold = sol;
                sold.x = double(sol.x);
                sold.y = double(sol.y);
                if nargout == 1
                    uxint = deval(sold,double(xint),idx);
                    uxint = cmu.unit(uxint,y.exponents,y.displaystring);
                else
                    [uxint,upxint] = deval(sold,double(xint),idx);
                    uxint = cmu.unit(uxint,y.exponents,y.displaystring);
                    dydx = (sol.y)./(sol.x);
                    upxint = cmu.unit(upxint,dydx.exponents,dydx.displaystring);
                end
            end
        end
        
        function varargout = ode23(varargin)
            error('cmu:unit','ode23 not implemented for units. Use ode45 or ode15s.')
        end
        
        function varargout = ode23s(varargin)
            error('cmu:unit','ode23s not implemented for units. Use ode45 or ode15s.')
        end
        
        function varargout = ode23t(varargin)
            error('cmu:unit','ode23t not implemented for units. Use ode45 or ode15s.')
        end
        
        function varargout = ode113(varargin)
            error('cmu:unit','ode113 not implemented for units. Use ode45 or ode15s.')
        end
        
        function varargout = ode23tb(varargin)
            error('cmu:unit','ode23tb not implemented for units. Use ode45 or ode15s.')
        end
        
        function varargout = ode15i(varargin)
            error('cmu:unit','ode15i not implemented for units. Use ode45 or ode15s.')
        end
        
        function varargout = dde23(varargin)
            error('cmu:unit','dde23 not implemented for units. Consider using dimensionless equations.')
        end
        
        function varargout = pdepe(varargin)
            error('cmu:unit','pdepe not implemented for units. Consider using dimensionless equations.')
        end
        
        function varargout = bvp4c(varargin)
            error('cmu:unit','bvp4c not implemented for units. Consider using dimensionless equations.')
        end
        
        function varargout = bvp5c(varargin)
            error('cmu:unit','bvp5c not implemented for units. Consider using dimensionless equations.')
        end
        
        function varargout = meshgrid(varargin)
            % Wrapper for meshgrid
            if (nargout < 2)
                error('meshgrid: wrong number of output arguments');
            end
            
            meshdata = cell(size(varargin));
            meshunit = meshdata;
            if (nargin >= 1)
                meshdata{1} = double(varargin{1});
                if isa(varargin{1},'cmu.unit')
                    meshunit{1} = varargin{1}.exponents;
                end
                if (nargout ~= 2)
                    error('meshgrid: wrong number of output arguments');
                end
                meshout = cell(1,2);
            end
            if (nargin >= 2)
                meshdata{2} = double(varargin{2});
                if isa(varargin{2},'cmu.unit')
                    meshunit{2} = varargin{2}.exponents;
                end
                if (nargout ~= 2)
                    error('meshgrid: wrong number of output arguments');
                end
                meshout = cell(1,2);
            end
            if (nargin == 3)
                meshdata{3} = double(varargin{3});
                if isa(varargin{2},'cmu.unit')
                    meshunit{3} = varargin{3}.exponents;
                end
                if (nargout ~= 3)
                    error('meshgrid: wrong number of output arguments');
                end
                meshout = cell(1,3);
            end
            if (nargin > 3)
                error('meshgrid takes up to 3 input arguments');
            end
            
            [meshout{:}] = meshgrid(meshdata{:});
            
            varargout = cell(size(meshout));
            for i = 1:nargin
                varargout{i} = cmu.unit(meshout{i},meshunit{i});
            end
            if (nargout > nargin)
                varargout{end} = cmu.unit(meshout{end},meshunit{end});
            end
        end
        
        function varargout = surf(X,Y,Z,varargin)
            % Wrapper for surf
            Xdata = double(X);
            Ydata = double(Y);
            Zdata = double(Z);
            
            h = surf(Xdata,Ydata,Zdata,varargin{:});
            
            if (nargout == 1)
                varargout = h;
            end
            if (nargout > 1)
                error('surf can take up to 1 output argument');
            end
        end
        
        function varargout = surfc(X,Y,Z,varargin)
            % Wrapper for surfc
            Xdata = double(X);
            Ydata = double(Y);
            Zdata = double(Z);
            
            h = surfc(Xdata,Ydata,Zdata,varargin{:});
            
            if (nargout == 1)
                varargout = h;
            end
            if (nargout > 1)
                error('surfc can take up to 1 output argument');
            end
        end
        
        function varargout = surfz(X,Y,Z,varargin)
            % Wrapper for surfz
            Xdata = double(X);
            Ydata = double(Y);
            Zdata = double(Z);
            
            h = surfz(Xdata,Ydata,Zdata,varargin{:});
            
            if (nargout == 1)
                varargout = h;
            end
            if (nargout > 1)
                error('surfz can take up to 1 output argument');
            end
        end
        
        function varargout = mesh(X,Y,Z,varargin)
            % Wrapper for mesh
            Xdata = double(X);
            Ydata = double(Y);
            Zdata = double(Z);
            
            h = mesh(Xdata,Ydata,Zdata,varargin{:});
            
            if (nargout == 1)
                varargout = h;
            end
            if (nargout > 1)
                error('mesh can take up to 1 output argument');
            end
        end
        
        function varargout = meshc(X,Y,Z,varargin)
            % Wrapper for meshc
            Xdata = double(X);
            Ydata = double(Y);
            Zdata = double(Z);
            
            h = meshc(Xdata,Ydata,Zdata,varargin{:});
            
            if (nargout == 1)
                varargout = h;
            end
            if (nargout > 1)
                error('meshc can take up to 1 output argument');
            end
        end
        
        function varargout = meshz(X,Y,Z,varargin)
            % Wrapper for meshz
            Xdata = double(X);
            Ydata = double(Y);
            Zdata = double(Z);
            
            h = meshz(Xdata,Ydata,Zdata,varargin{:});
            
            if (nargout == 1)
                varargout = h;
            end
            if (nargout > 1)
                error('meshz can take up to 1 output argument');
            end
        end
        
        function [q,varargout] = quad(fun,a,b,varargin)
            % Wrapper for quad
            % Convert inputs to double
            adata = double(a);
            bdata = double(b);
            
            function quadfunout = quadfun(fh,x)
                % Make sure there are units on x
                if isa(a,'cmu.unit')
                    x = cmu.unit(x,a.exponents,a.displaystring);
                end
                quadfunout = double(feval(fh,x));
            end
            
            % Get function units
            funit = feval(fun,a)*a;
            
            if ~isa(funit,'cmu.unit')
                % this can happen if the integral is dimensionless
                funit = cmu.unit(1,[0 0 0 0 0 0 0],'');
            end
            % Make function handle return double
            fh = @(x) quadfun(fun,x);
            
            % Call builtin quad
            if nargout == 0
                quadout = cell(1,1);
                quadout{1} = quad(fh,adata,bdata,varargin{:});
            else
                quadout = cell(1,nargout);
                [quadout{:}] = quad(fh,adata,bdata,varargin{:});
            end
            
            % Put units back on output
            if (nargin >= 1)
                q = cmu.unit(quadout{1},funit.exponents,funit.displaystring);
            end
            
            if(nargin == 2)
                varargout = quadout{2};
            end
        end
        
        function [x,varargout] = fsolve(varargin)
            % Wrapper for fsolve
            if nargout >= 1
                fsolveout = cell(1,nargout);
                varargout = cell(1,nargout - 1);
            end
            
            % Parse input arguments
            if (nargin == 1)
                fun = varargin.objective;
                x0 = varargin.x0;
                options = varargin.options;
            else
                if (nargin >= 2)
                    fun = varargin{1};
                    x0 = varargin{2};
                end
                
                options = [];
                if (nargin == 3)
                    options = varargin{3};
                end
            end
            
            % Objective function units
            fvalunit = feval(fun,x0);

            % Convert to double
            x0data = double(x0);
            
            function nlaout = nlafun(fh,x)
                % Make sure there are units on x
                x = cmu.unit(x,x0.exponents,x0.displaystring);
                nlaout = double(feval(fh,x));
            end
            
            fsolvefun = @(x) nlafun(fun,x);

            % Call builtin fsolve
            %fsolve(fsolvefun,x0data,options)
            if nargout == 0
                x = fsolve(fsolvefun,x0data,options);
                x = cmu.unit(x,x0.exponents,x0.displaystring);
            else
                [fsolveout{:}] = fsolve(fsolvefun,x0data,options);
                x = cmu.unit(fsolveout{1},x0.exponents,x0.displaystring);
            end
            
            % Process outputs
                    
            if (nargout >= 2)
                varargout{1} = cmu.unit(fsolveout{2},fvalunit.exponents,fvalunit.displaystring);
            end
            if (nargout >= 3)
                [varargout{2:end}] = fsolveout{3:end};
            end
        end
        
        function [uX, uR] = linsolve(A,b,opts)
            % wrapper for linsolve to preserve units

            if nargin == 2
                [X, R] = linsolve(double(A),double(b));
            elseif nargin == 3
                [X, R] = linsolve(double(A),double(b),opts);
            end

            %units on x, in Ax=b. the idea is to compute b(i)/A(1,i) for each i in the
            %output. that should define the units for each element of X.
            
            % check if well-posed problem
            [m1 n1] = size(A);
            [m2 n2] = size(b);
            if (m1 == n1) && (n1 == m2)
                % regular LA problem
                e = cell(size(b));
                for i=1:size(A,2) % operate on columns on A = rows in x
                    s = struct; s.type='()'; s.subs={i,1};
                    a = subsref(A,s);
                    s = struct; s.type='()'; s.subs={1};
                    bi = subsref(b,s);
                    % sprintf('a = %f',a)
                    % sprintf('bi = %f',bi)
                    u = bi/a;
                    
                    if isa(u,'cmu.unit')
                        e(i) = u.exponents;
                        ds{i} = u.displaystring;
                    else
                        %'setting dimensionless in linsolve'
                        e(i) = {[0 0 0 0 0 0 0]};
                        ds{i}='';
                    end
                end
            else
                % 'least squares problem'
                % rerun as square LA problem
                [uX uR] = linsolve(A'*A,A'*b);   
                return
            end
            
            if nargout == 0
                % no output, just print unit at command line
                cmu.unit(X,e,ds)
            end
            
            if nargout >= 1
                uX = cmu.unit(X,e,ds);
            end
            
            if nargout == 2
                uR = R;
            end
        end
        
        function [x,varargout] = fzero(varargin)
            % Wrapper for fzero
            if nargout > 0
                fzeroout = cell(1,nargout);
                varargout = cell(1,nargout - 1);
            end
            
            % Parse input arguments
            if (nargin == 1)
                fun = varargin.objective;
                x0 = varargin.x0;
                options = varargin.options;
            else
                if (nargin >= 2)
                    fun = varargin{1};
                    x0 = varargin{2};
                end
                
                options = [];
                if (nargin == 3)
                    options = varargin{3};
                end
            end
            
            % Objective function units
            if numel(x0) == 1
                fvalunit = feval(fun,x0);
            else
                 fvalunit = feval(fun,x0(1));
            end

            % Convert to double
            x0data = double(x0);
            
            function nlaout = nlafun(fh,x)
                % Make sure there are units on x
                x = cmu.unit(x,x0.exponents,x0.displaystring);
                nlaout = double(feval(fh,x));
            end
            
            fzerofun = @(x) nlafun(fun,x);

            % Call builtin fsolve
            if nargout == 0
                x = fzero(fzerofun,x0data,options);
                x = cmu.unit(x,x0.exponents,x0.displaystring);
            else
                [fzeroout{:}] = fzero(fzerofun,x0data,options);
                x = cmu.unit(fzeroout{1},x0.exponents,x0.displaystring);
            end
            if (nargout >= 2)
                varargout{1} = cmu.unit(fzeroout{2},fvalunit.exponents,fvalunit.displaystring);
            end
            if (nargout >= 3)
                [varargout{2:end}] = fzeroout{3:end};
            end
        end
        
        function y = linspace(a,b,varargin)
            % Wrapper for linspace
            % Get data
            adata = double(a);
            bdata = double(b);
            
            % Call builtin linspace
            y = linspace(adata,bdata,varargin{:});
            
            % Put units back on output
            y = cmu.unit(y,a.exponents,a.displaystring);
        end
        
        function y = logspace(a,varargin)
            % Wrapper for logspace
            % Get data
            adata = double(a);
            extraargs = varargin;
            if (varargin{1} ~= pi)
                extraargs{1} = double(b);
            end
            
            % Call builtin logspace
            y = logspace(adata,extraargs{:});
            
            % Put units back on output
            y = cmu.unit(y,a.exponents,a.displaystring);
        end
        
        function [B,BINT,R,RINT,STATS] = regress(Y, X, alpha)
            B = X\Y; % this should automatically have units.
            
            if nargin == 2
                alpha = 0.05;
            end
            
            [b bint r rint stats] = regress(double(Y), double(X), alpha);
            if nargout >= 2
                % assign units of B to rows of BINT
                for i = 1:size(B,1)
                    BINT(i,:) = bint(i,:)*B(i)/double(B(i));
                end
            end
            
            if nargout >= 3
                warning('units','units were not applied to R')
                R = r;
            end
            
            if nargout >= 4
                warning('units','units were not applied to RINT')
                RINT = rint;
            end
            
            if nargout == 5
                warning('units','units were not applied to STATS')
                STATS == stats
            end
            
        end
        
        function [BETA R J COVB MSE] = nlinfit(X,y,model,beta0,options)
            warning('units',{'if the units on the parameters are '
                'algebraically connected, there is a good chance they will'
                'be wrong. the units on the parameter output will be the'
                'same as the units of the initial guess'})
            
            if nargin == 4
                options = [];
            end
            
            thismodel = @(pars, x) double(model(pars,x));
            [beta,r,jacobian,covb,mse] = nlinfit(double(X), double(y), thismodel, double(beta0), options);
            
            % put units onto beta from beta0
            for i=1:length(beta)
                BETA(i) = beta(i)*beta0(i)/double(beta0(i));
            end
            
            if nargout >= 2
                % these are residuals
                R = model(BETA,X) - y;
            end
            
            if nargout >= 3
                warning('units: no units applied to J')
                J = jacobian;
            end
            
            if nargout >= 4
                warning('units: no units applied to covb')
                COVB = covb;
            end
            
            if nargout == 5
                warning('units: no units applied to covb')
                MSE = mse;
            end
            
            
        end
        
        function CI = nlparci(beta,resid,varargin)
            ci = nlparci(double(beta),double(resid),varargin{:});
            for i=1:length(ci)
                CI(i,:) = ci(i,:)*beta(i)/double(beta(i));
            end
            
        end
        
        function [YPRED DELTA] = nlpredci(model,X,beta,resid,varargin)
             thismodel = @(pars, x) double(model(pars,x));
            [ypred delta] = nlpredci(thismodel,double(X),double(beta),double(resid),varargin{:});
            
            Y = model(beta,X);
            for i=1:length(Y)
                YPRED(i,:) = ypred(i,:)*Y(i,:)/double(Y(i,:));
                DELTA(i,:) = delta(i,:)*Y(i,:)/double(Y(i,:));
            end         
        end
        
        function det(~)
            error('det with units is not defined')
        end
        
        function eig(~)
            error('eig with units is not defined')
        end
        
        function eigs(U)
            error('eigs with units is not defined')
        end
        
        function svd(U)
            error('svd with units is not defined')
        end
        
        function gsvd(U)
            error('gsvd with units is not defined')
        end
        
        function svds(U)
            error('svds with units is not defined')
        end
     
        function inv(~)
            error('inverse with units is not defined')
        end
        
        %% Trig, log and exponential functions
        function sin(~)
            error('Dimensionless Argument Required')
        end
        
        function cos(~)
            error('Dimensionless Argument Required')
        end
        
        function tan(~)
            error('Dimensionless Argument Required')
        end
        
        function exp(~)
            error('Dimensionless Argument Required')
        end
        
        function log(~)
            error('Dimensionless Argument Required')
        end
        
        function log10(~)
            error('Dimensionless Argument Required')
        end
        
        function acos(~)
            error('Dimensionless Argument Required')
        end
        
        function acsc(~)
            error('Dimensionless Argument Required')
        end
        
        function asin(~)
            error('Dimensionless Argument Required')
        end
        
        function asec(~)
            error('Dimensionless Argument Required')
        end
        
        function atan(~)
            error('Dimensionless Argument Required')
        end
        
        function acot(~)
            error('Dimensionless Argument Required')
        end
        
        function sec(~)
            error('Dimensionless Argument Required')
        end
        
        function csc(~)
            error('Dimensionless Argument Required')
        end
        
        function cot(~)
            error('Dimensionless Argument Required')
        end
        
        function sinh(~)
            error('Dimensionless Argument Required')
        end
        
        function asinh(~)
            error('Dimensionless Argument Required')
        end
        
        function cosh(~)
            error('Dimensionless Argument Required')
        end
        
        function acosh(~)
            error('Dimensionless Argument Required')
        end
        
        function tanh(~)
            error('Dimensionless Argument Required')
        end
        
        function atanh(~)
            error('Dimensionless Argument Required')
        end
        
        function sech(~)
            error('Dimensionless Argument Required')
        end
        
        function asech(~)
            error('Dimensionless Argument Required')
        end
        
        function csch(~)
            error('Dimensionless Argument Required')
        end
        
        function acsch(~)
            error('Dimensionless Argument Required')
        end
        
        function coth(~)
            error('Dimensionless Argument Required')
        end
        
        function acoth(~)
            error('Dimensionless Argument Required')
        end
                
        function boolean = gt(U1,U2)
            % U1 > U2 or gt(U1,U2)
            % for this comparison to work, the units must be the same
            e1 = U1.exponents;
            e2 = U2.exponents;
            if any(abs(e1{:} - e2{:})) > cmu.unit.TOLERANCE
                error('units do not match')
            end
            boolean = double(U1) > double(U2);
        end
        
        function boolean = ge(U1,U2)
            % U1 >= U2 or gt(U1,U2)
            % for this comparison to work, the units must be the same
            e1 = U1.exponents;
            e2 = U2.exponents;
            if any(abs(e1{:} - e2{:})) > cmu.unit.TOLERANCE
                error('units do not match')
            end
            boolean = double(U1) >= double(U2);
        end
        
        function boolean = lt(U1,U2)
            % U1 < U2 or lt(U1,U2)
            % for this comparison to work, the units must be the same
            e1 = U1.exponents;
            e2 = U2.exponents;
            if any(abs(e1{:} - e2{:})) > cmu.unit.TOLERANCE
                error('units do not match')
            end
            boolean = double(U1) < double(U2);
        end

        function boolean = le(U1,U2)
            % U1 <= U2 or le(U1,U2)
            % for this comparison to work, the units must be the same
            e1 = U1.exponents;
            e2 = U2.exponents;
            if any(abs(e1{:} - e2{:})) > cmu.unit.TOLERANCE
                error('units do not match')
            end
            boolean = double(U1) <= double(U2);
        end

    end
    
    methods (Static)
        
        function NotImplemented(msg)
            % simple error message
            error('unit:NotImplemented','Not implemented: %s',msg)
        end
        
        function bool = isDimensionless(U)
            if ~isa(U,'cmu.unit')
                bool = true;
                return
            end
            
            % returns true if all the exponents are within tolerance of zero
            exps = U.exponents;
            
            % walk through the cell. only up to 2D is supported until I can
            % figure out how to do this recursively or by vectorization
            [m,n] = size(exps);
            for i = 1:m
                for j = 1:n
                    e = exps{i,j};
                    if any(abs(e)) > cmu.unit.TOLERANCE
                        bool = false;
                        return
                    end
                end
            end
            bool = true;
        end
        
        function out = base_units(varargin)
            % set the base units - these are defined as 1 for conversion
            % factors.
            % use a string: 'MKS', 'SI', 'CGS' or 'American'
            % or a cell of units {'Bohr', 's', 'kg'}
            % if no argument is used, the current base_units are returned.
            persistent base_units;
            if nargin == 1
                a = varargin{1};
                if strcmp(a,'MKS') || strcmp(a,'SI')
                    base_units = {'m','s','kg', 'K', 'mol','coul','candela'};
                elseif strcmp(a,'CGS')
                    base_units = {'cm','s','gm','K', 'mol','coul','candela'};
                elseif strcmp(a,'American')
                    base_units = {'in','s','lb','R', 'mol','coul','candela'};
                else
                    error('baseunits must be a string ''MKS'', ''SI'', ''CGS'',''American'' or a cell of strings {''m'',''s'',''kg''}')
                end
            elseif nargin > 1
                base_units = varargin;
            end
            out = base_units;
        end
        
        function u = units(varargin)
            %Units package for Matlab
            %
            % unit algebra is enforced:
            %  1. dimensional consistency is required for adding and
            %  subtracting. 
            %  2. unit powers are kept track of for multiplication and
            %  division.
            %  3. Many matlab functions are overloaded to keep units.
            %
            % base units for length, time, mass, temperature, moles and
            % charge are defined. Temperature is a special unit, only
            % Kelvin and Rankine are supported as base units. Helper
            % functions to convert Fahrenheit and Celcius to these absolute
            % scales are provided. Derived units for many quantities
            % including energy, force, pressure, volume, etc... are also
            % defined. 
            %
            % typical usage:
            % u = cmu.unit.units;      %returns SI units
            % u = cmu.unit.units('CGS')  % returns a units with cm, gm, s as
            % base units.
            % u = cmu.unit.units('mm','min','ton','R','mol','coul','candela') % sets the
            % arguments as the base units. Note that you must define a base
            % unit for every unit type. You cannot specify lbmol as a
            % base_unit because it is a derived unit.
            %
            % units are normally displayed in the base units. 
            %
            % >> a = 2*u.J
            % 2*m^2/s^2*kg
            % 
            % to change the display use the as function
            % a.as(u.J)
            % 2.000 J
            %
            % Printing formatted units in strings
            % >> a=5*u.kg;
            % >> sprintf('The total mass is %1.2fs', a)
            %
            % plotting units works as expected, but the magnitude of the
            % number is determined by the base units. For example, in SI
            % units the molar flow would have units of mol/m^3 and volume
            % would be in m^3. to plot the molar flow in mol/min vs. the
            % volume in L, use this syntax:
            % >> plot(V/u.L,Fa/(u.mol/u.min))
            %
            % Defining your own units
            % If the unit can be derived from existing units, then define
            % your new unit as products of existing units. For example, 
            % 760mmHg = 1atm, so we define a new unit as:
            % >> u.mmHg = 1/760*u.atm;
            % 
            % it is not currently possible for users to define new base
            % units, e.g. dollars. Defining new base units requires
            % modification of this code including adding a new
            % containers.Map for the new base unit, and adding the new type
            % to all_units.            
            
            if nargin >= 1
                base_units = cmu.unit.base_units(varargin{:});
            else
                if isempty(cmu.unit.base_units)
                    %if base_units havent been set, we use SI by default
                    base_units = cmu.unit.base_units('SI');
                else
                    base_units = cmu.unit.base_units;
                end
            end
            
            %% Length units
            LENGTH = containers.Map();
            
            % conversion factors for length
            LENGTH('km') = 1000;
            LENGTH('m') = 1;
            LENGTH('dm')=1e-1;
            LENGTH('cm') = 1e-2;
            LENGTH('mm') = 1e-3;
            LENGTH('um') = 1e-6;
            LENGTH('nm') = 1e-9;
            LENGTH('angstrom') = 1e-10;
            LENGTH('a0') = 0.529e-10*LENGTH('m');
            LENGTH('Bohr') = LENGTH('a0');
            LENGTH('in') = 2.54*LENGTH('cm');
            LENGTH('mil') = 1e-3*LENGTH('in');
            LENGTH('ft') = 12*LENGTH('in');
            LENGTH('yd') = 3*LENGTH('ft');
            LENGTH('mile') = 5280*LENGTH('ft');
            LENGTH('furlong') = 660*LENGTH('ft');
            LENGTH('chain') = 66*LENGTH('ft');
            
            MASS = containers.Map();
            MASS('kg') = 1e3;
            MASS('gm') = 1;
            MASS('mg') = 1e-3;
            MASS('lb') = 0.45359237*MASS('kg');
            MASS('lbm') = MASS('lb');
            MASS('oz') = (1/16)*MASS('lb');
            MASS('amu') = 1.660538782e-27*MASS('kg');
            MASS('ton') = 2000*MASS('lb');
            MASS('tonne') = 1000*MASS('kg');
            MASS('longton') = 2240*MASS('lb');
            
            TIME = containers.Map();
            TIME('s') = 1;
            TIME('min') = 60;
            TIME('hr') = 60*TIME('min');
            TIME('day') = 24*TIME('hr');
            TIME('week') = 7*TIME('day');
            TIME('year') = 365.242199*TIME('day');
            
            TEMPERATURE = containers.Map();
            TEMPERATURE('K') = 1;
            TEMPERATURE('R') = 5/9*TEMPERATURE('K'); % I do not understand 
                                                     % why this is 5/9. I
                                                     % think it should be
                                                     % 9/5, but then no
                                                     % conversions work.
            TEMPERATURE('dC') = TEMPERATURE('K'); % relative degree C
            TEMPERATURE('dF') = TEMPERATURE('R'); % relative degree F
            
            MOL = containers.Map();
            % you cannot define lbmol or kgmol here because the masses may
            % not be normalized to the chosen base unit above. these units
            % are defined after the base-units are defined and all other
            % units are normalized. 
            MOL('mol') = 1;
            MOL('kmol') = 1000;
            MOL('mmol') = 1e-3;
            
            CHARGE = containers.Map();
            CHARGE('coul') = 1;
            
            LUMINOSITY = containers.Map();
            LUMINOSITY('cd') = 1;
                        
            all_units = {LENGTH, TIME, MASS, TEMPERATURE, MOL, CHARGE, LUMINOSITY};
            
            % now we need to normalize each category by the base_units
            if length(base_units) ~= length(all_units)
                error('user-defined base_units do not equal defined units')
            end
            
            % this computes the index for each base unit to figure out how
            % to normalize each one
            for i=1:length(base_units)
                for j=1:length(all_units)
                    map = all_units{j};
                    if isKey(map,base_units{i})
                        k = keys(map);
                        v = cell2mat(values(map))/map(base_units{i});
                        for jj = 1:length(map)
                            % u.kg = 1
                            % u(1).(k{i}) = v(i); %simple structure
                            data = v(jj);
                            disp_string = k{jj};
                            u(1).(disp_string) = cmu.unit(data,base_units{i},disp_string); 
                        end
                    end
                end
            end
                     
            % Only derived units should go after here. These all depend on
            % the base units having been normalized 
            
            % some derived mole units
            u.lbmol = u.lb/u.gm*u.mol;
            u.gmmol = u.gm*u.mol;
            u.kgmol = u.kg/u.gm*u.mol;
            u.mmol = u.mol/1000;
            u.umol = u.mol/1e6;
            
            %------- Volume -------
            u.cc = (u.cm)^3;           
            u.L = 1000*u.cc;           
            u.mL = u.cc;               
            u.floz = 29.5735297*u.cc;  
            u.pint = 473.176475*u.cc;  
            u.quart = 946.35295*u.cc;  
            u.gal = 3.78541197*u.L;    
            
            %---- frequency ----
            u.Hz = 1/u.s;       
            u.kHz = 1e3 *u.Hz;  
            u.MHz = 1e6 *u.Hz;
            u.GHz = 1e9 *u.Hz;
            
            %---- force -------
            u.N = u.kg*u.m/u.s^2;   
            u.dyne = 1e-5*u.N;      
            u.lbf = 4.44822*u.N;    
            
            %----- energy -----
            u.J = u.kg*u.m^2/u.s^2; 
            u.MJ = 1e6*u.J;         
            u.kJ = 1e3*u.J;         
            u.mJ = 1e-3*u.J;        
            u.uJ = 1e-6*u.J;        
            u.nJ = 1e-9*u.J;        
            u.eV = 1.6022e-19*u.J;    
            u.BTU = 1.0550559e3*u.J;  
            u.kWh = 3.6e6*u.J;        
            u.cal = 4.1868*u.J;       
            u.kcal = 1e3*u.cal;       
            u.erg = 1e-7*u.J;
            
            %---- pressure -----
            u.Pa = u.N/u.m^2;
            u.kPa = 1000*u.Pa;
            u.MPa = 1e6*u.Pa;
            u.GPa = 1e9*u.Pa;
            u.torr = 133.322*u.Pa;
            u.mtorr = 1e-3*u.torr;
            u.bar = 1e5*u.Pa;
            u.mbar = 1e-3*u.bar;
            u.atm = 1.013e5*u.Pa;
            u.psi = 6.895e3*u.Pa;
            u.mmHg = 1/760*u.atm;
            
            %----- power --- ---
            u.W = u.J/u.s;
            u.MW = 1e6*u.W;
            u.kW = 1e3*u.W;
            u.mW = 1e-3*u.W;
            u.uW = 1e-6*u.W;
            u.nW = 1e-9*u.W;
            u.pW = 1e-12*u.W;
            u.hp = 745.69987*u.W;
            
            %------ Voltage -----
            u.V = u.J/u.coul;
            u.kV = 1e3*u.V;
            u.mV = 1e-3*u.V;
            u.uV = 1e-6*u.V;
            
            %----- Current ------
            u.A = u.coul/u.s;
            u.mA = 1e-3*u.A;
            u.uA = 1e-6*u.A;
            u.nA = 1e-9*u.A;
            
            %----magnetic field -----
            u.T = u.V*u.s/u.m^2;
            u.tesla = u.T;
            
            u.gauss = 1e-4*u.T;
            
            %----area----------------
            u.acre = 4840*u.yd^2;
            u.hectare = 10000*u.m^2;
            
            %----electromagnetic units-----
            u.ohm = u.V/u.A;
            u.H = u.ohm*u.s; % Henry
            u.Wb = u.V*u.s; % Weber
            u.S = 1/u.ohm; % siemens
            u.siemens = u.S;
            
            u.F = u.coul/u.V; % farad
            u.farad = u.F;
            
            % now we set the display name for each object in the structure
            % displaystring is the same as the field name of the structure
            fnames = fieldnames(u);
            for i=1:length(fnames)
                fn = fnames{i};
                s = sprintf('u.%s.displaystring = ''%s'';',fn,fn);
                eval(s);
            end
            
            u.degC = @cmu.unit.degC;
            u.degF = @cmu.unit.degF;
            u.degF2C = @cmu.unit.degF2C;
            u.degC2F = @cmu.unit.degC2F;
            u.degF2R = @cmu.unit.degF2R;
            u.degC2R = @cmu.unit.degC2R;
            
            %----dimensionless-----------
            % sometimes you need this to force a number to be a unit
            % class so the right function will be called from the unit
            % class, e.g. when the limits for integration in quad are
            % dimensionless, or for an ode, or for a nonlinear solver.
            u.dimensionless = cmu.unit(1,[0 0 0 0 0 0 0],'');
        end
        
        function u = simple_units(base)
            % returns only the structure of conversion factors as doubles.
            % useful for functions where unit objects do not work. However,
            % unit algebra is not followed. 
            if nargin == 1
                u = cmu.unit.units(base);
            else
                u = cmu.unit.units;
            end
            
            fn = fieldnames(u);
            for i = 1:length(fn)
                if strcmp(fn{i},'degC') || strcmp(fn{i},'degF')
                    cmd = sprintf('u.%s = 1;',fn{i});
                    eval(cmd);
                elseif strcmp(fn{i},'degC2F') || strcmp(fn{i},'degF2C')
                    continue
                elseif strcmp(fn{i},'degC2R') || strcmp(fn{i},'degF2R')
                    continue
                else
                    cmd = sprintf('u.%s = double(u.%s);',fn{i},fn{i});
                    eval(cmd);
                end
            end
        end
        
        function c = constants(base)
            % common constants with units. the default base units are the
            % same as the base units for cmu.unit.units.
            %
            % c = cmu.unit.constants;
            % c.R  % the gas constant
            if nargin == 1
                u = cmu.unit.units(base);
            else
                u = cmu.unit.units;
            end
            
            c.R = 8.314472*u.J/(u.mol*u.K); % gas constant
            c.c = 299792458*u.m / u.s;      % speed of light
            c.Na = u.gm/u.amu;              % Avagadro's number
        end
        
        % temperature is especially problematic as a unit because the
        % conversions between units involve additive constants. for now we
        % provide these static functions that just convert Celcius and
        % Fahrenheit to Kelvin.
        function K = degC(C)
           % Celcius to Kelvin
           % >> u = cmu.units;
           % >> T1 = u.degC(100)
           u = cmu.unit.units;
           K = (C + 273.15)*u.K;
        end
       
        function K = degF(F)
           %Fahrenheit to Kelvin
           % >> u = cmu.units;
           % >> T1 = u.degF(100)
           C = cmu.unit.degF2C(F);
           K = cmu.unit.degC(C);
        end
       
        function C = degF2C(F)
           %Fahrenheit to Celcius
           % no units attached to output!
           C = (F - 32)*5/9;
       end
       
       function F = degC2F(C)
           %Celsius to Fahrenheit
           % no units attached to output!
           F = C*9/5 + 32;
       end
             
       function R = degF2R(F)
           % Fahrenheit to Rankine
           % u = cmu.units;
           % R = degF2R(212) % 212F in Rankine
           R = (F + 459.67)*u.R;
       end
       
       function R = degC2R(C)
           % Celcius to Rankine
           % >> u = cmu.units;
           % >> R = u.degC2R(100)  %100 degC in Rankine
           K = cmu.unit.degC(C);
           R = 5/9*double(K)*u.R;
       end
       
    end
end