function outstr = data2str(vdata)

% data2str( vdata )
% Smart transform of variable data to string representation - tries its
% best to make the shortest possible representation of the variable vdata
%

% Make string representation of the data
switch(class(vdata))
    
    case 'char'
        % Do nothing to str inputs - already formatted
        outstr = vdata;
        
    case 'function_handle'
        
        outstr = func2str(vdata);
        
    case {'double'}
        
        % Numbers - can be scalar or vectors, matrices
        if isscalar(vdata)
            outstr = sprintf('%.15g',vdata);
        elseif isempty(vdata)
            outstr = '[]';
        else
            
            % Check whether vector representation can be shortened
            if isvector(vdata) && ...
                    allvaleq( diff(vdata) ) && ...
                    all( diff(vdata) > 0 )
                % Currently looking at a vector, where all differences
                % are equal and positive: linspace
                startval = vdata(1);
                endval = vdata(end);
                numvals = length(vdata);
                outstr = sprintf('linspace(%g,%g,%i)', startval,endval,numvals);
            elseif isvector(vdata) && ...
                    allvaleq( diff( log10(vdata) ) ) && ...
                    all( diff( log10(vdata) ) > 0 )
                % Check for log-spaced vector
                startval = log10(vdata(1));
                endval = log10(vdata(end));
                numvals = length(vdata);
                outstr = sprintf('logspace(%g,%g,%i)',startval,endval,numvals);
            elseif vdata(1) ~= mean(vdata(:)) && ...
                    isvector(vdata) && length(vdata) > 2 && ...
                    allvaleq( diff( log10(vdata(2:end)) ) ) && ...
                    all( diff( log10(vdata(2:end)) ) > 0 )
                % Check for log-spaced vector with a different first
                % value! (Used in definition of time vector sometimes)
                firstval = vdata(1);
                startval = log10(vdata(2));
                endval = log10(vdata(end));
                numvals = length(vdata)-1;
                outstr = sprintf('[%g logspace(%g,%g,%i)]',...
                    firstval,startval,endval,numvals);
            elseif vdata(1) ~= mean(vdata(:)) && ...
                    isvector(vdata) && length(vdata) > 2 && ...
                    allvaleq( diff( vdata(2:end) ) ) && ...
                    all( diff( vdata(2:end) ) > 0 )
                % Check for lin-spaced vector with a different first
                % value! (Used in definition of time vector sometimes)
                firstval = vdata(1);
                startval = vdata(2);
                endval = vdata(end);
                numvals = length(vdata)-1;
                outstr = sprintf('[%g linspace(%g,%g,%i)]',...
                    firstval,startval,endval,numvals);
            elseif all(vdata == 0)
                % Check if all values are equal to zero
                outstr = sprintf( 'zeros(%i,%i)',size(vdata) );
            elseif length(vdata) > 1 && allvaleq( vdata )
                % All values are equal - print short form
                outstr = sprintf( '%g*ones(%i,%i)', mean(vdata(:)), size(vdata) );
            else
                % Not a regular vector, print whole matrix/vector
                outstr = mat2str(vdata);
            end % if else
            
        end % if else
        
    otherwise
        % For everything which is not recognized, try mat2str
        outstr = mat2str(vdata);
        
end % switch

end % function