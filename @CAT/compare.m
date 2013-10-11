function [results,same,different] = compareCATs(CAT1,CAT2,varargin)

p = properties(CAT1);
fprintf('Property-->Result \n-------------------\n')
result = zeros(length(p),1);
for i = 1:length(p)
    
    if isa(CAT1.(p{i}),'function_handle') && isa(CAT2.(p{i}),'function_handle')
        str1 = func2str(CAT1.(p{i}));
        str2 = func2str(CAT1.(p{i}));
    elseif isa(CAT1.(p{i}),'Distribution') && isa(CAT2.(p{i}),'Distribution')
        str1 = Dist2str(CAT1.(p{i}));
        str2 = Dist2str(CAT1.(p{i}));
    end
        result(i) = isequal(str1,str2);
    
    
    if nargout == 0
        if result(i)
            resultstr = 'EQUAL';
        else
            resultstr = 'DIFFERENT';
        end

        fprintf(strcat(p{i},'-->',resultstr,'\n'))
    end
    
end

if nargout>0
    results = result;
    issame = 1;isdiff = 1;
    same = [];different = [];
    for i = 1:length(result)
        if result(i)
            same{issame} = p{i};
            issame = issame +1;
        else
            different{isdiff} = p{i};
            isdiff = isdiff + 1;
        end
    end
end
