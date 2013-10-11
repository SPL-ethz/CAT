function fillOutForm(kitsetup)

for i = 1:length(kitsetup)
    if isa(kitsetup{i},'function_handle')
        kitsetup{i} = func2str(kitsetup{i});
    elseif isnumeric(kitsetup{i})
        kitsetup{i} = mat2str(kitsetup{i});
    elseif isa(kitsetup{i},'Distribution')
        kitsetup{i} = Dist2str(kitsetup{i});
    end
end

n = 1;
while exist(strcat('CAT_form',num2str(n)),'file')
    n = n+1;
end
fidnew = fopen(strcat('CAT_form',num2str(n),'.m'),'w+');
fid = fopen('protocat.m');


tline = fgets(fid);
i = 1;
while ischar(tline)
    disp(tline)
    tline = fgets(fid);
    
    if ischar(tline) && ~isempty(strfind(tline,'XXX'))
        if isempty(kitsetup{i})
            tnew = strrep(tline, tline, strcat('%% ',tline));
        else
            tnew = strrep(tline, 'XXX', kitsetup{i});
        end
        i = i + 1;
    elseif ischar(tline)
        tnew = tline;
    else
        tnew = '\n \n %% This file was generated for you by the fillOutForm function.';
    end

    fprintf(fidnew,tnew);

end

fclose(fid);
fclose(fidnew);
