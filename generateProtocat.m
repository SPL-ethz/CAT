function generateProtocat(kitsetup)

for i = 1:length(kitsetup)
    if isa(kitsetup{i},'function_handle')
        kitsetup{i} = func2str(kitsetup{i});
    elseif isnumeric(kitsetup{i})
        kitsetup{i} = num2str(kitsetup{i});
    end
end

n = 1;
while exist(strcat('test',num2str(n)),'file')
    n = n+1;
end
fidnew = fopen(strcat('test',num2str(n),'.m'),'w+');
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
        tnew = '\n %% This file was generated for you.';
    end

    fprintf(fidnew,tnew);

end

fclose(fid);
fclose(fidnew);

% fprintf(strcat('%%%% PROTOCAT %%%% \n \n %% This is an empty template for setting up and running CAT problems \n\n kitty = CAT; %% define kitty as CAT object',...
%     '\n \n %%-------------------------------------------------------------------------- \n ',...
%     '%%%% THERMODYNAMICS & KINETICS \n %% Solubility \n %% FIELD:    solubility \n %% STATUS:   REQUIRED',...
%     '%% UNITS:    [g/g] \n %% CLASS:    FUNCTION_HANDLE \n %% INPUTS:   1. Temperature [°C],  2.(optional) Mass fraction antisolvent [-] \n',...
%     '%% EXP.:     @(T,xm) 1e-4*(1-xm)+(0.0056*(T).^2+0.0436.*(T)+3.7646)/1000 \n %% COMMENT:  - \n',...
%     'kitty.solubility = %s'),kitsetup{1});

