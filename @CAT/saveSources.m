%% -- saveSources

function saveSources(O,sourceStruct)

if nargin<2
    sourceStruct = [];
end

fieldnames = properties(O);
fieldnames(strcmp(fieldnames,'gui')) = [];
fieldnames(~cellfun(@isempty,strfind(fieldnames,'calc'))) = [];

namestr = strcat('CATsource','_',datestr(now,'YYYYmmdd_HHMMSS'),'.m');
fid = fopen(namestr,'w+');

fprintf(fid,'kitty = CAT; \n\n');

for i = 1:length(fieldnames)
    y = explode(help(strcat('CAT.',fieldnames{i})),{'\n','\r'});
    pcsstr = repmat({'%'},[length(y')-2 1]);
    pcsstr = ['%%';pcsstr];
    y = [pcsstr(:) y(1:end-1)'];
    
    for j = 1:length(y(:,1))
        fprintf(fid,'%s %s \n',y{j,:});
    end
    
    if ~strcmp(fieldnames{i},'init_dist')
        if ~isempty(sourceStruct)  && isfield(sourceStruct,fieldnames{i})
            valstr = sourceStruct.(fieldnames{i});
        else
            valstr = data2str(O.(fieldnames{i}));
            if ischar(O.(fieldnames{i}))
                valstr = ['''',valstr,''''];
            end
        end
        
        fprintf(fid,strcat('kitty.',fieldnames{i},' = ',valstr,'; \n\n\n'));
    else
        if ~isempty(sourceStruct) && isfield(sourceStruct,'init_dist') && isfield(sourceStruct.init_dist,'y')
            valstr{1} = sourceStruct.init_dist.y;
        else
            if isa(O.init_dist,'Distribution')
                valstr{1} = data2str(O.init_dist.y);
            else
                valstr{1} = '[]';
            end
        end
        
        if ~isempty(sourceStruct) && isfield(sourceStruct,'init_dist') && isfield(sourceStruct.init_dist,'F')
            valstr{2} = sourceStruct.init_dist.F;
        else
            if isa(O.init_dist,'Distribution')
                valstr{2} = data2str(O.init_dist.F);
            else
                valstr{2} = '[]';
            end
        end
        
        fprintf(fid,strcat('kitty.init_dist = Distribution(',valstr{1},',',valstr{2},');\n\n\n'));
    end
    
end

fclose(fid);

end