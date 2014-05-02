%% Method browseVars

function browseVars(O,~,~,classvarname,displayfield,classfilter)

% BrowseVars is a callback function called when a browse button is pressed
% next to a variable input box. It lists the variables in the workspace and
% allows the choice of one
% Assign the chosen variable to the variable in the class
if nargin < 6 || isempty(classfilter)
    classfilter = {};
else
    % Check if classfilter is char or cell - convert to cell
    if ischar(classfilter)
        classfilter = {classfilter};
    else
        classfilter = [classfilter];
    end % if
end % if else

% Create a new window for the list
glb.fighandle = figure(...
    'MenuBar','none',...
    'Name','Import a variable',...
    'NumberTitle','off',...
    'Position',[200 200 400 400],...
    'Resize','off');

% Create a list box
glb.lbox = uicontrol(glb.fighandle,...
    'Style','listbox',...
    'String','',...
    'Max',1,...
    'Min',0,...
    'Units','pixels',...
    'Callback',@(hObject,Event)varDetails(hObject,Event),...
    'Position',[10 10 260 380]);

% Update its contents
glb.V = getVarList;
updateVarList();

% Create import button
glb.import = uicontrol(glb.fighandle,...
    'Style','pushbutton',...
    'String','Import',...
    'Units','pixels',...
    'Callback',@(hObject,Event)setVar(hObject,Event),...
    'Position',[280 360 110 30]);

% Create text box for variable details
glb.text = uicontrol(glb.fighandle,...
    'Style','text',...
    'String','',...
    'Units','pixels',...
    'BackgroundColor',get(glb.fighandle,'Color'),...
    'Position',[280 320 110 30]);

% Create refresh button
glb.import = uicontrol(glb.fighandle,...
    'Style','pushbutton',...
    'String','Refresh',...
    'Units','pixels',...
    'Callback',@(hObject,Event)updateVarList(hObject,Event),...
    'Position',[280 260 110 30]);

    function V = getVarList
        
        % Get variable list
        varList = evalin('base','who');
        varList = varList(~strcmp('varList',varList) & ~strcmp('activeList',varList));
        varList = [varList repmat({''},size(varList))];
        activeList = varList;
        % if there are structures, whose fields may contain the necessary
        % data, we should be aware of this... therefore let's first explore
        % all structures in the workspace
        while ~isempty(activeList)
            if isa(evalin('base',[activeList{1,2},activeList{1,1}]),'struct')
                for i = 1:length(evalin('base',[activeList{1,2},activeList{1,1}]))
                    if length(evalin('base',[activeList{1,2},activeList{1,1}]))>1
                        indstr = ['(',num2str(i),')'];
                    else
                        indstr = '';
                    end
                    fields = fieldnames(evalin('base',[activeList{1,2},activeList{1,1},indstr]));
                    varList = [varList; fields repmat({[activeList{1,2},activeList{1,1},indstr,'.']},[length(fields) 1])];
                    activeList = [activeList; fields repmat({[activeList{1,2},activeList{1,1},indstr,'.']},[length(fields) 1])];
                end
            end

            activeList(1,:) = [];

        end
        
        
        % If filters defined, keep only those class types listed
        if ~isempty(classfilter)
            Vname = cell(0);
            for i = 1:length(classfilter)
                for j = 1:length(varList(:,1))
                   if isa(evalin('base',[varList{j,2},varList{j,1}]),classfilter{i})
                       Vname = [Vname;[varList{j,2},varList{j,1}]];
                   end
                end
            end % for
        else
            Vname = cell(0);
            for j = 1:length(varList(:,1))
                Vname = [Vname;[varList{2,j},varList{1,j}]];
            end
        end % if
        
        for i = 1:length(Vname)
           V(i).name = Vname{i};
           V(i).class = class(evalin('base',Vname{i}));
           V(i).size = size(evalin('base',Vname{i}));
        end
        
        
    end % function getVarList

    function updateVarList(~,~)
        
        glb.V = getVarList;
        
        
        varnames = [{glb.V.name}']; %#ok<NBRAK>
        
        set(glb.lbox,'String',varnames);
        
    end % function updateVarList

    function varDetails(hObject,~)
        
        % Find selected variable
        varnum = get(hObject,'Value');
        
        % Make string for description of this variable
        % Output something like:
        % a
        % 2x3 double
        
        vname = sprintf('%s\n',glb.V(varnum).name);
        if length(glb.V(varnum).size) == 2
            vsize = sprintf('%ix%i',glb.V(varnum).size);
        else %multidimensional
            vsize = sprintf('%ix%ix...',glb.V(varnum).size(1:2));
        end % if else
        
        vtext = sprintf('%s %s %s',vname,vsize,glb.V(varnum).class);
        
        % Print variable details
        set(glb.text,'String',vtext)
        
%         if (strcmp(get(glb.fighandle, 'SelectionType'), 'open')) % if double click
%             expandStruct(hObject,[])
%         end
        
    end % function

%     function expandStruct(hObject,~)
%         strind = get(glb.lbox,'value');
%         strname = structInfos.names{structInfos.values == strind};
%         fields = fieldnames(evalin('base',strname));
%         updateVarList(hObject,[],fields);
%         keyboard
%     end

    function setVar(hObject,~)
        
        % Get the chosen variable number
        varnum = get(glb.lbox,'Value');
        
        % Assign the corresponding CAT variable to the chosen variable
        % Get variable data
        if varnum <= length(glb.V)
            
            if ~isempty(classvarname)
                vardata = evalin('base',[glb.V(varnum).name]);
                set(O,classvarname,vardata);
            end % if
            
            % Close the variable list window
            close(get(hObject,'Parent'))
            
            % Update the field for the current variable
            set(displayfield,'String',data2str(vardata));
            
        end % if
        
    end % function

end % function