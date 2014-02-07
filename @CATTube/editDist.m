%% Method editDist

function editDist(O,~,~)

% Open new GUI window which lets user define the distribution
% by choosing a grid and

% Calculate settings for currently defined distribution
if allvaleq( diff( log10(O.init_dist.y) ) )
    % Log spacing
    orig_gmin = min(log10(O.init_dist.y));
    orig_gmax = max(log10(O.init_dist.y));
    orig_gspacing = mean(diff( log10( O.init_dist.y )));
else
    % Linear spacing
    orig_gmin = min(O.init_dist.y);
    orig_gmax = max(O.init_dist.y);
    orig_gspacing = mean(diff( O.init_dist.y ));
end % if else

orig_gnumpoints = length(O.init_dist.y);

% Get distribution density definition
orig_densfnc = O.init_dist.getFunction();
orig_densval = O.init_dist.F;
if ~isempty(orig_densfnc)
    % Cut @(x) from front of definition
    orig_densfnc = func2str(orig_densfnc);
    orig_densfnc = orig_densfnc(5:end);
    orig_densval = [];
    % Set active field
    Dgui.density.activefield = 'function';
elseif ~isempty(orig_densval)
    orig_densfnc = '';
    orig_densval = O.init_dist.F;
    % Set active field
    Dgui.density.activefield = 'values';
else
    % No distribution defined
    orig_densfnc = '';
    orig_densval = [];
    % Set active field to function by default
    Dgui.density.activefield = 'function';
end % if else
    
% Create the figure handle
Dgui.fighandle = figure(...
    'MenuBar','none',...
    'Name','Edit initial distribution',...
    'NumberTitle','off',...
    'Position',[200 200 500 445],...
    'Resize','off');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Grid
%

% Panel
Dgui.grid.panel = uipanel('Parent',Dgui.fighandle,...
    'Title','Grid',...
    'Units','pixels',...
    'Position',[20 355 460 70]);

% Button group for the radio buttons
Dgui.grid.spacingtype = uibuttongroup('Parent',Dgui.grid.panel,...
    'Title','',...
    'Units','pixels',...
    'BorderType','none',...
    'Position',[10 10 100 40],...
    'SelectionChangeFcn',@(hObject,Eventdata)changegridtype(hObject,Eventdata));

% Radio button 'linear'
Dgui.grid.spacingtype_lin = uicontrol(Dgui.grid.spacingtype,...
    'Style','radiobutton',...
    'String','Linear',...
    'Position',[0 20 60 20]);

% Radio button 'logarithmic'
Dgui.grid.spacingtype_lin = uicontrol(Dgui.grid.spacingtype,...
    'Style','radiobutton',...
    'String','Log10',...
    'Position',[0 0 60 20]);

% Opening bracket
Dgui.grid.bracket_left = uicontrol(Dgui.grid.panel,...
    'Style','text',...
    'String','[',...
    'Units','pixels',...
    'FontSize',20,...
    'HorizontalAlignment','right',...
    'Position',[80 10 60 30]);

% Min field
Dgui.grid.min_text = uicontrol(Dgui.grid.panel,...
    'Style','text',...
    'String','Min.',...
    'Units','pixels',...
    'Position',[140 35 70 20]);

Dgui.grid.min = uicontrol(Dgui.grid.panel,...
    'Style','edit',...
    'String',num2str(orig_gmin),...
    'TooltipString','Starting point of calculation grid',...
    'Units','pixels',...
    'Callback',@(hObject,Eventdata)changegrid(hObject,Eventdata,'min'),...
    'Position',[140 10 70 30],...
    'BackgroundColor','w');

% Spacing field
Dgui.grid.spacing_text = uicontrol(Dgui.grid.panel,...
    'Style','text',...
    'String','Spacing',...
    'Units','pixels',...
    'Position',[215 35 70 20]);

Dgui.grid.spacing = uicontrol(Dgui.grid.panel,...
    'Style','edit',...
    'String',num2str(orig_gspacing),...
    'Enable','off',...
    'TooltipString','Spacing of grid points',...
    'Units','pixels',...
    'Position',[215 10 70 30]);

% Max field
Dgui.grid.max_text = uicontrol(Dgui.grid.panel,...
    'Style','text',...
    'String','Max.',...
    'Units','pixels',...
    'Position',[290 35 70 20]);

Dgui.grid.max = uicontrol(Dgui.grid.panel,...
    'Style','edit',...
    'String',num2str(orig_gmax),...
    'TooltipString','End point of calculation grid',...
    'Units','pixels',...
    'Callback',@(hObject,Eventdata)changegrid(hObject,Eventdata,'max'),...
    'Position',[290 10 70 30],...
    'BackgroundColor','w');

% Closing bracket
Dgui.grid.bracket_right = uicontrol(Dgui.grid.panel,...
    'Style','text',...
    'String',']',...
    'Units','pixels',...
    'FontSize',20,...
    'HorizontalAlignment','left',...
    'Position',[360 10 30 30]);

% Number of gridpoints field
Dgui.grid.numpoints_text = uicontrol(Dgui.grid.panel,...
    'Style','text',...
    'String','Gridpoints',...
    'Units','pixels',...
    'Position',[380 35 70 20]);

Dgui.grid.numpoints = uicontrol(Dgui.grid.panel,...
    'Style','edit',...
    'String',num2str(orig_gnumpoints),...
    'TooltipString','Number of grid points',...
    'Units','pixels',...
    'Callback',@(hObject,Eventdata)changegrid(hObject,Eventdata,'numpoints'),...
    'Position',[380 10 70 30],...
    'BackgroundColor','w');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Population density
%

% Panel
Dgui.density.panel = uipanel('Parent',Dgui.fighandle,...
    'Title','Population density',...
    'Units','pixels',...
    'Position',[20 275 460 70]);

% Function field
Dgui.density.function_text = uicontrol(Dgui.density.panel,...
    'Style','text',...
    'String','Function',...
    'Units','pixels',...
    'Position',[40 35 180 20]);
% @(x) label
Dgui.density.function_text2 = uicontrol(Dgui.density.panel,...
    'Style','text',...
    'String','@(x)',...
    'Units','pixels',...
    'Position',[10 10 30 20]);
Dgui.density.function = uicontrol(Dgui.density.panel,...
    'Style','edit',...
    'String',orig_densfnc,...
    'TooltipString','Function defining population density',...
    'Units','pixels',...
    'Callback',@(hObject,Eventdata)changedensity(hObject,Eventdata,'function'),...
    'Position',[40 10 205 30]);

% Values field
% Text
Dgui.density.values_text = uicontrol(Dgui.density.panel,...
    'Style','text',...
    'String','Values',...
    'Units','pixels',...
    'Position',[255 35 180 20]);
% Left bracket
Dgui.density.values_bracketleft = uicontrol(Dgui.density.panel,...
    'Style','text',...
    'String','[',...
    'Units','pixels',...
    'FontSize',20,...
    'HorizontalAlignment','right',...
    'Position',[255 10 5 30]);
% Input box
Dgui.density.values = uicontrol(Dgui.density.panel,...
    'Style','edit',...
    'String',regexprep(data2str(orig_densval),'\[|\]',''),...
    'TooltipString','Values defining population density (separate with space). You can use workspace variables as well.',...
    'Units','pixels',...
    'Callback',@(hObject,Eventdata)changedensity(hObject,Eventdata,'values'),...
    'Position',[260 10 180 30]);
% Right bracket
Dgui.density.values_bracketright = uicontrol(Dgui.density.panel,...
    'Style','text',...
    'String',']',...
    'Units','pixels',...
    'FontSize',20,...
    'HorizontalAlignment','left',...
    'Position',[440 10 10 30]);

% Set white background for currently active field
set(Dgui.density.(Dgui.density.activefield),'BackgroundColor','w')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preview
%

% Panel
Dgui.preview.panel = uipanel('Parent',Dgui.fighandle,...
    'Title','Preview',...
    'Units','pixels',...
    'Position',[20 45 460 220]);

% Update button
Dgui.preview.update = uicontrol(Dgui.preview.panel,...
    'Style','pushbutton',...
    'String','Update',...
    'Value',0,...
    'Callback',@plotDist,...
    'Position',[405 185 50 20]...
    );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OK, Reset, Cancel buttons
%

% OK Button
Dgui.buttons.ok = uicontrol(Dgui.fighandle,...
    'Style','pushbutton',...
    'String','OK',...
    'Value',0,...
    'Callback',@returnVals,...
    'Position',[70 10 100 30]...
    );

% Reset Button
Dgui.buttons.ok = uicontrol(Dgui.fighandle,...
    'Style','pushbutton',...
    'String','Reset',...
    'Value',0,...
    'Callback',@resetVals,...
    'Position',[200 10 100 30]...
    );

% Cancel Button
Dgui.buttons.ok = uicontrol(Dgui.fighandle,...
    'Style','pushbutton',...
    'String','Cancel',...
    'Value',0,...
    'Callback',@cancelDgui,...
    'Position',[330 10 100 30]...
    );

%% % - Subfunctions

%% Function cancelDgui

    function cancelDgui(~,~)
        
        % Close the editing GUI without doing anything
        close(Dgui.fighandle)
        
        clear Dgui
        
    end % function cancelDgui

%% Function changedensity

    function changedensity(hObject,~,field)
        
        % Change background of the field not currently edited
        % Save active field - function or values
        
        % Set both background to default, then change current field to
        % white
        set([Dgui.density.function Dgui.density.values],'BackgroundColor','default');
        set(hObject,'BackgroundColor','w')
        
        % Save active field
        Dgui.density.activefield = field;
        
    end % function

%% Function changegrid

    function changegrid(hObject,~,field)
        
        % Update spacing or gridpoints field when min, spacing, max or
        % gridpoints field changed
        
        newval = str2double(get(hObject,'String'));
        
        switch field
            
            case 'min'
                % New minimum point
                
                gmax = str2double(get(Dgui.grid.max,'String'));
                gnumpoints = str2double(get(Dgui.grid.numpoints,'String'));
                gmin = newval;
                
            case 'max'
                % New maximum point
                
                gmin = str2double(get(Dgui.grid.min,'String'));
                gnumpoints = str2double(get(Dgui.grid.numpoints,'String'));
                gmax = newval;
                
            case 'numpoints'
                % New number of points
                
                gmin = str2double(get(Dgui.grid.min,'String'));
                gmax = str2double(get(Dgui.grid.max,'String'));
                gnumpoints = newval;
                
        end % switch
        
        newspacing = (gmax - gmin) / (gnumpoints-1);
        set(Dgui.grid.spacing,'String',num2str(newspacing));
        
    end % function changespacing

%% Function changegridtype

    function changegridtype(hObject,Eventdata)
        
        % Change text infront of bracket
        
        if strcmp( get(Eventdata.NewValue,'String') , 'Linear')
            set(Dgui.grid.bracket_left,'String','[')
            set(Dgui.grid.bracket_right,'String',']')
            
            % Recalculate min, max, spacing
            gmin = 10^str2double(get(Dgui.grid.min,'String'));
            if ~isfinite(gmin)
                gmin = 0;
            end % if
            gmax = 10^str2double(get(Dgui.grid.max,'String'));
            gnumpoints = str2double(get(Dgui.grid.numpoints,'String'));
            gspacing = (gmax - gmin) / (gnumpoints-1);
        else
            set(Dgui.grid.bracket_left,'String','10^(')
            set(Dgui.grid.bracket_right,'String',')')
            
            % Recalculate min, max, spacing
            gmin = log10(str2double(get(Dgui.grid.min,'String')));
            if isinf(gmin)
                gmin = [];
            end % if
            gmax = log10(str2double(get(Dgui.grid.max,'String')));
            gnumpoints = str2double(get(Dgui.grid.numpoints,'String'));
            gspacing = (gmax - gmin) / (gnumpoints-1);
        end % if else
        
        % Set new values
        set(Dgui.grid.min,'String',num2str(gmin));
        set(Dgui.grid.spacing,'String',num2str(gspacing));
        set(Dgui.grid.max,'String',num2str(gmax));
        
    end % function

%% Function getDgui_current

    function [x,f] = getDgui_current(type)
        
        % Return currently defined x vector and f values or function, if
        % defined.
        %
        % [x,f] = getDgui_current('func')
        % returns f as a function, if defined (default)
        %
        % [x,f] = getDgui_current('values')
        % returns f as values evaluated at x
        
        if nargin < 1 || isempty(type)
            type = 'func';
        end % if
        
        % Get currently defined x value
        % Use min:spacing:max, assume that number of gridpoints has been
        % calculated correctly
        xmin = str2double( get(Dgui.grid.min,'String') );
        xmax = str2double( get(Dgui.grid.max,'String') );
        xnumpoints = str2double( get( Dgui.grid.numpoints,'String') );
        
        if strcmp( get(get(Dgui.grid.spacingtype,'SelectedObject'),'String') , 'Log10' )
            x = logspace(xmin,xmax,xnumpoints);
        else
            x = linspace(xmin,xmax,xnumpoints);
        end % if
        
        % Get currently defined distribution function
        
        % Check for defined function
        if strcmp(Dgui.density.activefield,'function')
            
            fncstr = get(Dgui.density.function,'String');
            if ~isempty(fncstr)
                f = str2func([get(Dgui.density.function_text2,'String') fncstr]);
            else
                f = [];
            end % if else
            
            % If values wanted, translate to values
            if strcmpi(type,'values')
                try
                    f = f(x);
                catch ME
                    warndlg(...
                        sprintf('Distribution function could not be evaluated. Error: %s',ME.message),...
                        'Error evaluating given function','modal');
                end % try-catch
            end % if
            
        else
            
            % Get values field
            vals = get(Dgui.density.values,'String');
            
            try
                f = evalin('base',[ '[' vals ']' ]);
                if length(f) ~= length(x)
                    warndlg('Vector of values must have the right number of entries',...
                    'Error using given values','modal');
                    f = zeros(size(x));
                end % if
            catch ME
                warndlg(...
                    sprintf('Distribution values could not be applied. Error: %s',ME.message),...
                    'Error using given values','modal');
            end % try-catch
            
        end % if else
        
    end % function

%% Function plotDist

    function plotDist(~,~)
        
        % Check for the plot axis, create if it does not exist
        if ~isfield(Dgui.preview,'axes') || ~ishandle(Dgui.preview.axes)
            Dgui.preview.axes = axes('Parent',Dgui.preview.panel,...
                'Position',[0.13 0.19 0.75 0.73]);
            xlabel(Dgui.preview.axes,'x')
            ylabel(Dgui.preview.axes,'f(x)')
            box(Dgui.preview.axes,'on')
        else
            % If axes exists, clear any previous plots
            if isfield(Dgui.preview,'lines') && ishandle(Dgui.preview.lines)
                delete(Dgui.preview.lines)
            end % if
        end % if
        
        % Create new plot
        [x,f] = getDgui_current('values');
        
        Dgui.preview.lines = line(x,f);
        
    end % function plotDist

%% Function returnVals

    function returnVals(~,~)
        
        % Set values to distribution, close GUI
        
        % Get currently defined values
        [x,f] = getDgui_current;
        % Set these values
        O.init_dist = Distribution(x,f);
        
        % Update the field for the current variable
        set(O.gui.init.init_dist,'String',data2str(O.init_dist));
        
        % Close the GUI
        close(Dgui.fighandle)
        
        clear Dgui
        
    end % function

%% Function resetVals

    function resetVals(~,~)
        
        % Reset fields to the values which were valid when the window was
        % opened
        
        set(Dgui.grid.min,'String',num2str(orig_gmin));
        set(Dgui.grid.max,'String',num2str(orig_gmax));
        set(Dgui.grid.spacing,'String',num2str(orig_gspacing));
        set(Dgui.grid.numpoints,'String',num2str(orig_gnumpoints));
        
        set(Dgui.density.function,'String',orig_densfnc);
        set(Dgui.density.values,'String',regexprep(data2str(orig_densval),'\[|\]',''));
        
    end % function resetVals

end % function editDist
