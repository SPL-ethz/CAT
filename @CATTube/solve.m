function O = solve(O)

% CATTube.gui_solve

% Disable 'Solve' button
set(O.gui.run.start,'Enable','off')

% Show something graphical?
O.gui.run.solving_text = uicontrol(O.gui.run.panel,...
    'Style','text',...
    'String','Solving...',...
    'Position',[20 50 450 30],...
    'Fontsize',12,...
    'FontWeight','bold',...
    'ForegroundColor','r'...
    );

% Force drawnow to actually display the text..
drawnow

% Run the real solve method
solve@CAT(O);

% Remove text again
delete(O.gui.run.solving_text)

% Activate the plot button
set(O.gui.run.plot,'Enable','on');

end % function