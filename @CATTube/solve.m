function O = solve(O)

% CATTube.gui_solve

% Disable 'Solve' button
set(O.gui.run.start,'Enable','off')

% Run the real solve method
solve@CAT(O);

% Activate the plot button and the overwrite checkbox
set(O.gui.run.plot,'Enable','on');
set(O.gui.run.plot_overwrite,'Enable','on');

end % function