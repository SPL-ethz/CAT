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

h=get(gcf,'children');
axes('Parent',findall(h,'bordertype','none'),...
    'Position',[0.13 0.19 0.75 0.53],'xtick',[],'ytick',[]);
set(gca,'tag','1')
xlim([0 1])
axis off
hold on

drawnow

% Run the real solve method
solve@CAT(O);

% Remove text again
delete(O.gui.run.solving_text)
delete(gca)

% Activate the plot button
set(O.gui.run.plot,'Enable','on');

end % function