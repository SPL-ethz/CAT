classdef CATTube < CAT
    
    % GUI for CAT
    %
    %
    
    % Martin Iggland, first version 2013-12-04
    
    properties
        
        % GUI handles and stuff
        gui
        
    end % properties
    
    methods
        
        %% Constructor method
        
        function O = CATTube(varargin)
            
            % Call CAT constructor
            O = O@CAT( varargin{:} );
            
            % Start constructing GUI
            
            % Initialise the main GUI, define all buttons and graphics
            % objects
            
            %
            % Initialise the GUI
            %
            % Create the figure handle
            O.gui.fighandle = figure(...
                'MenuBar','none',...
                'Name','Looking at CATs',...
                'NumberTitle','off',...
                'Position',[200 200 720 540],...
                'Resize','off');
            
            % Create GUI ui controls
            
            % Panel for thermodynamics settings
            O.gui.td.panel = uipanel('Parent',O.gui.fighandle,...
                'Title','Thermodynamics',...
                'Units','pixels',...
                'Position',[20 310 335 205]);
            
            % Panel for kinetics settings
            O.gui.kin.panel = uipanel('Parent',O.gui.fighandle,...
                'Title','Kinetics',...
                'Units','pixels',...
                'Position',[365 310 335 205]);
            
            % Panel for process settings
            O.gui.proc.panel = uipanel('Parent',O.gui.fighandle,...
                'Title','Process settings',...
                'Units','pixels',...
                'Position',[20 95 335 205]);
            
            % Panel for solver settings
            O.gui.solv.panel = uipanel('Parent',O.gui.fighandle,...
                'Title','Solver settings',...
                'Units','pixels',...
                'Position',[365 95 335 205]);
            
            % Panel for running programs
            O.gui.Run.panel = uipanel('Parent',O.gui.fighandle,...
                'Title','Run',...
                'Units','pixels',...
                'Position',[20 20 680 70]);
            
        end % function
        
        
    end % methods
    
end % classdef