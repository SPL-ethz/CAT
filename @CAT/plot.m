%% Method plot

function out = plot(Oin,varargin)

% Plot function for the results
%
% Use plot(PD,plotwhat) to plot the results of the
% simulation. plotwhat is a string that defines what exactly
% should be plotted. Possible input:
% 'results'         -   plot everything
% 'detailed_results'-   plot everything and more
%
% plotted using results
% 'distributions'   -   plot distributions
% 'distoverlap'     -   only plot overlapping distributions (2D)
% 'dist3D'          -   only plot 3D surf plot of distributions
% 'cumprop'         -   plot cumulative properties (moments)
% 'process'         -   plot process variables (T, conc)
%
% additionally plotted in detailed mode
% 'moments'         -   plots of the first four moments
% 'integration'     -   details from the integration
% (massbalance over time,...)
%
% and any combination thereof.
%
% The syntax
%   pl = plot(PD,plotwhat)
% returns the handles to the plot objects created. (If
% several plots are created the sequence of handles is the
% following: PSDs overlapping, PSDs3D, cumulative properties
% (moments), process variables (temperature, concentration)
%
% To close all the plots created by this plotting function,
% use:
%   plot(PD,'close')
%
% Graphs can currently not be plotted in existing figures !!
%

%% Define list of available plots

available_plots = {...
    'PSD_3D',... Three dimensional plots of particle size distributions vs time
    'PSD_cumul',... Cumulative properties (moments) of PSDs vs time
    'proc_var',... Process variables vs time: concentration, supersaturation, temperature, total solvent mass
    'operating',... Operating diagram: concentration vs temperature
    'massbal'... Mass balance
    };

%% Figure out what to plot

% Check if requesting list of available figures - this returns cell
% available_plots and exits
if any(strcmpi(varargin,'available_plots'))
    out = available_plots;
    return
end % if

% Check if requesting to close all plotted figures - close all figures,
% reset handles structure to empty
if any(strcmpi(varargin,'close'))
    all_handles = [Oin.handles_figures];
    close( all_handles( ishandle(all_handles) ) );
    [O.handles_figures] = deal([]);
    return
end % if

% Check for calculated data
noresults = 0;
for i1 = 1:length(Oin)
    if isempty(Oin(i1).calc_time)
        warning('CAT:plot:noCalcs',...
            'No calculation results available in entry %i',i1)
        noresults = noresults + 1;
    end % if
end % for
if noresults == length(Oin)
    error('CAT:plot:noCalcs',...
        'No calculation results available to plot')
end % if

% Check for existing plots for each entry, update existing_plots vector
%
% Wanted behaviour:
%   - no figures present: create new figures (all of them), plot all data
%   into same figrues
%   - figures present for one series but not for others: use these figures,
%   plot all data into them
%   - figures available for two or more series: use these, create new
%   figures where they do not exist but are requested

% Check inputs
if ~isempty(varargin) && iscell(varargin{1})
    varargin = [varargin{1} varargin{2:end}];
end % if

% Extract names of plots from varargin - these are plotwhat, the list of
% explicitly requested plots
plotwhat = intersect(lower(varargin),lower(available_plots));

% Keep rest for varargin - are these really used?
varargin = setdiff(lower(varargin),lower(available_plots));

% Find which figures already exist
% Set empty handles_figures vectors to NaN
for O = Oin
    if isempty(O.handles_figures)
        O.handles_figures = NaN(1,length(available_plots));
    end % if
end % for
existing_plots = ishandle( ...
    reshape([Oin(:).handles_figures]',length(available_plots),length(Oin))' ...
    );

% Count how many series have existing plots
existing_plots_series = any(existing_plots,2);

switch sum(existing_plots_series)
    % Do different things depending on how many series have available plots
    
    case 0
        % No plots exist - All series should be plotted into same figure
        %
        % Check if any figures are requested explicitly,
        % or if all figures should be created
        
        if nargin == 1 || isempty(plotwhat) || ...
                ischar(plotwhat) && any(strcmpi(plotwhat,{'detailed_results','all'}))
            requested_plots = true(1,length(available_plots));
        elseif exist('plotwhat','var') && iscell(plotwhat)
            % Go through cell list of what to plot - add each requested one to the
            % list
            requested_plots = false(1,length(available_plots));
            for j = 1:length(plotwhat)
                requested_plots = requested_plots | strcmpi(plotwhat{j},available_plots);
            end
        end % if
        
    case 1
        % Plots exist for one series - all series should be plotted into
        % the existing figures.
        %
        % Update the existing figures, add explicitly requested plots
        
        % Already available figures
        requested_plots = existing_plots(existing_plots_series == 1,:);
        
        % Apply the handles of the existing series to all the input series
        [Oin.handles_figures] = deal(Oin(existing_plots_series == 1).handles_figures);
        
        if exist('plotwhat','var') && iscell(plotwhat)
            
            % A specific list of plots was requested - add these to the list of
            % existing plots
            
            % Go through cell list of what to plot - add each requested one to the
            % list
            for j = 1:length(plotwhat)
                requested_plots = requested_plots | strcmpi(plotwhat{j},available_plots);
            end
            
        end % if
        
    otherwise
        % Plots exist for 2 or more series - each series should be plotted
        % into its own figure.
        %
        % In this case, plots are created/updated if they exist for any
        % series, or if they are explicitly requested
        
        % Already available figures
        requested_plots = any(existing_plots,1);
        
        if any( sum(existing_plots,1) == 1 )
            % For at least one plot, the plot only exists for one series -
            % plot all series into the same figure
            
            singleplots = find( sum(existing_plots,1) == 1);
            
            for S = singleplots
                % Go through the list of plots which are available only for
                % one series and set the handle for the existing plot to
                % the ones where this plot does not exist yet
                for O = Oin(~existing_plots(:,S))
                    O.handles_figures(S) = Oin(existing_plots(:,S)).handles_figures(S);
                end % for
            end % for
            
            keyboard
            
        end % if
        
        if exist('plotwhat','var') && iscell(plotwhat)
            
            % A specific list of plots was requested - add these to the list of
            % existing plots
            
            % Go through cell list of what to plot - add each requested one to the
            % list
            for j = 1:length(plotwhat)
                requested_plots = requested_plots | strcmpi(plotwhat{j},repmat(available_plots,length(Oin),1));
            end
            
        end % if
        
end % switch

% The existing_plots variable is no longer needed, as it will be checked
% agains the actual figure handles later
clear existing_plots

% requested_plots now has one row
%
% The handles_figures vector in each input series has the figure handle of
% the figure to use to plot its values, if the figure exists
%
% Now, we need to make or update all the plots

%% Make plots

% Loop through the plots

for p = find(requested_plots)
    
    % Call the subfunction which does the work
    makeplot(p)
    
end % for

%% Plotting subfunction

    function makeplot(p)
        
        % Makeplot makes one plot, according to which one is requested, for
        % every series.
        %
        % For each series, the plot is put into the given figure handle, if
        % it exists, otherwise, a new one is created
        
        % Increase default text sizes
        set(0,'defaultaxesfontsize',14,'defaulttextfontsize',16);
        
        % Identify which plot to make by the list of available plots and
        % the number
        whichplot = available_plots{p};
        
        % Loop on O
        for O = Oin %#ok<FXUP>
            
            % Check for existance of figure - or create it
            if ~ishandle(O.handles_figures(p))
                % The figure needs to be created
                O.handles_figures(p) = figure;
            end % if
            
            % Everything else is figure-specific, so go do specific stuff
            % now
            
            switch whichplot
                
                case 'PSD_3D'
                    
                    %% Plot: PSD_3D
                    % Three dimensional plots of particle size distributions vs time
                    
                    % Extract handles for simplicity
                    figh = O.handles_figures(p);
                    if ~isempty(O.handles_axes) && length(O.handles_axes) >= p
                        axh = O.handles_axes{p};
                    else
                        % Use NaN here - subplot checks for existance of
                        % any axes when it is called
                        axh = [NaN NaN];
                    end % if else
                    
                    % Setup the figure title
                    set(figh,'numbertitle','off','name','PSDs (3D time evolution)');
                    
                    % Check for axes - create and change them only if they
                    % don't exist
                    if ~ishandle(axh(1))
                        axh(1) = subplot(1,2,1,'Parent',figh);
                        
                        % Grid
                        grid(axh(1),'on')
                        
                        % Axes labels
                        xlabel(axh(1),'Time')
                        ylabel(axh(1),'Mean Char. Length')
                        zlabel(axh(1),'Normalized Number Distribution')
                        
                        % Default view
                        view(axh(1),[-32 30])
                        
                    end % if
                    if ~ishandle(axh(2))
                        axh(2) = subplot(1,2,2,'Parent',figh);
                        
                        % Grid
                        grid(axh(2),'on')
                        
                        % Axes labels
                        xlabel(axh(2),'Time')
                        ylabel(axh(2),'Mean Char. Length')
                        zlabel(axh(2),'Normalized Volume Distribution')
                        
                        % Default view
                        view(axh(2),[-32 30])
                        
                    end % if
                    % Save these handles
                    O.handles_axes{p} = axh;
                    
                    % Calculate some stuff necessary for plotting - how
                    % this is done depends on whether the size vectors are
                    % equal or not.
                    %
                    % In for example moving pivot, the size vector changes
                    if isequal(O.calc_dist.y)
                        
                        % Preallocate - put times in rows, size in columns
                        [Fmat,FmatV] = deal( zeros(length(O.calc_time),length(O.calc_dist(1).y)) );
                        
                        for i2 = 1:length(O.calc_time)
                            % Normalize only the first distribution
                            Fmat(i2,:) = O.calc_dist(i2).F/moments(O.calc_dist(1),0) ;
                            FmatV(i2,:) = O.calc_dist(i2).F.*O.calc_dist(i2).y.^3 / moments(O.calc_dist(1),3) ;
                        end % for
                        
                        % Assign times and sizes for simplicity
                        time = O.calc_time;
                        y = O.calc_dist(1).y;
                        
                        % Make the plots - use surface, not surf, as this does
                        % not affect the existing axes
                        %
                        % If plots are to be replaced, the old objects should
                        % be deleted before plotting the new ones
                        O.handles_objects{p} = {...
                            surface(...
                            time,...
                            y,...
                            Fmat',...
                            'Parent',axh(1)),...
                            surface(...
                            time,...
                            y,...
                            FmatV',...
                            'Parent',axh(2))...
                            };
                        
                    else
                        % The size vectors are not equal - surf will not
                        % work in the same way (size vectors can have
                        % different lengths)
                        
                        % Plot each distribution as a line in 3D space
                        % Connect all the distributions by NaN values so
                        % that they are not visually connected but are in
                        % the same object
                        
                        time = [];
                        y = [];
                        Fmat = [];
                        FmatV = [];
                        
                        for i3 = 1:length(O.calc_time)
                            
                            currsize = O.calc_dist(i3).y;
                            currdist = O.calc_dist(i3).F/moments(O.calc_dist(1),0);
                            currdistV = O.calc_dist(i3).F.*O.calc_dist(i3).y.^3 / moments(O.calc_dist(1),3) ;
                            
                            y = [y; currsize(:); NaN]; %#ok<AGROW>
                            time = [time; repmat(O.calc_time(i3),length(currsize),1); NaN]; %#ok<AGROW>
                            Fmat = [Fmat; currdist(:); NaN]; %#ok<AGROW>
                            FmatV = [FmatV; currdistV(:); NaN]; %#ok<AGROW>
                            
                        end % for
                        
                        % Make the plots
                        % Use line instead of plot3 because this does not
                        % affect the axes (which have been changed outside
                        % of CAT manually).
                        %
                        % If plots are to be replaced, the old objects should
                        % be deleted before plotting the new ones
                        O.handles_objects{p} = {...
                            line(...
                            time,...
                            y,...
                            Fmat,...
                            'Parent',axh(1)),...
                            line(...
                            time,...
                            y,...
                            FmatV,...
                            'Parent',axh(2))...
                            };
                        
                    end % if else
                    
                case 'PSD_cumul'
                    
                    %% Plot: PSD_cumul
                    % Cumulative properties (moments) of PSDs vs time
                    
                    % Extract handles for simplicity
                    figh = O.handles_figures(p);
                    if ~isempty(O.handles_axes) && length(O.handles_axes) >= p
                        axh = O.handles_axes{p};
                    else
                        % Use NaN here - subplot checks for existance of
                        % any axes when it is called
                        axh = [NaN NaN NaN];
                    end % if else
                    
                    % Setup the figure title
                    set(figh,'numbertitle','off','name','PSD cumulative properties');
                    
                    % Check for axes - create and change them only if they
                    % don't exist
                    if ~ishandle(axh(1))
                        axh(1) = subplot(3,1,1,'Parent',figh);
                        
                        % Grid, box
                        grid(axh(1),'on')
                        box(axh(1),'on')
                        
                        % Axes labels
                        ylabel(axh(1),'0^{th} moment')
                        
                    end % if
                    if ~ishandle(axh(2))
                        axh(2) = subplot(3,1,2,'Parent',figh);
                        
                        % Grid, box
                        grid(axh(2),'on')
                        box(axh(2),'on')
                        
                        % Axes labels
                        ylabel(axh(2),'3^{rd} moment')
                        
                    end % if
                    if ~ishandle(axh(3))
                        axh(3) = subplot(3,1,3,'Parent',figh);
                        
                        % Grid, box
                        grid(axh(3),'on')
                        box(axh(3),'on')
                        
                        % Axes labels
                        xlabel(axh(3),'Time')
                        ylabel(axh(3),'Weight average length')
                        
                    end % if
                    % Save these handles
                    O.handles_axes{p} = axh;
                    
                    % Link the x axes so that they behave the same on
                    % zoom/pan
                    linkaxes(axh,'x')
                    
                    % Make the plots
                    %
                    % Delete current objects first?
                    O.handles_objects{p} = {...
                        line(O.calc_time,moments(O.calc_dist,0),'Parent',axh(1)),...
                        line(O.calc_time,moments(O.calc_dist,3),'Parent',axh(2)),...
                        line(O.calc_time,moments(O.calc_dist,4)./moments(O.calc_dist,3),'Parent',axh(3))...
                        };
                    
                    % Update the limits - use maximum of current limit or
                    % x/y values
                    
                    % Since axes are linked, only do this for the first x
                    % axes, but for all y axes
                    xl = xlim(axh(1));
                    xlim(axh(1), [ ...
                        min([xl(1) min(O.calc_time)]) ...
                        max([xl(2) max(O.calc_time)]) ...
                        ] );
                    
                    yl = ylim(axh(1));
                    ylim(axh(1), [ ...
                        min([yl(1) min(moments(O.calc_dist,0))]) ...
                        max([yl(2) max(moments(O.calc_dist,0))]) ...
                        ] );
                    
                    yl = ylim(axh(2));
                    ylim(axh(2), [ ...
                        min([yl(1) min(moments(O.calc_dist,3))]) ...
                        max([yl(2) max(moments(O.calc_dist,3))]) ...
                        ] );
                    
                    yl = ylim(axh(3));
                    ylim(axh(3), [ ...
                        min([yl(1) min(moments(O.calc_dist,4)./moments(O.calc_dist,3))]) ...
                        max([yl(2) max(moments(O.calc_dist,4)./moments(O.calc_dist,3))]) ...
                        ] );
                    
                case 'proc_var'
                    
                    %% Plot: proc_var
                    % Process variables vs time: concentration, supersaturation, temperature, total solvent mass
                    
                    % Extract handles for simplicity
                    figh = O.handles_figures(p);
                    if ~isempty(O.handles_axes) && length(O.handles_axes) >= p
                        axh = O.handles_axes{p};
                    else
                        % Use NaN here - subplot checks for existance of
                        % any axes when it is called
                        axh = [NaN NaN NaN NaN];
                    end % if else
                    
                    % Setup the figure title
                    set(figh,'numbertitle','off','name','Process Variables');
                    
                    % Check for axes - create and change them only if they
                    % don't exist
                    if ~ishandle(axh(1))
                        axh(1) = subplot(2,2,1,'Parent',figh);
                        
                        % Grid, box
                        grid(axh(1),'on')
                        box(axh(1),'on')
                        
                        % Axes labels
                        xlabel(axh(1),'Time')
                        ylabel(axh(1),'Concentration')
                        
                    end % if
                    if ~ishandle(axh(2))
                        axh(2) = subplot(2,2,2,'Parent',figh);
                        
                        % Grid, box
                        grid(axh(2),'on')
                        box(axh(2),'on')
                        
                        % Axes labels
                        xlabel(axh(2),'Time')
                        ylabel(axh(2),'Supersaturation [-]')
                        
                    end % if
                    if ~ishandle(axh(3))
                        axh(3) = subplot(2,2,3,'Parent',figh);
                        
                        % Grid, box
                        grid(axh(3),'on')
                        box(axh(3),'on')
                        
                        % Axes labels
                        xlabel(axh(3),'Time')
                        ylabel(axh(3),'Temperature')                        
                        
                    end % if
                    if ~ishandle(axh(4))
                        axh(4) = subplot(2,2,4,'Parent',figh);
                        
                        % Grid, box
                        grid(axh(4),'on')
                        box(axh(4),'on')
                        
                        % Axes labels
                        xlabel(axh(4),'Time')
                        ylabel(axh(4),'Total mass Solvent + Antisolvent')
                        
                    end % if
                    
                    % Save these handles
                    O.handles_axes{p} = axh;
                    
                    % Link the x axes so that they behave the same on
                    % zoom/pan
                    linkaxes(axh,'x')
                    
                    % Calculate some extra data
                    Tprof = evalanonfunc(O.Tprofile,O.calc_time);
                    
                    if length(Tprof) == 1
                        Tprof = Tprof*ones(size(O.calc_time));
                    end % if
                    
                    % Make the plots
                    %
                    % Delete current objects first?
                    O.handles_objects{p} = {...
                        line(O.calc_time,O.calc_conc,'Parent',axh(1)),...
                        line(O.calc_time,O.calc_conc./evalanonfunc(O.solubility,evalanonfunc(O.Tprofile,O.calc_time),O.massmedium),'Parent',axh(2)),...
                        line(O.calc_time,Tprof,'Parent',axh(3))...
                        line(O.calc_time,O.massmedium,'Parent',axh(4))...
                        };
                    
                    % Since axes are linked, only do this for the first x
                    % axes, but for all y axes
                    xl = xlim(axh(1));
                    xlim(axh(1), [ ...
                        min([xl(1) min(O.calc_time)]) ...
                        max([xl(2) max(O.calc_time)]) ...
                        ] );
                    
                    yl = ylim(axh(1));
                    ylim(axh(1), [ ...
                        min([yl(1) min(O.calc_conc)]) ...
                        max([yl(2) max(O.calc_conc)]) ...
                        ] );
                    
                    yl = ylim(axh(2));
                    ylim(axh(2), [ ...
                        min([yl(1) min(O.calc_conc./evalanonfunc(O.solubility,evalanonfunc(O.Tprofile,O.calc_time),O.massmedium))]) ...
                        max([yl(2) max(O.calc_conc./evalanonfunc(O.solubility,evalanonfunc(O.Tprofile,O.calc_time),O.massmedium))]) ...
                        ] );
                    
                    yl = ylim(axh(3));
                    ylim(axh(3), [ ...
                        min([yl(1) min(evalanonfunc(O.Tprofile,O.calc_time))]) ...
                        max([yl(2) max(evalanonfunc(O.Tprofile,O.calc_time))]) ...
                        ] );
                    
                    yl = ylim(axh(4));
                    ylim(axh(4), [ ...
                        min([yl(1) min(O.massmedium)]) ...
                        max([yl(2) max(O.massmedium)]) ...
                        ] );
                    
                case 'operating'
                    
                    %% Plot: operating
                    % Operating diagram: concentration vs temperature
                    
                    % Extract handles for simplicity
                    figh = O.handles_figures(p);
                    if ~isempty(O.handles_axes) && length(O.handles_axes) >= p
                        axh = O.handles_axes{p};
                    else
                        % Use NaN here - subplot checks for existance of
                        % any axes when it is called
                        axh = NaN;
                    end % if else
                    
                    % Setup the figure title
                    set(figh,'numbertitle','off','name','Operating diagram');
                    
                    % Check for axes - create and change them only if they
                    % don't exist
                    if ~ishandle(axh(1))
                        axh(1) = axes('Parent',figh);
                        
                        % Grid
                        grid(axh(1),'on')
                        box(axh(1),'on')
                        
                        % Axes labels
                        xlabel(axh(1),'Temperature')
                        ylabel(axh(1),'Concentration')
                        
                        % Setup the legend
                        legend(axh(1),'hide','Location','NorthWest')
                        
                    end % if
                    
                    % Save these handles
                    O.handles_axes{p} = axh;
                    
                    % Calculate some extra stuff
                    Tprof = evalanonfunc(O.Tprofile,O.calc_time);
                    Tvec = linspace(min(Tprof)-5,max(Tprof)+5);
                    sol = evalanonfunc(O.solubility,Tvec);
                    
                    % Make sure Tprof is a vector for a proper plot
                    if length(Tprof) == 1
                        Tprof = Tprof*ones(size(O.calc_time));
                    end % if
                    % Make sure sol is a vector for a proper plot
                    if length(sol) == 1
                        sol = sol*ones(size(Tvec));
                    end % if
                    
                    % Make the plots
                    %
                    % Delete current objects first?
                    O.handles_objects{p} = {...
                        line(Tvec,sol,'Displayname','Solubility','Parent',axh(1)),...
                        line(Tprof,O.calc_conc,'Parent',axh(1))...
                        };
                    
                    % Turn the legend on
                    legend(O.handles_objects{p}{1})
                    
                case 'massbal'
                    
                    %% Plot: massbal
                    % Mass balance of integration
                    
                    % Extract handles for simplicity
                    figh = O.handles_figures(p);
                    if ~isempty(O.handles_axes) && length(O.handles_axes) >= p
                        axh = O.handles_axes{p};
                    else
                        % Use NaN here - subplot checks for existance of
                        % any axes when it is called
                        axh = NaN;
                    end % if else
                    
                    % Setup the figure title
                    set(figh,'numbertitle','off','name','Details from Integration');
                    
                    % Check for axes - create and change them only if they
                    % don't exist
                    if ~ishandle(axh(1))
                        axh(1) = axes('Parent',figh);
                        
                        % Grid
                        grid(axh(1),'on')
                        box(axh(1),'on')
                        
                        % Axes labels
                        xlabel(axh(1),'Time')
                        ylabel(axh(1),'Mass balance [% error]')
                        
                        % Set Yscale to log
                        set(axh(1),'YScale','log');
                        
                    end % if
                    
                    % Save these handles
                    O.handles_axes{p} = axh;
                    
                    % Make the plots
                    %
                    % Delete current objects first?
                    O.handles_objects{p} = {...
                        line(O.calc_time,O.massbal,'Parent',axh(1))...
                        };
                    
                    
                otherwise
                    % This point should never be reached, unless
                    % available_plots is changed at the beginning of this
                    % file without adding a block for it here
                    error('CAT:plot:unknownfigure','Something has gone wrong. A non-existant figure was requested');
                    
            end % switch
            
        end % for
        
        % Reset the default text sizes
        set(0,'defaultaxesfontsize','default','defaulttextfontsize','default');
        
    end % function

end % function plot

% lineProps = {'b-','k--','g-d','m-s','c-o','r+'}; % line properties for series
% 
