%% Method plot

function out = plot(O,varargin)

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

% Check for calculated data
noresults = 0;
for i = 1:length(O)
    if isempty(O(i).calc_time)
        warning('CAT:plot:noCalcs',...
            'No calculation results available in entry %i',i)
        noresults = noresults + 1;
    end % if
end % if
if noresults == length(O)
    error('CAT:plot:noCalcs',...
        'No calculation results available to plot')
end % if

% Check for existing plots
if isempty(O.handles)
    existing_plots = false(size(available_plots));
else
    existing_plots = ishandle(O.handles);
end % if

% Check inputs
if iscell(varargin{1})
    varargin = [varargin{1} varargin{2:end}];
end % if

% Check if requesting to close all plotted figures - close all figures,
% reset handles structure to empty
if any(strcmpi(varargin,'close'))
    close(O.handles.figures(ishandle(O.handles.figures)));
    O.handles = [];
    return
end % if


% Extract names of plots from varargin - these are plotwhat
plotwhat = intersect(lower(varargin),lower(available_plots));

% Keep rest for varargin
varargin = setdiff(lower(varargin),lower(available_plots));

if nargin == 1 || isempty(plotwhat)
    
    % If no specific plots are mentioned, replot those which are existing -
    % unless none exist, in which case all figures should be plotted
    %
    % This way figures are created the first time plot is called, and any
    % remaining figures are updated the next time plot is called.
    %
    % New plots can be forced by calling C.plot('plotname')
    
    if all(existing_plots == 0)
        requested_plots = true(size(available_plots));
    else
        requested_plots = existing_plots;
    end % if else
    
elseif ischar(plotwhat) && any(strcmpi(plotwhat,{'detailed_results','all'}))
    
    % All the plots are chosen explicitly
    requested_plots = true(size(available_plots));

elseif iscell(plotwhat)
    
    % A specific list of plots was requested - add these to the list of
    % existing plots
    
    requested_plots = existing_plots;
    
    % Go through cell list of what to plot - add each requested one to the
    % list
    for i = 1:length(plotwhat)
        requested_plots = requested_plots | strcmpi(plotwhat{i},available_plots);
    end
    
else
    error('CAT:plot:wronginput','The input given to plot is misformed. See the help for information')
end % if


%%

set(0,'defaultaxesfontsize',14,'defaulttextfontsize',16);


if ishandle(21)
    h = get(21,'children');
    defactoSeries = length(get(h(1),'children'))+1;
else
    defactoSeries = 1;
end

serLen = length(O); % series length
lineProps = {'b-','k--','g-d','m-s','c-o','r+'}; % line properties for series

plothandles = [];
fighandles = [];

if ~isempty(O(1).calc_time)
    % 3D plot of distributions over time
    if serLen==1 && (~isempty(find(strcmp(plotwhat,'distributions'), 1)) ...
            || ~isempty(find(strcmp(plotwhat,'dist3D'), 1)) ...
            || ~isempty(find(strcmp(plotwhat,'results'), 1))...
            || ~isempty(find(strcmp(plotwhat,'detailed_results'), 1))...
            && ~strcmp(O.sol_method,'movingpivot'))
        
        
        for i = 1:length(O.calc_dist)
            Fmat(:,i) = O.calc_dist(i).F;
        end % for
        
        figure(12)
        fighandles = [fighandles 12];
        set(gcf,'numbertitle','off','name','PSDs (3D time evolution)')
        
        % Handles for plots
        PDpl_local = zeros(2,1);
        
        subplot(1,2,1)
        PDpl_local(1) = surf(O.calc_time(:),O.calc_dist(1).y(:),...
            Fmat./repmat(moments(O.calc_dist,0),length(O.calc_dist(1).y),1),varargin{:});
        ylabel('Mean Char. Length')
        xlabel('Time')
        zlabel('Normalized Number Distribution')
        
        subplot(1,2,2)
        PDpl_local(2) = surf(O.calc_time,O.calc_dist(1).y,...
            Fmat.*repmat(O.calc_dist(1).y(:).^3,1,length(O.calc_time))...
            ./repmat(moments(O.calc_dist,3),length(O.calc_dist(1).y),1),...
            varargin{:});
        ylabel('Mean Char. Length')
        xlabel('Time')
        zlabel('Normalized Volume Distribution')
        
        plothandles = [plothandles; PDpl_local];
        
    elseif serLen==1 &&(~isempty(find(strcmp(plotwhat,'distributions'), 1)) ...
            || ~isempty(find(strcmp(plotwhat,'dist3D'), 1)) ...
            || ~isempty(find(strcmp(plotwhat,'results'), 1))...
            || ~isempty(find(strcmp(plotwhat,'detailed_results'), 1))...
            && strcmp(O.sol_method,'movingpivot')...
            && serLen==1)
        
        figure(12)
        fighandles = [fighandles 12];
        set(gcf,'numbertitle','off','name','PSDs (3D time evolution)')
        
        % Handles for plots
        PDpl_local = zeros(2,1);
        
        subplot(1,2,1)
        for i = 1:length(O.calc_time)
            PDpl_local(i) = plot3(repmat(O.calc_time(i),size(O.calc_dist(i).y)),O.calc_dist(i).y(:),...
                O.calc_dist(i).F/moments(O.calc_dist,0,i),varargin{:});
            hold on
        end
        grid on
        hold off
        ylabel('Mean Char. Length')
        xlabel('Time')
        zlabel('Normalized Number Distribution')
        
        subplot(1,2,2)
        for i = 1:length(O.calc_time)
            PDpl_local(i+length(O.calc_time)) = plot3(repmat(O.calc_time(i),size(O.calc_dist(i).y)),O.calc_dist(i).y,...
                O.calc_dist(i).F(:).*O.calc_dist(i).y(:).^3./...
                moments(O.calc_dist,3,i),...
                varargin{:});
            hold on
        end
        grid on
        hold off
        ylabel('Mean Char. Length')
        xlabel('Time')
        zlabel('Normalized Volume Distribution')
        set(PDpl_local,'linewidth',1.5,'color','k')
        plothandles = [plothandles; PDpl_local];
        
    end % if
    
    for ii = 1:serLen
        % Cumulative Properties
        if (~isempty(find(strcmp(plotwhat,'results'), 1)) || ...
                ~isempty(find(strcmp(plotwhat,'detailed_results'), 1)) || ...
                ~isempty(find(strcmp(plotwhat,'cumprop'), 1)))
            
            figure(21)
            fighandles = [fighandles 21];
            set(gcf,'numbertitle','off','name','PSD cumulative properties')
            
            % Handles for plots
            PDpl_local = zeros(3,1);
            
            subplot(3,1,1)
            hold on
            PDpl_local(1) = plot(O(ii).calc_time,moments(O(ii).calc_dist,0),lineProps{ii+defactoSeries},'linewidth',1.5);
            ylabel('0^{th} moment')
            grid on
            hold off
            
            subplot(3,1,2)
            hold on
            PDpl_local(2) = plot(O(ii).calc_time,moments(O(ii).calc_dist,3),lineProps{ii+defactoSeries},'linewidth',1.5);
            ylabel('3^{rd} moment')
            grid on
            hold off
            
            subplot(3,1,3)
            hold on
            PDpl_local(3) = plot(O(ii).calc_time,moments(O(ii).calc_dist,4)./moments(O(ii).calc_dist,3),lineProps{ii+defactoSeries},'linewidth',1.5);
            ylabel('Weight average length')
            xlabel('Time')
            grid on
            hold off
            plothandles = [plothandles; PDpl_local];
        elseif (~isempty(find(strcmp(plotwhat,'detailed_results'), 1)) || ...
                ~isempty(find(strcmp(plotwhat,'moments'), 1)))
            
            figure(22)
            fighandles = [fighandles 22];
            set(gcf,'numbertitle','off','name','Moments Only')
            PDpl_local = zeros(4,1);
            
            subplot(2,2,1)
            hold on
            PDpl_local(1) = plot(O(ii).calc_time,moments(O(ii).calc_dist,0),lineProps{ii+defactoSeries},'linewidth',1.5);
            ylabel('0^{th} moment')
            xlabel('Time')
            grid on
            hold off
            
            subplot(2,2,2)
            hold on
            PDpl_local(1) = plot(O(ii).calc_time,moments(O(ii).calc_dist,1),lineProps{ii+defactoSeries},'linewidth',1.5);
            ylabel('1^{st} moment')
            xlabel('Time')
            grid on
            hold off
            
            subplot(2,2,3)
            hold on
            PDpl_local(1) = plot(O(ii).calc_time,moments(O(ii).calc_dist,2),lineProps{ii+defactoSeries},'linewidth',1.5);
            ylabel('2^{nd} moment')
            xlabel('Time')
            grid on
            hold off
            
            subplot(2,2,4)
            hold on
            PDpl_local(1) = plot(O(ii).calc_time,moments(O(ii).calc_dist,3),lineProps{ii+defactoSeries},'linewidth',1.5);
            ylabel('3^{th} moment')
            xlabel('Time')
            grid on
            hold off
        end % if
        
        % Process Variables
        if (~isempty(find(strcmp(plotwhat,'results'), 1)) || ...
                ~isempty(find(strcmp(plotwhat,'detailed_results'), 1)) ||...
                ~isempty(find(strcmp(plotwhat,'process'), 1)))
            
            figure(31)
            fighandles = [fighandles 31];
            set(gcf,'numbertitle','off','name','Process Variables (I)')
            
            % Handles for plots
            PDpl_local = zeros(1,1);
            
            nopvit = 1;
            if ~isempty(O(ii).calc_conc)
                
                subplot(2,2,nopvit)
                hold on
                PDpl_local = plot(O(ii).calc_time,O(ii).calc_conc,lineProps{ii+defactoSeries},'linewidth',1.5);
                xlim([min(O(ii).calc_time) max(O(ii).calc_time)])
                xlabel('Time')
                ylabel('Concentration')
                grid on
                plothandles = [plothandles; PDpl_local];
                
                nopvit = nopvit + 1;
                hold off
            end % if
            
            if ~isempty(O(ii).calc_conc)
                subplot(2,2,nopvit);
                hold on
                PDpl_local = plot(O(ii).calc_time(:),O(ii).calc_conc(:)./O(ii).solubility(O(ii).Tprofile(O(ii).calc_time(:))),lineProps{ii+defactoSeries},'linewidth',1.5);
                xlabel('Time')
                xlim([min(O(ii).calc_time) max(O(ii).calc_time)])
                ylabel('Supersaturation [-]')
                grid on
                plothandles = [plothandles; PDpl_local(:)];
                hold off
                nopvit = nopvit + 1;
            end % if
            
            subplot(2,2,nopvit)
            hold on
            PDpl_local = plot(O(ii).calc_time,O(ii).Tprofile(O(ii).calc_time),lineProps{ii+defactoSeries},'linewidth',1.5);
            xlabel('Time')
            xlim([min(O(ii).calc_time) max(O(ii).calc_time)])
            ylabel('Temperature')
            grid on
            plothandles = [plothandles; PDpl_local(:)];
            hold off
            nopvit = nopvit + 1;
            
            
            subplot(2,2,nopvit);
            hold on
            PDpl_local = plot(O(ii).calc_time,massmedium(O(ii)),lineProps{ii+defactoSeries},'linewidth',1.5);
            xlabel('Time')
            xlim([min(O(ii).calc_time) max(O(ii).calc_time)])
            ylabel('Total mass Solvent + Antisolvent')
            grid on
            plothandles = [plothandles; PDpl_local(:)];
            
            nopvit = nopvit + 1;
            hold off
            
            if ~isempty(O(ii).calc_conc) && ...
                    (~isempty(find(strcmp(plotwhat,'detailed_results'), 1)) || ...
                    ~isempty(find(strcmp(plotwhat,'process'), 1)))
                PDpl_local = zeros(1,1);
                
                figure(32)
                fighandles = [fighandles 32];
                hold on
                set(gcf,'numbertitle','off','name','Process Variables (II)')
                Tvec = linspace(min(O(ii).Tprofile(O(ii).calc_time))-5,max(O(ii).Tprofile(O(ii).calc_time))+5);
                plot(Tvec,O(ii).solubility(Tvec),'r--','linewidth',1.5)
                legend('Solubility','location','southeast')
                PDpl_local = plot(O(ii).Tprofile(O(ii).calc_time),O(ii).calc_conc,lineProps{ii+defactoSeries},'linewidth',1.5);
                xlabel('Temperature')
                ylabel('Concentration')
                grid on
                plothandles = [plothandles; PDpl_local(:)];
                hold off
            end
            
        end % if
        
        if (~isempty(find(strcmp(plotwhat,'detailed_results'), 1)) || ...
                ~isempty(find(strcmp(plotwhat,'integration'), 1)))
            
            if ~isempty(O(ii).calc_conc)
                figure(41)
                fighandles = [fighandles 41];
                hold on
                set(gcf,'numbertitle','off','name',...
                    'Details from Integration')
                PDpl_local = semilogy(O(ii).calc_time,massbal(O(ii)),lineProps{ii+defactoSeries},'linewidth',1.5);
                xlabel('Time')
                ylabel('Mass balance [% error]')
                grid on
                plothandles = [plothandles; PDpl_local];
                hold off
            end % if
            
            
        end
    end
else
    warning('CAT:plot:noCalcs',...
        'Simulations haven''t run yet, there is nothing to plot');
end
set(fighandles,'units','normalized','position',[0.2 0.3 0.5 0.6]);
set(0,'defaultaxesfontsize','default','defaulttextfontsize','default');

end % function