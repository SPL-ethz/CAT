%% Method plot
        
        function PDpl = plot(O,plotwhat,varargin)
            
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
            % PLOT returns the handles to the plot objects created. (If
            % several plots are created the sequence of handles is the
            % following: PSDs overlapping, PSDs3D, cumulative properties
            % (moments), process variables (temperature, concentration)
            %
            % Graphs can currently not be plotted in existing figures !!
            %
            
            set(0,'defaultaxesfontsize',14,'defaulttextfontsize',16);
            
            if nargin == 1
                plotwhat = 'detailed_results';
            end
            
            serLen = length(O); % series length
            lineProps = {'b-','r-','k-','g-d','m-s','c-o'}; % line properties for series
            
            PDpl = [];
            fhandle = [];

  
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
                fhandle = [fhandle 12];
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

                PDpl = [PDpl; PDpl_local];
                
            elseif serLen==1 &&(~isempty(find(strcmp(plotwhat,'distributions'), 1)) ...
                    || ~isempty(find(strcmp(plotwhat,'dist3D'), 1)) ...
                    || ~isempty(find(strcmp(plotwhat,'results'), 1))...
                    || ~isempty(find(strcmp(plotwhat,'detailed_results'), 1))...
                    && strcmp(O.sol_method,'movingpivot')...
                    && serLen==1)
                
                figure(12)
                fhandle = [fhandle 12];
                set(gcf,'numbertitle','off','name','PSDs (3D time evolution)')
                
                % Handles for plots
                PDpl_local = zeros(2,1);

                subplot(1,2,1)
                for i = 1:length(O.calc_time)
                    PDpl_local(i) = plot3(repmat(O.calc_time(i),size(O.calc_dist(i).y)),O.calc_dist(i).y(:),...
                        O.calc_dist(i).F/moments(O.calc_dist,0),varargin{:});
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
                PDpl = [PDpl; PDpl_local];
                
            end % if
            
            for ii = 1:serLen
                % Cumulative Properties
                if (~isempty(find(strcmp(plotwhat,'results'), 1)) || ...
                    ~isempty(find(strcmp(plotwhat,'detailed_results'), 1)) || ...
                    ~isempty(find(strcmp(plotwhat,'cumprop'), 1)))

                    figure(21)
                    fhandle = [fhandle 21];
                    set(gcf,'numbertitle','off','name','PSD cumulative properties')  

                    % Handles for plots
                    PDpl_local = zeros(3,1);

                    subplot(3,1,1)
                    hold on
                    PDpl_local(1) = plot(O(ii).calc_time,moments(O(ii).calc_dist,0),lineProps{ii});                
                    ylabel('0^{th} moment')
                    grid on
                    hold off
                    
                    subplot(3,1,2)
                    hold on
                    PDpl_local(2) = plot(O(ii).calc_time,moments(O(ii).calc_dist,3),lineProps{ii});
                    ylabel('3^{rd} moment')
                    grid on
                    hold off
                    
                    subplot(3,1,3)
                    hold on
                    PDpl_local(3) = plot(O(ii).calc_time,moments(O(ii).calc_dist,4)./moments(O(ii).calc_dist,3),lineProps{ii});
                    ylabel('Weight average length')
                    xlabel('Time')
                    grid on
                    hold off
                    PDpl = [PDpl; PDpl_local];
                elseif (~isempty(find(strcmp(plotwhat,'detailed_results'), 1)) || ...
                    ~isempty(find(strcmp(plotwhat,'moments'), 1)))

                    figure(22)
                    fhandle = [fhandle 22];
                    set(gcf,'numbertitle','off','name','Moments Only')  
                    PDpl_local = zeros(4,1);

                    subplot(2,2,1)
                    hold on
                    PDpl_local(1) = plot(O(ii).calc_time,moments(O(ii).calc_dist,0),lineProps{ii});
                    ylabel('0^{th} moment')
                    xlabel('Time')
                    grid on
                    hold off
                    
                    subplot(2,2,2)
                    hold on
                    PDpl_local(1) = plot(O(ii).calc_time,moments(O(ii).calc_dist,1),lineProps{ii});
                    ylabel('1^{st} moment')
                    xlabel('Time')
                    grid on
                    hold off
                    
                    subplot(2,2,3)
                    hold on
                    PDpl_local(1) = plot(O(ii).calc_time,moments(O(ii).calc_dist,2),lineProps{ii});
                    ylabel('2^{nd} moment')
                    xlabel('Time')
                    grid on
                    hold off
                    
                    subplot(2,2,4)
                    hold on
                    PDpl_local(1) = plot(O(ii).calc_time,moments(O(ii).calc_dist,3),lineProps{ii});
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
                    fhandle = [fhandle 31];
                    set(gcf,'numbertitle','off','name','Process Variables (I)')

                % Handles for plots
                    PDpl_local = zeros(1,1);

                    nopvit = 1;
                    if ~isempty(O(ii).calc_conc)

                        subplot(2,2,nopvit)
                        hold on
                        PDpl_local = plot(O(ii).calc_time,O(ii).calc_conc,lineProps{ii},'linewidth',1.5);
                        xlim([min(O(ii).calc_time) max(O(ii).calc_time)])
                        xlabel('Time')
                        ylabel('Concentration')
                        grid on
                        PDpl = [PDpl; PDpl_local];

                        nopvit = nopvit + 1;
                        hold off
                    end % if

                    if ~isempty(O(ii).calc_conc)
                        subplot(2,2,nopvit);              
                        hold on
                        PDpl_local = plot(O(ii).calc_time(:),O(ii).calc_conc(:)./O(ii).solubility(O(ii).Tprofile(O(ii).calc_time(:))),lineProps{ii},'linewidth',1.5);
                        xlabel('Time')
                        xlim([min(O(ii).calc_time) max(O(ii).calc_time)])
                        ylabel('Supersaturation [-]')
                        grid on
                        PDpl = [PDpl; PDpl_local(:)];
                        hold off
                        nopvit = nopvit + 1;
                    end % if

                    subplot(2,2,nopvit)
                    hold on
                    PDpl_local = plot(O(ii).calc_time,O(ii).Tprofile(O(ii).calc_time),lineProps{ii},'linewidth',1.5);
                    xlabel('Time')
                    xlim([min(O(ii).calc_time) max(O(ii).calc_time)])
                    ylabel('Temperature')
                    grid on
                    PDpl = [PDpl; PDpl_local(:)];
                    hold off
                    nopvit = nopvit + 1;


                    subplot(2,2,nopvit);               
                    hold on
                    PDpl_local = plot(O(ii).calc_time,massmedium(O(ii)),lineProps{ii},'linewidth',1.5);
                    xlabel('Time')
                    xlim([min(O(ii).calc_time) max(O(ii).calc_time)])
                    ylabel('Total mass Solvent + Antisolvent')
                    grid on
                    PDpl = [PDpl; PDpl_local(:)];

                    nopvit = nopvit + 1;
                    hold off

                    if ~isempty(O(ii).calc_conc) && ...
                            (~isempty(find(strcmp(plotwhat,'detailed_results'), 1)) || ...
                            ~isempty(find(strcmp(plotwhat,'process'), 1)))
                        PDpl_local = zeros(1,1);

                        figure(32)
                        fhandle = [fhandle 32];
                        hold on
                        set(gcf,'numbertitle','off','name','Process Variables (II)')
                        Tvec = linspace(min(O(ii).Tprofile(O(ii).calc_time))-5,max(O(ii).Tprofile(O(ii).calc_time))+5);
                        plot(Tvec,O(ii).solubility(Tvec),'r--','linewidth',1.5)
                        legend('Solubility','location','southeast')
                        PDpl_local = plot(O(ii).Tprofile(O(ii).calc_time),O(ii).calc_conc,lineProps{ii},'linewidth',1.5);
                        xlabel('Temperature')
                        ylabel('Concentration')
                        grid on
                        PDpl = [PDpl; PDpl_local(:)];
                        hold off
                    end

                end % if

                if (~isempty(find(strcmp(plotwhat,'detailed_results'), 1)) || ...
                    ~isempty(find(strcmp(plotwhat,'integration'), 1)))

                    if ~isempty(O(ii).calc_conc)
                        figure(41)
                        fhandle = [fhandle 41];
                        hold on
                        set(gcf,'numbertitle','off','name',...
                            'Details from Integration')
                        PDpl_local = semilogy(O(ii).calc_time,massbal(O(ii)));
                        xlabel('Time')
                        ylabel('Mass balance [% error]')
                        grid on
                        PDpl = [PDpl; PDpl_local];
                        hold off
                    end % if


                end
            end
         
            set(fhandle,'units','normalized','position',[0.2 0.3 0.5 0.6]);
            set(0,'defaultaxesfontsize','default','defaulttextfontsize','default');
  
            
        end % function