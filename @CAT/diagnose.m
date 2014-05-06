%% Mehod diagnose
        % This function diagnoses whether the set values for individual
        % fields are valid and whether the CAT instance itself is runnable.
        function [validInput,solvable] = diagnose(O,fieldnames,values,varargin)
            
            if nargin>1 && ~isempty(fieldnames) && ~iscell(fieldnames)
                fieldnames = {fieldnames};
            end
            
            if nargin == 1
                fieldnames = properties(O);
                
                values = NaN;
            end
            Iquery = length(fieldnames);
            
                       
            if  ~isempty(fieldnames) && ~iscell(values)
                values = {values};
            end
            
            if nargin>3 && any(strcmp(varargin,'solvable'))
                fieldnames = [fieldnames {'ASprofile','Tprofile','init_dist','init_conc','solubility','init_seed','init_massmedium','growthrate','rhoc','kv','sol_time'}];    
                values = [values num2cell(ones(size(fieldnames))*NaN)];
            end
            
            validInput = ones(length(fieldnames),1);
            for i = 1:length(fieldnames)
                
                fieldname = fieldnames{i};
                if (~isempty(values{i}) && ~isa(values{i},'function_handle') && (isscalar(values{i}) && isnan(values{i}))) || (~isempty(values) && length(values)<i)
                    value = O.(fieldnames{i});
                else
                    value = values{i};
                end
                
                % analyzing individual fields
                if strcmp(fieldname,'rhoc')
                    if ~(isnumeric(value) && length(value)==1) || value<=0
                        
                        if i<=Iquery
                            warning('CAT:Setrhoc:WrongType',...
                            'The rhoc property must be a scalar.');
                        end

                    validInput(i) = 0;
                    end
                elseif strcmp(fieldname,'kv')
                    
                    if ~(isnumeric(value) && length(value)==1) || value<=0 || value>1
                        if i<=Iquery
                            warning('CAT:Setkv:WrongType',...
                            'The kv property must be a scalar in the range (0,1].');
                        end
                    validInput(i) = 0;
                    end
                    
                elseif strcmp(fieldname,'init_dist')
                    
                    if ~isa(value,'Distribution')
                        if i<=Iquery
                            warning('CAT:SetInit_Dist:WrongType',...
                                'The init_dist property must be a Distribution object');
                        end

                    validInput(i) = 0;
                    end
                
                elseif strcmp(fieldname,'growthrate')
                    if isempty(value) || (~isnumeric(value) && ~ischar(value) && ~isa(value,'function_handle'))
                        if i<=Iquery
                            warning('Distribution:setgrowthrate:Wrongtype',...
                            'The growth rate must be a non-negative, finite value, or a function handle with up to three inputs');
                        end
                        validInput(i) = 0;
                    end
                
                elseif strcmp(fieldname,'nucleationrate')
                    
                    if ~isnumeric(value) && ~ischar(value) && ~isa(value,'function_handle') && ~isempty(value)
                        if i<=Iquery
                            warning('Distribution:setnucleationrate:Wrongtype',...
                                'The nucleation rate must be a non-negative, finite value, empty, or a function handle with up to three inputs');
                        end
                        validInput(i) = 0;
                    end
                
                elseif strcmp(fieldname,'sol_time')
                    if (isnumeric(value) && length(value) > 1  && any(diff(value)<=0)) || ~isnumeric(value) || any(value<0) || isempty(value)
                        
                        if i<=Iquery
                        warning('CAT:SetSol_Time:WrongType',...
                            'The sol_time property must be a vector of monotonically increasing values or a positive scalar');
                        end
                        validInput(i) = 0;
                        
                    end
                    
                elseif strcmp(fieldname,'solubility')
                    if ((isnumeric(value) && ~isscalar(value)) && ~isa(value,'function_handle') && ~ischar(value)) || (isnumeric(value) && value<=0) || isempty(value)
                        
                        if i<=Iquery
                        warning('CAT:SetSolubility:WrongType',...
                            'The solubility property must be a positive, finite value or a function handle with one or two inputs');
                        end
                        validInput(i) = 0;
                    end
                    
                elseif strcmp(fieldname,'init_conc')
                    if isempty(value) || (ischar(value) && ~strcmp(value,'sat')) || (isscalar(value) && value <= 0 && ~isfinite(value))
                        if i<=Iquery
                        warning('CAT:SetInit_Conc:WrongType',...
                            'The init_conc property must be a positive, finite scalar (may be zero) or the string ''sat''');
                        end
                        validInput(i) = 0;
                    end
                    
                elseif strcmp(fieldname,'init_seed')
                    if isempty(value) || (isnumeric(value) && length(value)>1) || isnan(value) || value<0 || isinf(value)
                        if i<=Iquery
                        warning('CAT:SetInit_Seed:WrongType',...
                            'The init_seed property must be a non-negative,finite scalar');
                        end
                        validInput(i) = 0;
                    end
                    
                elseif strcmp(fieldname,'init_massmedium')
                    
                    if isempty(value) || (isnumeric(value) && length(value)>1) || isnan(value) || value<=0 || isinf(value)
                        if i<=Iquery
                        warning('CAT:SetInit_Massmedium:WrongType',...
                            'The init_conc property must be a positive, finite scalar');
                        end
                        validInput(i) = 0;
                
                    end
                    
                elseif strcmp(fieldname,'sol_options')
                    
                    if isempty(value)
                        value = {''};
                    end
                    
                    if ~iscell(value)
                        if i<=Iquery
                        warning('CAT:SetSol_Options:WrongType',...
                            'The sol_options property must be a cell object');
                        end
                        
                        validInput(i) = 0;
                    end
                    
                elseif strcmp(fieldname,'Tprofile')
                    
                    if  ~isa(value,'function_handle') && ~(isnumeric(value) && isscalar(value)) && ~(~isempty(value) && isnumeric(value) && ismatrix(value) && length(value(:,1)) == 2)
                        if i<=Iquery
                        warning('CAT:SetTprofile:WrongType',...
                            'The Tprofile property must be a positive, finite matrix (may be empty) or a function handle with one input');
                        end
                        validInput(i) = 0;
                    end
                    
                elseif strcmp(fieldname,'ASprofile')
                    
                    
                    if ~isa(value,'function_handle') && ~(isnumeric(value) && isscalar(value)) && ~(~isempty(value) && isnumeric(value) && ismatrix(value) && length(value(:,1)) == 2)
                        if i<=Iquery
                        warning('CAT:SetASprofile:WrongType',...
                            'The ASprofile property must be a positive, finite matrix (may be empty) or a function handle with one input');
                        end
                        validInput(i) = 0;
                    end
                    
                elseif strcmp(fieldname,'sol_method')
                    
                    if ~ischar(value) && ~isempty(value)
                        if i<=Iquery
                        warning('CAT:SetSol_Method:WrongType',...
                            'The sol_method property must be a string or empty (chooses default)');
                        end
                        validInput(i) = 0;
                    end
                    
                elseif ~isempty(strfind(fieldname,'calc_'))
                    ; % calculated fields may be empty
                    
                else
                    validInput(i) = NaN;
                end
            end
            
            if nargin>3 && any(strcmp(varargin,'solvable'))
                if any(validInput(1:2)) && all(validInput(3:end))
                    solvable = 1;
                else
                    solvable = 0;
                end
                
            else
                solvable = [];
            end
            
            if Iquery > 0
                validInput = validInput(1:Iquery);
            end
                
        end
        