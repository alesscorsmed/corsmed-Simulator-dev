function [anatomicalModel,expControl] = updateTissueProperties(expControl,...
    anatomicalModel, gamma, pulseSeqFamilyName)
%
% COILS.UPDATETISSUEPROPERTIES
%
%     Updates the selected coil, positioning and SAR values.
%
% INPUT
%   expControl      control struct with experiment info
%   anatomicalModel struct with model data
%
% OUTPUT
%   anatomicalModel updated struct with model data
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'anatomical.updateTissueProperties';
if (nargin < 2)
    ME = MException('Domain:wrongArgCount',...
        '%s : wrong argument count',functionName);
    throw(ME);
elseif (nargin < 3)
    gamma = 42.577478518e6;
end

%% info for debugging
if expControl.debug.debugMode
    try % open file if possible, otherwise dump to stdout
        fid = fopen(expControl.debug.debugFile,'a');
    catch
        fid = 1;
    end
    % time it
    tTotal = tic();
    fprintf(fid, '\n%s : start', functionName);
end

%% Run checks
% Deactivate CS if this is SE-based sequence
if expControl.simulation.activateCS == 1 && ...
        contains(lower(pulseSeqFamilyName), 'se')
    expControl.simulation.activateCS = 0;
    msg = 'This is a SE-based sequence. CS will be deactivated.';
    fprintf(1, '\nWARNING: %s\n',msg);
    if expControl.advancedNotifications
        eduTool.frontend.updateExperimentProgress(expControl,'','info',msg)
    end
end

if (expControl.anatomicalID ~= 6) && ...
        expControl.simulation.activateSusc && ...
        expControl.advancedNotifications
    msg = 'Susceptibility is only available for the brain model.';
    fprintf(1, '\nWARNING: %s\n',msg);
    eduTool.frontend.updateExperimentProgress(expControl,'','info',msg)
    expControl.simulation.activateSusc = 0;
elseif (expControl.anatomicalID == 6) && ...
        expControl.simulation.activateSusc && ...
        contains(lower(pulseSeqFamilyName), 'se') && ...
        expControl.advancedNotifications
    msg = ['Susceptibility is not currently available for SE-based sequences.',...
        ' This feature will be deactivated for this experiment.'];
    fprintf(1, '\nWARNING: %s\n',msg);
    eduTool.frontend.updateExperimentProgress(expControl,'','info',msg)
    expControl.simulation.activateSusc = 0;
end

%% grid size
if ~expControl.model.isotropicGrid
    % use native model resolution
    expControl.model.gridStep(:) = anatomicalModel.resolution(:);
end

%%
if ~isfield(expControl,'connLocalDB') || isempty(expControl.connLocalDB)
    
    tissueValues = expControl.tissueData.tissueValues;
    
    % assign (with correct number of tissues)
    anatomicalModel.tissueValues(1:size(tissueValues,1),:) = tissueValues(:,1:6);
    
    % deactivate CS if needed
    if expControl.simulation.activateCS == 0
        anatomicalModel.tissueValues(:,4) = 0.0;
    end
    
    % set PD of zero T1/T2 to zero
    anatomicalModel.tissueValues( ...
        (anatomicalModel.tissueValues(:,1) < 1e-4) ...
        |(anatomicalModel.tissueValues(:,2) < 1e-4), 3 ) = 0;
    
    % Compute Susceptibility map for MIDA model
    if (expControl.anatomicalID == 6) && ...
            expControl.simulation.activateSusc && ...
            ~contains(lower(pulseSeqFamilyName), 'se')

        % Take susc. values from db
        tissueSusc = anatomicalModel.susceptibility.X_Water*...
            ones(size(anatomicalModel.tissueValues,1),1);
        tissueSusc(1:size(tissueValues,1),1) = tissueValues(:,7);
        bi = tissueSusc(anatomicalModel.tissueType)*10^-6;
        bi = reshape(bi,anatomicalModel.dimensions);

        % Numerical fft method:
        fftsz   = 512;
        bi      = simModel.susceptibility.getBZfromSuscPaddedFFT(bi, ...
            expControl.mrSystem.b0, anatomicalModel.susceptibility.X_Water, ...
            fftsz); % Head is currently surrounded by water

        anatomicalModel.b0Inhomogeneity = bi(:);
    end
    

else
    
    %% query for tissue properties
    sqlQuery = ['SELECT tissue_id,tissue_name,T1,T2,PD,CS,density,electric_sigma,susc FROM ',...
        'edt_tool_local.anatomical_model_tissues WHERE',...
        ' course_id=',num2str(expControl.courseID),' ORDER BY tissue_id ASC'];
    sqlQueryResults = exec(expControl.connLocalDB, sqlQuery);
    sqlQueryResults = fetch(sqlQueryResults);
    anamTissueData  = sqlQueryResults.Data;
    if ~strcmp(anamTissueData,'No Data')
        
        tissueValues = zeros(max(cell2mat(anamTissueData(:,1))),6);
        
        % get the new Tissue properties: convert from ms to s in T1/T2
        tissueValues(cell2mat(anamTissueData(:,1)),1:6) = [...
            cell2mat(anamTissueData(:,3:4))/1e3,...
            cell2mat(anamTissueData(:,5:6)),...
            cell2mat(anamTissueData(:,7:8))];
        
        % deactivate CS is needed
        if expControl.simulation.activateCS == 0
            tissueValues(:,4) = 0.0;
        end
        
        % set PD of zero T1/T2 to zero
        tissueValues((tissueValues(:,1) < 1e-4)|(tissueValues(:,2) < 1e-4), 3) = 0;
        
        % Compute Susceptibility map for MIDA model
        if (expControl.anatomicalID == 6) && ...
                expControl.simulation.activateSusc && ...
                ~contains(lower(pulseSeqFamilyName), 'se')

            % Take susc. values from db
            tissueSusc = anatomicalModel.susceptibility.X_Water*ones(max(cell2mat(anamTissueData(:,1))),1);
            tissueSusc(cell2mat(anamTissueData(:,1)),1) = cell2mat(anamTissueData(:,9));
            bi = tissueSusc(anatomicalModel.tissueType)*10^-6;
            bi = reshape(bi,anatomicalModel.dimensions);

            % Numerical fft method:
            fftsz   = 512;
            bi      = simModel.susceptibility.getBZfromSuscPaddedFFT(bi, ...
                expControl.mrSystem.b0, anatomicalModel.susceptibility.X_Water, ...
                fftsz); % Head is currently surrounded by water
            
            anatomicalModel.b0Inhomogeneity = bi(:);
        end

		% assign (with correct number of tissues)
		anatomicalModel.tissueValues(1:size(tissueValues,1),:) = tissueValues(:,:);
        
    end
end

%% final message
if expControl.debug.debugMode
    tTotal = toc(tTotal);
    fprintf(fid, '\n%s : done for experiment %d',...
        functionName, expControl.experimentID);
    fprintf(fid, '\n  Elapsed Time      %.3fs', tTotal);
    fprintf(fid, '\n');
    if fid ~=1
        fclose(fid);
    end
end

