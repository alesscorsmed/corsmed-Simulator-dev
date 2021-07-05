function [phasesIds, numOfPhases] = CINEcalculatePhases(acquisition, anatomicalModel)

cardiacCycleDuration = 60/anatomicalModel.HR; % [sec]

if isfield(acquisition.data,'interleaveSlice') && ...
        acquisition.data.interleaveSlice
    
    slicesPerGroup = acquisition.data.numSlices;
    
    % while number of phases is less than 10, split the slices into groups
    while min(cardiacCycleDuration/...
            (acquisition.data.vps*acquisition.data.TR*slicesPerGroup),...
            anatomicalModel.numPhases)<10
        % A specific number of slices can be acquired per cardiac cycle. If
        % more, then the remaining slices can be acquired in the next pass.
        slicesPerGroup  = round(slicesPerGroup/2);
        if slicesPerGroup == 1
            break;
        end
    end
    
    numOfPhases     = floor(min(cardiacCycleDuration/...
        (acquisition.data.vps*acquisition.data.TR*slicesPerGroup),...
        anatomicalModel.numPhases));
    
%     groupsOfSlices  = round(acquisition.data.numSlices/slicesPerGroup);    
%     TIRL            = ceil(kspaceLines * groupsOfSlices * cardiacCycleDuration / vps);
    
else
    
    numOfPhases = floor(min(cardiacCycleDuration/...
        (acquisition.data.vps*acquisition.data.TR),anatomicalModel.numPhases));
%     TIRL        = ceil(kspaceLines * slices * cardiacCycleDuration / vps);
end
    
fprintf(1, '\nNum. of phases: %d', numOfPhases)
phasesIds = floor(linspace(1,anatomicalModel.numPhases,numOfPhases));

%     
%     
%     
%     
%     
%     
% 
% vps             = 4;
% multisliceMode  = 'interleaved'; % 'interleaved' or 'sequential'
% slices          = 5;
% kspaceLines     = 128;
% HR              = 60; % [bpm]
% TR              = 0.0054;
% 
% cardiacCycleDuration = 60/HR; % [sec]
% 
% if strcmp(multisliceMode,'sequential')
%     TIRL        = ceil(kspaceLines * slices * cardiacCycleDuration / vps);
%     numOfPhases = floor(min(cardiacCycleDuration/(vps*TR),20));
% else
%     slicesPerGroup = slices;
%     
%     % while number of phases is less than 10, split the slices into groups
%     while min(cardiacCycleDuration/(vps*TR*slicesPerGroup),20)<10
%         % A specific number of slices can be acquired per cardiac cycle. If
%         % more, then the remaining slices can be acquired in the next pass.
%         slicesPerGroup  = round(slicesPerGroup/2);
%         if slicesPerGroup == 1
%             break;
%         end
%     end
%     
%     numOfPhases     = floor(min(cardiacCycleDuration/(vps*TR*slicesPerGroup),20));
%     groupsOfSlices  = round(slices/slicesPerGroup);    
%     TIRL            = ceil(kspaceLines * groupsOfSlices * cardiacCycleDuration / vps);
% end
% 
% disp(['TIRL: ',num2str(TIRL),'sec'])
% disp(['Num. of phases: ',num2str(numOfPhases)])