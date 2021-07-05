function TIRL = CINEcalculateTIRL(acqData, anatomicalModel, pulseSequence)

cardiacCycleDuration    = 60/anatomicalModel.HR; % [sec]
kspaceLines             = pulseSequence.numEnc;

if isfield(acqData,'interleaveSlice') && ...
        acqData.interleaveSlice
    
    slicesPerGroup = acqData.numSlices;
    
    % while number of phases is less than 10, split the slices into groups
    while min(cardiacCycleDuration/...
            (acqData.vps*acqData.TR*slicesPerGroup),...
            anatomicalModel.numPhases)<10
        % A specific number of slices can be acquired per cardiac cycle. If
        % more, then the remaining slices can be acquired in the next pass.
        slicesPerGroup  = round(slicesPerGroup/2);
        if slicesPerGroup == 1
            break;
        end
    end
    
    groupsOfSlices  = round(acqData.numSlices/slicesPerGroup);    
    TIRL            = ceil(kspaceLines * groupsOfSlices * cardiacCycleDuration / acqData.vps);
    
else
    
    TIRL        = ceil(kspaceLines * acqData.numSlices * cardiacCycleDuration / acqData.vps);
    
end
    
fprintf(1, '\nTIRL: %.3fs', TIRL)