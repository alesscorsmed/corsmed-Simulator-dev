vps             = 4;
multisliceMode  = 'interleaved'; % 'interleaved' or 'sequential'
slices          = 5;
kspaceLines     = 128;
HR              = 60; % [bpm]
TR              = 0.0054;

cardiacCycleDuration = 60/HR; % [sec]

if strcmp(multisliceMode,'sequential')
    TIRL        = ceil(kspaceLines * slices * cardiacCycleDuration / vps);
    numOfPhases = floor(min(cardiacCycleDuration/(vps*TR),20));
else
    slicesPerGroup = slices;
    
    % while number of phases is less than 10, split the slices into groups
    while min(cardiacCycleDuration/(vps*TR*slicesPerGroup),20)<10
        % A specific number of slices can be acquired per cardiac cycle. If
        % more, then the remaining slices can be acquired in the next pass.
        slicesPerGroup  = round(slicesPerGroup/2);
        if slicesPerGroup == 1
            break;
        end
    end
    
    numOfPhases     = floor(min(cardiacCycleDuration/(vps*TR*slicesPerGroup),20));
    groupsOfSlices  = round(slices/slicesPerGroup);    
    TIRL            = ceil(kspaceLines * groupsOfSlices * cardiacCycleDuration / vps);
end

disp(['TIRL: ',num2str(TIRL),'sec'])
disp(['Num. of phases: ',num2str(numOfPhases)])