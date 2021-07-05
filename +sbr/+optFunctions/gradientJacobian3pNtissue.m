function Jac_d = gradientJacobian3pNtissue(x,auxdata)
%
% sbr.run.sbrSlicer
%
%	Creates the anatomical model with tissuesX x tissuesY tissues within
%	the single-slice anatomical model.
%
% INPUT
%   gridStep is a [1x3] array representing the size of the element in 3D
%
%
% OUTPUT
%

%%
gradStepT1  = auxdata.gradient.stepT1;
gradStepT2  = auxdata.gradient.stepT2;
gradStepPD  = auxdata.gradient.stepPD;

%% load the ground-truth kspace
kSpaceGT    = auxdata.simulation.kSpaceGT;

%% create the diagonal matrix
GRstepV0    = [gradStepT1,gradStepT2,gradStepPD];
GRstepV1    = repmat(GRstepV0,[1,fix(length(x)/3)]);
GRstep      = diag(GRstepV1);

%% Identify the size of the parpool
poolobj = gcp('nocreate'); % If no pool, do not create new one.
if isempty(poolobj)
    poolsize = 1;
else
    poolsize = poolobj.NumWorkers;
end

%% One-side gradient
InValues        = x+GRstep;
val1            = zeros(1,size(InValues,1));
pulseSequence   = auxdata.simulation.pulseSequence;
spinModel       = auxdata.simulation.spinModel;
motionModel     = auxdata.simulation.motionModel;
expControl      = auxdata.simulation.expControl;
encoding        = auxdata.simulation.encoding;
parfor i=1:size(InValues,1)
    
    InValuesGPU     = InValues(i,:);
    
    T1val           = InValuesGPU(1,1:3:end);
    T2val           = InValuesGPU(1,2:3:end);
    PDval           = InValuesGPU(1,3:3:end);
    
    spinModelGPU = spinModel;
    
    % Modify the spinModel based on the T1val, T2val and PDval.
    spinModelGPU.tissueValues  = [T1val(:),T2val(:),PDval(:),zeros(size(PDval(:),1),3)];
    spinModelGPU.pd            = transpose(PDval(spinModelGPU.tissueType));
    
    % Suppress warnings
    warning('off','parallel:gpu:device:DeviceDeprecated')
    
    % Run the simulation using the auxdata.simulation.pulseSequence and the
    % modified auxdata.simulation.spinModel
    GPUindex = mod(i,poolsize)+1;
    tic
    timeSolution    = sbr.run.sbrEngine(spinModelGPU,pulseSequence,...
        motionModel,expControl,GPUindex);
    toc
    
    timeSolution.noise  = zeros(size(timeSolution.Sy));
    
    % Reshape the outcome of the simulator
    [kSpace] = reconstructor.signal.mapKspace( timeSolution,...
        encoding, expControl );
    
    % Find the differences in kspaces
    val1(i) = sum(abs((kSpace(:)-kSpaceGT(:))./max(kSpaceGT(:))).^2);  %#ok
    
end
term1=val1;

%% Opposite-side gradient
if strcmp(auxdata.gradient.side,'double')

    InValues        = x-GRstep;
    val2            = zeros(1,size(InValues,1));    
    pulseSequence   = auxdata.simulation.pulseSequence;
    spinModel       = auxdata.simulation.spinModel;
    motionModel     = auxdata.simulation.motionModel;
    expControl      = auxdata.simulation.expControl;
    parfor i=1:size(InValues,1) 

        InValuesGPU     = InValues(i,:);

        T1val           = InValuesGPU(1,1:3:end);
        T2val           = InValuesGPU(1,2:3:end);
        PDval           = InValuesGPU(1,3:3:end);   

        spinModelGPU = spinModel;

        % Modify the spinModel based on the T1val, T2val and PDval.
        spinModelGPU.tissueValues  = [T1val(:),T2val(:),PDval(:),zeros(size(PDval(:),1),3)];
        spinModelGPU.pd            = transpose(PDval(spinModelGPU.tissueType));

        % Suppress warnings
        warning('off','parallel:gpu:device:DeviceDeprecated')

        % Run the simulation using the auxdata.simulation.pulseSequence and the
        % modified auxdata.simulation.spinModel
        GPUindex = mod(i,poolsize)+1;
        tic
        timeSolution    = sbr.run.sbrEngine(spinModelGPU,pulseSequence,...
            motionModel,expControl,GPUindex);
        toc

        timeSolution.noise = zeros(size(timeSolution.Sy));

        % Reshape the outcome of the simulator
        [kSpace] = reconstructor.signal.mapKspace(timeSolution,...
            encoding, expControl);

        % Find the differences in kspaces
        val2(i) = sum(abs((kSpace(:)-kSpaceGT(:))./max(kSpaceGT(:))).^2);  %#ok
    end
    term2 = val2;
end

%% Calculate the Jacobian
if strcmp(auxdata.gradient.side,'double')
    disp('Double-side gradient has chosen')
    for j = 1:length(x)
        Jac_d(j) = (term1(j)-term2(j))/(2*GRstepV1(j)); %#ok<AGROW>
    end
else
    disp('Single-side gradient has chosen')
    val = sbr.optFunctions.simulator3pNtissue(x,auxdata);
    for j = 1:length(x)
        Jac_d(j) = (term1(j)-val)/(GRstepV1(j)); %#ok<AGROW>
    end
end