function val = simulator3pNtissue(x,auxdata)

%% load the ground-truth kspace
kSpaceGT        = auxdata.simulation.kSpaceGT;
pulseSequence   = auxdata.simulation.pulseSequence;
spinModel       = auxdata.simulation.spinModel;
motionModel     = auxdata.simulation.motionModel;
expControl      = auxdata.simulation.expControl;

T1val = x(1:3:end);
T2val = x(2:3:end);
PDval = x(3:3:end);

spinModel.tissueValues  = [T1val(:),T2val(:),PDval(:),zeros(size(PDval(:),1),3)];
spinModel.pd            = transpose(PDval(spinModel.tissueType));

% Run the simulation using the auxdata.simulation.pulseSequence and the
% modified auxdata.simulation.spinModel
GPUindex = 1;
tic
timeSolution    = sbr.run.sbrEngine(spinModel,pulseSequence,...
    motionModel,expControl,GPUindex);
toc

timeSolution.noise = zeros(size(timeSolution.Sy));

% Reshape the outcome of the simulator
[kSpace] = reconstructor.signal.mapKspace( timeSolution,...
    auxdata.simulation.encoding, expControl );

val = sum(abs((kSpace(:)-kSpaceGT(:))./max(kSpaceGT(:))).^2);