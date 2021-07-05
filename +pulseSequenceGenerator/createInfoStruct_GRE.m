function info = createInfoStruct_GRE(kspace,FOV,dt,receiverBW,TR,TE,...
    kspace_times,f0,sliceThickness,acquiredkspace,acquiredFOV,...
    encodedkSpace)

acquiredLines   = encodedkSpace.lines;
acquiredRO      = encodedkSpace.columns;
numCO           = encodedkSpace.contrasts;

info.pulseSequence.Nx       = kspace(1,1);
info.pulseSequence.Ny       = kspace(1,2);
info.pulseSequence.FOVx     = FOV(1,2);
info.pulseSequence.FOVy     = FOV(1,2);
info.pulseSequence.dt       = dt;
info.pulseSequence.rBW      = receiverBW;

info.pulseSequence.RepetitionTime   = TR;
info.pulseSequence.EchoTime         = TE;

%% Structure the info.pulseSequence.kspace
% Old code
% info_kspace = [ones(size(kspace_times,1),14),zeros(size(kspace_times,1),2)];
% info_kspace(:,1) = [1:acquiredLines]';
% info_kspace(:,2) = kspace_times(:,1);
% info_kspace(:,3) = kspace_times(:,2)-kspace_times(:,1)+1;
% info_kspace(:,4) = zeros(size(info_kspace(:,4)));
% info_kspace(:,5) = kspace_times(:,2)-kspace_times(:,1)+1;
% info_kspace(:,6) = [1:acquiredLines]';

% New code
info_kspace = [ones(acquiredLines*numCO,14),...
    zeros(acquiredLines*numCO,2)];

info_kspace(:,1)    = [1:acquiredLines*numCO]';

kspace_times_new1   = permute(kspace_times(:,1,:),[3,1,2]);
kspace_times_new2   = permute(kspace_times(:,2,:),[3,1,2]);
info_kspace(:,2)    = kspace_times_new1(:);

info_kspace(:,3)    = kspace_times_new2(:)-kspace_times_new1(:)+1;
info_kspace(:,4)    = zeros(size(info_kspace(:,4)));
info_kspace(:,5)    = kspace_times_new2(:)-kspace_times_new1(:)+1;

kspaceEncodeStep1   = repmat(1:acquiredLines,numCO,1);
kspaceEncodeStep1   = kspaceEncodeStep1(:);
info_kspace(:,6)    = kspaceEncodeStep1;

kspace_contrasts    = 1:numCO;
kspace_contrasts    = repmat(kspace_contrasts,acquiredLines,1)';
info_kspace(:,10)   = kspace_contrasts(:);

info.pulseSequence.kspace = info_kspace;

%%
info.reconstruction.ismrmrd.header.experimentalConditions.H1resonanceFrequency_Hz = f0;

info.reconstruction.ismrmrd.header.encoding.trajectory = 'cartesian';
info.reconstruction.ismrmrd.header.encoding.encodedSpace.fieldOfView_mm.x = acquiredFOV(1,1)*1000;
info.reconstruction.ismrmrd.header.encoding.encodedSpace.fieldOfView_mm.y = acquiredFOV(1,2)*1000;
info.reconstruction.ismrmrd.header.encoding.encodedSpace.fieldOfView_mm.z = sliceThickness*1000;

info.reconstruction.ismrmrd.header.encoding.encodedSpace.matrixSize.x = acquiredRO;
info.reconstruction.ismrmrd.header.encoding.encodedSpace.matrixSize.y = acquiredLines;
info.reconstruction.ismrmrd.header.encoding.encodedSpace.matrixSize.z = 1;

info.reconstruction.ismrmrd.header.encoding.reconSpace.fieldOfView_mm.x = FOV(1,1)*1000;
info.reconstruction.ismrmrd.header.encoding.reconSpace.fieldOfView_mm.y = FOV(1,2)*1000;
info.reconstruction.ismrmrd.header.encoding.reconSpace.fieldOfView_mm.z = sliceThickness*1000;

info.reconstruction.ismrmrd.header.encoding.reconSpace.matrixSize.x = kspace(1,1);
info.reconstruction.ismrmrd.header.encoding.reconSpace.matrixSize.y = kspace(1,2);
info.reconstruction.ismrmrd.header.encoding.reconSpace.matrixSize.z = 1;

info.reconstruction.ismrmrd.header.encoding.encodingLimits.kspace_encoding_step_0.minimum   = 0;
info.reconstruction.ismrmrd.header.encoding.encodingLimits.kspace_encoding_step_0.maximum   = ...
    acquiredRO - 1;
info.reconstruction.ismrmrd.header.encoding.encodingLimits.kspace_encoding_step_0.center    = ...
    ((acquiredkspace(1,1)/receiverBW/dt)-1)/2;  %(kspace_times(1,2)-kspace_times(1,1)+1)/2;

info.reconstruction.ismrmrd.header.encoding.encodingLimits.kspace_encoding_step_1.minimum   = 0;
info.reconstruction.ismrmrd.header.encoding.encodingLimits.kspace_encoding_step_1.maximum   = acquiredLines-1;
info.reconstruction.ismrmrd.header.encoding.encodingLimits.kspace_encoding_step_1.center    = (acquiredkspace(1,2)-1)/2;

info.reconstruction.ismrmrd.header.encoding.encodingLimits.repetition.minimum   = 0;
info.reconstruction.ismrmrd.header.encoding.encodingLimits.repetition.maximum   = acquiredLines-1;
info.reconstruction.ismrmrd.header.encoding.encodingLimits.repetition.center    = (acquiredkspace(1,2)-1)/2;

% If reconstruction3D is equal to 1, a 3D FT will be utilized. If it is 
% equal to 0, a 2D FT will be utilized per slice (in cases of a slab that 
% includes several slices).  
info.reconstruction.reconstruction3D = 1;