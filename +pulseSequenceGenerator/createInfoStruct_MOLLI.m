function info = createInfoStruct_MOLLI(kspace,FOV,dt,receiverBW,TR,TE,...
    kspace_times,f0,sliceThickness,receiver_phase,RFmatrix,...
    acquiredkspace,acquiredFOV,encodedkSpace,infoMOLLI,TI)

acquiredLines   = encodedkSpace.lines;
acquiredRO      = encodedkSpace.columns;

info.pulseSequence.Nx       = kspace(1,1);
info.pulseSequence.Ny       = kspace(1,2);
info.pulseSequence.FOVx     = FOV(1,1);
info.pulseSequence.FOVy     = FOV(1,2);
info.pulseSequence.dt       = dt;
info.pulseSequence.rBW      = receiverBW;
info.pulseSequence.TI       = TI;

info.pulseSequence.RepetitionTime   = TR;
info.pulseSequence.EchoTime         = TE;

info.pulseSequence.RF.duration      = RFmatrix(1,1);
info.pulseSequence.RF.cycles        = RFmatrix(2,1);
info.pulseSequence.RF.angle         = RFmatrix(3,1);
info.pulseSequence.SliceThickness   = RFmatrix(4,1);

info.pulseSequence.kspace           = infoMOLLI.pulseSequence.kspace;
info.reconstruction.imagesPerSlice  = size(TI,2);

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
    ((acquiredkspace(1,1)/receiverBW/dt)-1)/2; %(kspace_times(1,2)-kspace_times(1,1)+1)/2;

info.reconstruction.ismrmrd.header.encoding.encodingLimits.kspace_encoding_step_1.minimum   = 0;
info.reconstruction.ismrmrd.header.encoding.encodingLimits.kspace_encoding_step_1.maximum   = acquiredLines-1;
info.reconstruction.ismrmrd.header.encoding.encodingLimits.kspace_encoding_step_1.center    = (acquiredkspace(1,2)-1)/2;

info.reconstruction.ismrmrd.header.encoding.encodingLimits.repetition.minimum   = 0;
info.reconstruction.ismrmrd.header.encoding.encodingLimits.repetition.maximum   = acquiredLines-1;
info.reconstruction.ismrmrd.header.encoding.encodingLimits.repetition.center    = (acquiredkspace(1,2)-1)/2;

% If reconstruction3D is equal to 1, a 3D FT will be utilized. If it is 
% equal to 0, a 2D FT will be utilized per slice (in cases of a slab that 
% includes several slices).  
info.reconstruction.reconstruction3D = 0;