function [solution] = runSDW( solution, ...
    sequence, simModel, motionModel, expControl, GPUindex)
%
%
% SIMULATOR.DIFFUSION.RUNSDW
%
%	Runs a Self-Diffusion Weighted Bloch Simulation, calling CUDA kernels.
%
% INPUT
%   solution        solution struct with initial data
%   sequence        sequence to simulate
%   simModel           structure with slice data from spinModel 
%   motionModel     struct with motion model
%   expControl      experiment control struct
%
% OUTPUT
%   solution        solution struct with filled data
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'simulator.diffusion.kernel.runSDW';
if (nargin < 4)
    ME = MException('simulator:wrongArgCount',...
        '%s : wrong argument count',functionName);
    throw(ME);
end

%% open debugging
if expControl.debug.debugMode
    try % open file if possible, otherwise dump to stdout
        fid = fopen(expControl.debug.debugFile,'a');
    catch
        fid = 1;
    end
    tTotal = tic();
    fprintf(fid, '\n%s : start', functionName);
end

%% select current GPU
currentGPU = gpuDevice;

%% extract sequence data and into the GPU
try
    timeDiff    = gpuArray(single(sequence.timeDiff));
    if isfield(sequence,'gdwSignal') && ~isempty(sequence.gdwSignal)
        % update gradients with diffusion encoding
        gxSignal = gpuArray(single( sequence.gxSignal(:) ...
            + sequence.gdwSignal(:,1) ));
        gySignal = gpuArray(single( sequence.gySignal(:)...
            + sequence.gdwSignal(:,2) ));
        gzSignal = gpuArray(single( sequence.gzSignal(:)...
            + sequence.gdwSignal(:,3) ));
    else %% no diffusion encodings
        gxSignal = gpuArray(single( sequence.gxSignal(:) ));
        gySignal = gpuArray(single( sequence.gySignal(:) ));
        gzSignal = gpuArray(single( sequence.gzSignal(:) ));
    end
    rfmSignal   = gpuArray(single(sequence.rfmSignal(:)));
    rfpSignal   = gpuArray(single(sequence.rfpSignal(:)));
    rffSignal   = gpuArray(single(sequence.rffSignal(:)));
    rxSignal    = gpuArray(uint32(sequence.rxSignal(:)));
    swcSignal   = gpuArray(uint32(sequence.swcSignal(:)));
    gamma       = gpuArray(single(sequence.gamma));
    numSteps    = uint32(sequence.numSteps);
    numRxs      = int32(sequence.numRxs);
    seqType     = sequence.type;
catch
    ME = MException('simulator:transfToGPU',...
        '%s : error transfering sequence to GPU',functionName);
    throw(ME);
end

%% convert model to GPU data
try
    % indexes of voxels to simulate (for partition, if needed)
    if isempty(solution.indexes)
        indexes = 1:simModel.numIsochromats;
        numIso  = simModel.numIsochromats;
    else
        indexes = solution.indexes;
        numIso  = length(indexes);
    end
    % convert to correct data types
    numIso          = uint32(numIso);
    b0              = single(simModel.b0);
    mu              = single(simModel.mu);
    
    x               = gpuArray(single(simModel.x(indexes)));
    y               = gpuArray(single(simModel.y(indexes)));
    z               = gpuArray(single(simModel.z(indexes)));
    
    bi              = gpuArray(single(simModel.bi(indexes)));
    pd              = gpuArray(single(simModel.pd(indexes)));
    
    tissueValues    = gpuArray(single(transpose(simModel.tissueValues)));
    tissueType      = gpuArray(uint32(simModel.tissueType(indexes)));
    
    % generate the diffusion arrays
    simModel.xDiffusion = simModel.tissueDiff(simModel.tissueType,1);
    simModel.yDiffusion = simModel.tissueDiff(simModel.tissueType,2);
    simModel.zDiffusion = simModel.tissueDiff(simModel.tissueType,3);
    
    % transfer
    xDiff           = gpuArray(single(simModel.xDiffusion(indexes)));
    yDiff           = gpuArray(single(simModel.yDiffusion(indexes)));
    zDiff           = gpuArray(single(simModel.zDiffusion(indexes)));
    
    coilMapsX       = gpuArray(single(simModel.rxCoilMapsX(indexes,:)));
    coilMapsY       = gpuArray(single(simModel.rxCoilMapsY(indexes,:)));
    nCoils          = gpuArray(int32(simModel.numRxCoils));
    
    voxLx           = gpuArray(single(simModel.resolution(1)));
    voxLy           = gpuArray(single(simModel.resolution(2)));
    voxLz           = gpuArray(single(simModel.resolution(3)));
catch
    ME = MException('simulator:transfToGPU',...
        '%s : error transfering model to GPU',functionName);
    throw(ME);
end

%% convert motion model to GPU data
try
    % indexes of voxels to simulate (for partition, if needed)
    if isempty(motionModel) || strcmpi(motionModel.type,'none')
        %% no motion
        applyMotion     = int32(0);
        angleRotXY      = gpuArray(single(0));
        angleRotXZ      = gpuArray(single(0));
        angleRotYZ      = gpuArray(single(0));
        xCenterRotXY    = gpuArray(single(0));
        yCenterRotXY    = gpuArray(single(0));
        xCenterRotXZ    = gpuArray(single(0));
        zCenterRotXZ    = gpuArray(single(0));
        yCenterRotYZ    = gpuArray(single(0));
        zCenterRotYZ    = gpuArray(single(0));
        transX          = gpuArray(single(0));
        transY          = gpuArray(single(0));
        transZ          = gpuArray(single(0)); 
    else
        %% apply motion
        applyMotion     = int32(1);
        angleRotXY      = gpuArray(single(motionModel.angleRotXY));
        angleRotXZ      = gpuArray(single(motionModel.angleRotXZ));
        angleRotYZ      = gpuArray(single(motionModel.angleRotYZ));
        xCenterRotXY    = gpuArray(single(motionModel.xCenterRotXY));
        yCenterRotXY    = gpuArray(single(motionModel.yCenterRotXY));
        xCenterRotXZ    = gpuArray(single(motionModel.xCenterRotXZ));
        zCenterRotXZ    = gpuArray(single(motionModel.zCenterRotXZ));
        yCenterRotYZ    = gpuArray(single(motionModel.yCenterRotYZ));
        zCenterRotYZ    = gpuArray(single(motionModel.zCenterRotYZ));
        transX          = gpuArray(single(motionModel.transX));
        transY          = gpuArray(single(motionModel.transY));
        transZ          = gpuArray(single(motionModel.transZ));
    end    
catch
   ME = MException('simulator:transfToGPU',...
       '%s : error transfering motion model to GPU',functionName);
   throw(ME);
end

%% Allocate space for magnetization (initialize) and Signal solution
try
    Mm = gpuArray(single(solution.Mm));
    Mp = gpuArray(single(solution.Mp));
    Mz = gpuArray(single(solution.Mz));
    
    dMpDx = gpuArray(single(solution.dMpDx));
    dMpDy = gpuArray(single(solution.dMpDy));
    dMpDz = gpuArray(single(solution.dMpDz));
    
    Sx = gpuArray(single(solution.Sx));
    Sy = gpuArray(single(solution.Sy));
    Sz = gpuArray(single(solution.Sz));
catch
    ME = MException('simulator:transfToGPU',...
        '%s : error transfering solution to GPU',functionName);
    throw(ME);
end

%% check status
if strcmp(expControl.application,'edutool') && ...
        isfield(expControl, 'connLocalDB') ...
        && ~isempty(expControl.connLocalDB)
    eduTool.frontend.checkExperimentStatus(expControl,expControl.experimentID)
end

%% prepare to run the kernels
kernelPtx  = expControl.simulation.kernelPtx;
numThreads = expControl.simulation.threads;
numBlocks  = ceil(single(numIso)/numThreads);
% kernel proto is common to all cases
kernelProto = [ 'float *, float *, float *, ',... % Mm, Mp, Mz
                'float *, float *, float *, ',... % dMpdX,dMpdY,dMpdZ
                'float *, float *, float *, ',... % Sx, Sy, Sz
                'float *, unsigned int *, unsigned int *, ',... % Tdiff, RX, SWC
                'float *, float *, float *, ',... % pulse Gx, Gy and Gz
                'float *, float *, float *, ',... % pulse RFmag, RFphase and RFfreq
                'float *, float *, float *, '... % x, y, z
                'float *, float *, ',... % Bi, PD
                'float *, unsigned int *, '... % tissueValues, tissueType
                'float *, float *, float *, '... % diffusion for x, y, z
                'int, float *, float *, ',... % nCoils, Cx, Cy
                'unsigned int, unsigned int, unsigned int, ',... % nIso, tStart, tEnd
                'float, float, float, ',... % gamma, b0, mu
                'float, float, float, ',... % dx, dy, dz
                'int, ',... % applyMotion
                'float *, float, float, '... % angleRotXY, xCenterRotXY, yCenterRotXY
                'float *, float, float, '... % angleRotXZ, xCenterRotXZ, zCenterRotXZ
                'float *, float, float, '... % angleRotYZ, yCenterRotYZ, zCenterRotYZ
                'float *, float*, float*' ]; % transX, transY, transZ

%% RF kernel -- kernel for RF parts
kernelFunction = 'swdExpRF';
% prepare the kernel
KRF = parallel.gpu.CUDAKernel(kernelPtx, kernelProto, kernelFunction);
KRF.ThreadBlockSize = numThreads;
KRF.GridSize = numBlocks;

%% GR kernel -- kernel for only Gradients (no Readouts, no SWC)
kernelFunction = 'swdExpGR';
% prepare the kernel
KGR = parallel.gpu.CUDAKernel(kernelPtx, kernelProto, kernelFunction);
KGR.ThreadBlockSize = numThreads;
KGR.GridSize = numBlocks;

%% RO kernel -- gradient kernel with potential readouts and/or SWC
kernelFunction = 'swdExpRO';
% prepare the kernel
KRO = parallel.gpu.CUDAKernel(kernelPtx, kernelProto, kernelFunction);
KRO.ThreadBlockSize = numThreads;
KRO.GridSize = numBlocks;

%% loop on sequence parts and run the kernels
if  expControl.debug.debugMode
    fprintf(fid, '\nStaring Phasor simulation using Kernel File %s',kernelPtx);
    % for collecting data
    tRF = 0.0; nRF = 0;
    tGR = 0.0; nGR = 0;
    tRO = 0.0; nRO = 0;
    tCUDA = 0.0;
    
    % back end progress bar
    if expControl.debug.waitBarBE
        WB = waitbar(0, sprintf('Diffusion Simulation: %d%%', 0 ), ...
            'Name','corsmed', 'CreateCancelBtn','setappdata(gcbf,''stop'',1)');
        setappdata(WB,'stop',0);
    end
end

%% check status
if strcmp(expControl.application,'edutool') && ...
        isfield(expControl, 'connLocalDB') ...
        && ~isempty(expControl.connLocalDB)
    eduTool.frontend.checkExperimentStatus(expControl,expControl.experimentID)
end

for ss = 1:sequence.numParts
    % define type or part to simulate
    partType = sequence.partType{ss};
    % part limits
    timeIni = sequence.partLimits(ss,1); 
    timeEnd = sequence.partLimits(ss,2);
    % reset resolutions: during GR and RO they are voxel size
    dx = voxLx;
    dy = voxLy;
    dz = voxLz;
    switch lower(partType)
        case 'rf'
            % use standard RF kernel
            K = KRF;
            
            % estimated of phase accrued during part
            accPhaseX = sum(gxSignal(timeIni:timeEnd).*timeDiff(timeIni:timeEnd));
            accPhaseY = sum(gySignal(timeIni:timeEnd).*timeDiff(timeIni:timeEnd));
            accPhaseZ = sum(gzSignal(timeIni:timeEnd).*timeDiff(timeIni:timeEnd));
            % add accrued phase in voxel due to current derivatives of phase
            accPhaseX = 2*pi*gamma*abs(accPhaseX) + max(abs(dMpDx));
            accPhaseY = 2*pi*gamma*abs(accPhaseY) + max(abs(dMpDy));
            accPhaseZ = 2*pi*gamma*abs(accPhaseZ) + max(abs(dMpDz));
            
            % modify the voxel resolution to improve accuracy
            % scale so that maximum phase accrued in spatial delta is 0.5*PI
            if accPhaseX > 0
                dx = min(0.5*pi/accPhaseX, 5e-5);
            end
            if accPhaseY > 0
                dy = min(0.5*pi/accPhaseY, 5e-5);
            end
            if accPhaseZ > 0
                dz = min(0.5*pi/accPhaseZ, 5e-5);
            end
        case 'gr'
            K = KGR;
        case 'ro'
            K = KRO;
        otherwise
            ME = MException('simulator:kernelPrep',...
                '%s : sequence type %s not supported',functionName,partType);
            throw(ME);
    end

    % part limits: move it to zero based for CUDA
    timeIni = uint32(sequence.partLimits(ss,1)-1); 
    timeEnd = uint32(sequence.partLimits(ss,2)-1);
    
    tkernel = tic();
    % launch kernel
    [Mm, Mp, Mz, Sx, Sy, Sz, dMpDx, dMpDy, dMpDz ] = feval( K, ...
        Mm, Mp, Mz, Sx, Sy, Sz, dMpDx, dMpDy, dMpDz, ...
        timeDiff, rxSignal, swcSignal, ...
        gxSignal, gySignal, gzSignal, ...
        rfmSignal, rfpSignal, rffSignal, ...
        x, y, z, bi, pd, tissueValues, tissueType, ...
        xDiff, yDiff, zDiff,...
        nCoils, coilMapsX, coilMapsY,...
        numIso, timeIni, timeEnd, gamma, b0, mu, ...
        dx, dy, dz, ...
        applyMotion, ...
        angleRotXY, xCenterRotXY, yCenterRotXY, ...
        angleRotXZ, xCenterRotXZ, zCenterRotXZ, ...
        angleRotYZ, yCenterRotYZ, zCenterRotYZ, ...
        transX, transY, transZ );
    
    % force to wait for feval (asynchronous)
    wait(currentGPU);
    tkernel = toc(tkernel);
        
    if  expControl.debug.debugMode
                
        % back end progress bar
        if expControl.debug.waitBarBE
            if getappdata(WB,'stop')
                delete(WB);
                clear WB;
                ME = MException('simulator:Stopped',...
                    '%s : Stopped by user',functionName);
                throw(ME);
            end
            progress = round(single(timeEnd)/single(numSteps)*100);
            waitbar(single(timeEnd)/single(numSteps), WB, ...
                sprintf('Diffusion Simulation: %d%%', progress ));
        end
        %fprintf(fid, '\n  %s part done in %.3fs -- progress %d%%',...
        %    sequence.partType{ss}, tkernel, progress );
        % data collection
        tCUDA = tCUDA + tkernel;
        switch lower(partType)
            case 'rf'
                tRF = tRF + tkernel;
                nRF = nRF + single(timeEnd) - single(timeIni) + 1;
            case 'gr'
                tGR = tGR + tkernel;
                nGR = nGR + single(timeEnd) - single(timeIni) + 1;
            case 'ro'
                tRO = tRO + tkernel;
                nRO = nRO + single(timeEnd) - single(timeIni) + 1;
            otherwise % default
                ME = MException('simulator:kernelPrep',...
                    '%s : sequence type %s not supported (2nd catch)',functionName,partType);
                throw(ME);
        end
    end
    
end

%% collect data from GPU
solution.Mm = gather(Mm);
solution.Mp = gather(Mp);
solution.Mz = gather(Mz);
solution.Mx = solution.Mm.*cos(solution.Mp);
solution.My = solution.Mm.*sin(solution.Mp);
solution.dMpDx = gather(dMpDx);
solution.dMpDy = gather(dMpDy);
solution.dMpDz = gather(dMpDz);
solution.Sx = gather(Sx);
solution.Sy = gather(Sy);
solution.Sz = gather(Sz);

%% report
if  expControl.debug.debugMode
    tTotal = toc(tTotal);
    fprintf(fid, '\n%s : done for sequence type %s, %d parts, IRL time %.3fs',...
        functionName, seqType, sequence.numParts,sequence.time(end));
    numIso = double(numIso);
    fprintf(fid, '\n  Current GPU       %d (%s)', ...
        currentGPU.Index, currentGPU.Name);
    fprintf(fid, '\n  # GPU Blocks      %d', numBlocks);
    fprintf(fid, '\n  # GPU Threads     %d', numThreads);
    fprintf(fid, '\n  # Isochromas      %d', numIso);
    fprintf(fid, '\n  # Time Steps      %d', numSteps);
    fprintf(fid, '\n  # Readouts        %d', numRxs);
    fprintf(fid, '\n  Elapsed Time      %.3fs / %d steps (%.1fps/voxel/step)',...
        tTotal, numSteps, 1e12*tTotal/numIso/numSteps);
    fprintf(fid, '\n  Kernel  Time      %.3fs / %d steps (%.1fps/voxel/step)',...
        tCUDA, numSteps, 1e12*tCUDA/numIso/numSteps);
    fprintf(fid, '\n  RF part Time      %.3fs / %d steps (%.1fps/voxel/step)',...
        tRF, nRF, 1e12*tRF/numIso/max(nRF,1));
    fprintf(fid, '\n  GR part Time      %.3fs / %d steps (%.1fps/voxel/step)',...
        tGR, nGR, 1e12*tGR/numIso/max(nGR,1));
    fprintf(fid, '\n  RO part Time      %.3fs / %d steps (%.1fps/voxel/step)',...
        tRO, nRO, 1e12*tRO/numIso/max(nRO,1));
    fprintf(fid, '\n');
    if fid ~=1
        fclose(fid);
    end
    % clear progress bar
    if expControl.debug.waitBarBE
        delete(WB);
        clear WB;
    end
end
