function [solution,stats] = runDiffusion( solution, ...
    sequence, simModel, motionModel, simControl, dbgControl )
%
% SPINTWIN.FWDBLOCH.RUNDIFFUSION
%
%	Runs Diffusion Bloch simulation, calling CUDA kernels:
%       Choice of Numerical/Analytical for RF
%       Analytical exponential integration for remainder
%       Allows for single/double precision
%       Phase based signal integration
%       Incorporates intra-voxel diffusion effects
%
% INPUT
%   solution        solution struct with initial data
%   sequence        sequence to simulate
%   simModel        structure with slice data from spinModel 
%   motionModel     struct with motion model
%   simControl      control struct for the simulation
%
% OUTPUT
%   solution        solution struct with filled data
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'spinTwin:fwdBloch:runDiffusion';
% check args
if (nargin < 5) ...
        || isempty(solution) || isempty(simModel) ...
        || isempty(sequence) || isempty(simControl)
    ME = MException(['error:',functionName],...
        'Wrong arguments.');
    throw(ME);
end
% optional arguments
if (nargin < 6) || isempty(dbgControl)
    dbgControl.mode = 0;
    dbgControl.file = [];
end
% info for debugging
if dbgControl.mode
    try % open file if possible, otherwise dump to stdout
        fid = fopen(dbgControl.file,'a');
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
    msg = sprintf( 'Failed to transfer sequence to GPU' );
    ME = MException(['error:',functionName], '%s', msg);
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
    numIso	= uint32(numIso);
    b0      = single(simModel.b0);
    mu      = single(simModel.mu);
    % positions / gradient
    x       = gpuArray(single(simModel.x(indexes)));
    y       = gpuArray(single(simModel.y(indexes)));
    z       = gpuArray(single(simModel.z(indexes)));
    % generate voxel data if needed
    if ~isfield(simModel, 'r1') || isempty(simModel.r1)
        simModel.r1 = reshape(1./simModel.tissueValues(simModel.tissueType(:),1),[],1);
        simModel.r1(isnan(simModel.r1)) = 1e5;
    end
    if ~isfield(simModel, 'r2') || isempty(simModel.r2)
        simModel.r2 = reshape(1./simModel.tissueValues(simModel.tissueType(:),2),[],1);
        simModel.r2(isnan(simModel.r2)) = 1e6 + simModel.r1(isnan(simModel.r2));
    end
    if ~isfield(simModel, 'pd') || isempty(simModel.pd)
        simModel.pd = reshape(simModel.tissueValues(simModel.tissueType(:),3),[],1);
        simModel.pd(isnan(simModel.pd)) = 0.0;
    end
    if ~isfield(simModel, 'cs') || isempty(simModel.cs)
        simModel.cs = reshape(simModel.tissueValues(simModel.tissueType(:),4),[],1);
        simModel.cs(isnan(simModel.cs)) = 0.0;
    end
    if ~isfield(simModel, 'bi') || isempty(simModel.bi)
        simModel.bi = 0*simModel.pd;
    end
    % transfer voxel properties to GPU
    pd      = gpuArray(single(simModel.pd(indexes)));
    r1      = gpuArray(single(simModel.r1(indexes)));
    r2      = gpuArray(single(simModel.r2(indexes)));
    cs      = gpuArray(single(1.5e-6*simModel.cs(indexes)));
    bi      = gpuArray(single(simModel.bi(indexes)));
    % coil sens
    coilMapsX	= gpuArray(single(simModel.rxCoilMapsX(indexes,:)));
    coilMapsY	= gpuArray(single(simModel.rxCoilMapsY(indexes,:)));
    nCoils      = gpuArray(int32(simModel.numRxCoils));
    % voxel volume
    voxLx       = gpuArray(single(simModel.resolution(1)));
    voxLy       = gpuArray(single(simModel.resolution(2)));
    voxLz       = gpuArray(single(simModel.resolution(3)));
    % generate the diffusion arrays
    simModel.xDiffusion = simModel.tissueDiff(simModel.tissueType,1);
    simModel.yDiffusion = simModel.tissueDiff(simModel.tissueType,2);
    simModel.zDiffusion = simModel.tissueDiff(simModel.tissueType,3);
    % transfer
    xDiffusion      = gpuArray(single(simModel.xDiffusion(indexes)));
    yDiffusion      = gpuArray(single(simModel.yDiffusion(indexes)));
    zDiffusion      = gpuArray(single(simModel.zDiffusion(indexes)));
catch
    msg = sprintf( 'Failed to transfer tissue model to GPU' );
    ME = MException(['error:',functionName], '%s', msg);
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
    msg = sprintf( 'Failed to transfer motion model to GPU' );
    ME = MException(['error:',functionName], '%s', msg);
    throw(ME);
end

%% Allocate space for magnetization (initialize) and Signal solution
try
    Mx = gpuArray(single(solution.Mx));
    My = gpuArray(single(solution.My));
    Mz = gpuArray(single(solution.Mz));
    
    dPx = gpuArray(single(solution.dMpDx));
    dPy = gpuArray(single(solution.dMpDy));
    dPz = gpuArray(single(solution.dMpDz));
    
    Sx = gpuArray(single(solution.Sx));
    Sy = gpuArray(single(solution.Sy));
    Sz = gpuArray(single(solution.Sz));
catch
    msg = sprintf( 'Failed to transfer solution to GPU' );
    ME = MException(['error:',functionName], '%s', msg);
    throw(ME);
end

%% prepare to run the kernels
kernelPtx  = simControl.ptxSIM;
numThreads = simControl.threads;
numBlocks  = ceil(single(numIso)/numThreads);
% kernel proto is common to all cases
kernelProto = [ 'float *, float *, float *, ',... % Mx, My, Mz
                'float *, float *, float *, ',... % dPx,dPY,dPz
                'float *, float *, float *, ',... % Sx, Sy, Sz
                'const float *,' ...    % timeDiff
                'const unsigned int *, const unsigned int *,' ... % RX, SWC
                'const float *, const float *, const float *,' ... % pulse Gx, Gy and Gz
                'const float *, const float *, const float *,' ... % pulse RFmag, RFphase and RFfreq
                'const float *, const float *, const float *,' ... % x, y, z
                'const float *, const float *, const float *,' ... % pd r1, r2
                'const float *, const float *,' ... % cs, bi
                'const float *, const float *, const float *, '... % diffusion for x, y, z
                'const int, const float *, const float *,' ... % nCoils, Cx, Cy
                'const unsigned int, const unsigned int, const unsigned int,' ... % nIso, tStart, tEnd
                'const float, const float, const float,' ... % gamma, b0, mu
                'const float, const float, const float,' ... % dx, dy, dz
                'const int, ',... % applyMotion
                'const float *, const float, const float, '... % angleRotXY, xCenterRotXY, yCenterRotXY
                'const float *, const float, const float, '... % angleRotXZ, xCenterRotXZ, zCenterRotXZ
                'const float *, const float, const float, '... % angleRotYZ, yCenterRotYZ, zCenterRotYZ
                'const float *, const float*, const float*' ]; % transX, transY, transZ

%% RF kernel -- kernel for RF parts
switch lower(simControl.odeMethod)
    case 'analytical'
        kernelFunction = 'diffusionExpRF';
        rfMode = 'Analytical';
    otherwise
        % Explicit 
        kernelFunction = 'diffusionExplicitODE';
        rfMode = 'Explicit';
end
% define precision
switch lower(simControl.precision)
    case 'double'
        kernelFunction = sprintf('%sId', kernelFunction);
        simPrecision = 'Double Precision';
    otherwise
        kernelFunction = sprintf('%sIf', kernelFunction);
        simPrecision = 'Single Precision';
end
% prepare the kernel
KRF = parallel.gpu.CUDAKernel(kernelPtx, kernelProto, kernelFunction);
KRF.ThreadBlockSize = numThreads;
KRF.GridSize = numBlocks;

%% GR kernel -- kernel for only Gradients (no Readouts, no SWC)
kernelFunction = 'diffusionExpGR';
% define precision
switch lower(simControl.precision)
    case 'double'
        kernelFunction = sprintf('%sId', kernelFunction);
    otherwise
        kernelFunction = sprintf('%sIf', kernelFunction);
end
% prepare the kernel
KGR = parallel.gpu.CUDAKernel(kernelPtx, kernelProto, kernelFunction);
KGR.ThreadBlockSize = numThreads;
KGR.GridSize = numBlocks;

%% RO kernel -- gradient kernel with potential readouts and/or SWC
kernelFunction = 'diffusionExpRO';
switch lower(simControl.precision)
    case 'double'
        kernelFunction = sprintf('%sId', kernelFunction);
    otherwise
        kernelFunction = sprintf('%sIf', kernelFunction);
end
% prepare the kernel
KRO = parallel.gpu.CUDAKernel(kernelPtx, kernelProto, kernelFunction);
KRO.ThreadBlockSize = numThreads;
KRO.GridSize = numBlocks;
            
%% prepare stats
stats.tRF = 0.0; stats.nRF = 0;
stats.tGR = 0.0; stats.nGR = 0;
stats.tRO = 0.0; stats.nRO = 0;
stats.tDF = 0.0; stats.nDF = 0;
stats.tKernel   = 0.0;
stats.gpuType   = currentGPU.Name;
stats.kernel    = kernelPtx;
stats.precision = simPrecision;
stats.rfMode    = rfMode;
stats.threads   = numThreads;
stats.blocks    = numBlocks;
stats.numSteps  = numSteps;
stats.numRxs    = numRxs;
stats.numVoxels = numIso;

if dbgControl.mode
    fprintf(fid, '\nStaring Diffusion simulation using Kernel File %s',kernelPtx); 
    fprintf(fid, '\n      0%%');
    progress = 0;
end            

%% Loop on the sequence parts and simulate
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
            accPhaseX = 2*pi*gamma*abs(accPhaseX) + max(abs(dPx));
            accPhaseY = 2*pi*gamma*abs(accPhaseY) + max(abs(dPy));
            accPhaseZ = 2*pi*gamma*abs(accPhaseZ) + max(abs(dPz));
            
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
            msg = sprintf( 'Sequence type %s not supported', partType );
            ME = MException(['error:',functionName], '%s', msg);
            throw(ME);
    end

    % part limits: move it to zero based for CUDA
    timeIni = uint32(sequence.partLimits(ss,1)-1); 
    timeEnd = uint32(sequence.partLimits(ss,2)-1);
    
    tkernel = tic();
    %% launch kernel
    [Mx, My, Mz, Sx, Sy, Sz, dPx, dPy, dPz ] = feval( K, ...
        Mx, My, Mz, Sx, Sy, Sz, dPx, dPy, dPz, ...
        timeDiff, rxSignal, swcSignal, ...
        gxSignal, gySignal, gzSignal, ...
        rfmSignal, rfpSignal, rffSignal, ...
        x, y, z, ...
        pd, r1, r2, cs, bi, ...
        xDiffusion, yDiffusion, zDiffusion,...
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
    
    % data collection
    stats.tKernel = stats.tKernel + tkernel;
    switch lower(sequence.partType{ss})
        case 'rf'
            stats.tRF = stats.tRF + tkernel;
            stats.nRF = stats.nRF + single(timeEnd) - single(timeIni) + 1;
        case 'gr'
            stats.tGR = stats.tGR + tkernel;
            stats.nGR = stats.nGR + single(timeEnd) - single(timeIni) + 1;
        case 'ro'
            stats.tRO = stats.tRO + tkernel;
            stats.nRO = stats.nRO + single(timeEnd) - single(timeIni) + 1;
        otherwise % default
            msg = sprintf( 'Sequence type %s not supported', partType );
            ME = MException(['error:',functionName], '%s', msg);
            throw(ME);
    end
    
    %% store performance data and update progress
    if  dbgControl.mode
        % back end progress bar
        if round(single(timeEnd)/single(numSteps)*100) > progress+5
            progress = round(single(timeEnd)/single(numSteps)*100);
            fprintf(fid, '...%d%%', progress );
        end
    end
    
end

%% collect data from GPU
solution.Mx = gather(Mx);
solution.My = gather(My);
solution.Mz = gather(Mz);
solution.dMpDx = gather(dPx);
solution.dMpDy = gather(dPy);
solution.dMpDz = gather(dPz);
solution.Sx = gather(Sx);
solution.Sy = gather(Sy);
solution.Sz = gather(Sz);

%% report
if  dbgControl.mode
    tTotal = toc(tTotal);
    fprintf(fid, '\n');
    fprintf(fid, '\n%s : done for sequence %s, %d parts, IRL time %.3fs',...
        functionName, seqType, sequence.numParts,sequence.time(end));
    numIso = double(numIso);
    fprintf(fid, '\n  Current GPU       %d (%s)', ...
        currentGPU.Index, currentGPU.Name);
    fprintf(fid, '\n  # GPU Blocks      %d', numBlocks);
    fprintf(fid, '\n  # GPU Threads     %d', numThreads);
    fprintf(fid, '\n  # Isochromas      %d', numIso);
    fprintf(fid, '\n  # Time Steps      %d', numSteps);
    fprintf(fid, '\n  # Readouts        %d', numRxs);
    fprintf(fid, '\n  Elapsed Time      %.3fs / %d steps (%.1fps/voxel/step - %s)',...
        tTotal, numSteps, 1e12*tTotal/numIso/numSteps, simPrecision);
    fprintf(fid, '\n  Kernel  Time      %.3fs / %d steps (%.1fps/voxel/step - %s)',...
        stats.tKernel, numSteps, 1e12*stats.tKernel/numIso/numSteps, simPrecision);
    fprintf(fid, '\n  RF part Time      %.3fs / %d steps (%.1fps/voxel/step - %s %s )',...
        stats.tRF, stats.nRF, 1e12*stats.tRF/numIso/max(stats.nRF,1), rfMode, simPrecision);
    fprintf(fid, '\n  GR part Time      %.3fs / %d steps (%.1fps/voxel/step - %s)',...
        stats.tGR, stats.nGR, 1e12*stats.tGR/numIso/max(stats.nGR,1), simPrecision);
    fprintf(fid, '\n  RO part Time      %.3fs / %d steps (%.1fps/voxel/step - %s)',...
        stats.tRO, stats.nRO, 1e12*stats.tRO/numIso/max(stats.nRO,1), simPrecision);
    fprintf(fid, '\n');
    if fid ~=1
        fclose(fid);
    end
end
