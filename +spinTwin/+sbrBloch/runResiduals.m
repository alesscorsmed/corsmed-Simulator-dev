function [solution, stats] = runResiduals( solution, ...
    sequence, simModel, simControl, dbgControl )
%
% SPINTWIN.SBRBLOCH.RUNJACOBIANS
%
%	Runs SBR Phasor Bloch simulation, calling CUDA kernels:
%       Choice of Numerical/Analytical for RF
%       Analytical exponential integration for remainder
%       Allows for single/double precision
%       Phase based signal integration
%   Returns a solution structure (see below format)
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
%                    (Input / Output) Magnetizations at the center of the voxel
%                     solution.Mx  array with voxel's x component magnetization
%                     solution.My  array with voxel's y component magnetization
%                     solution.Mz  array with voxel's z component magnetization
%                    (Input / Output) Spatial derivatives of the Phase
%                     solution.dPx  d(Phase)/dx
%                     solution.dPy  d(Phase)/dy
%                     solution.dPz  d(Phase)/dz
%                    (Input / Output) Parameter derivatives: derivative w.r.t. R1 = 1/T1
%                     solution.dR1x    x component: d(Mx)/dR1
%                     solution.dR1y    y component: d(My)/dR1
%                     solution.dR1z    z component: d(Mz)/dR1
%                    (Input / Output) Parameter derivatives: derivative w.r.t. R2 = 1/T2
%                     solution.dR2x    x component: d(Mx)/dR2
%                     solution.dR2y    y component: d(My)/dR2
%                     solution.dR2z    z component: d(Mz)/dR2
%                    (Input / constant) Reference (acquired) signal: flatten array with order numRxCoils x numRxs
%                     solution.Srefx   x component of signal
%                     solution.Srefy   y component of signal
%                    (Output) Integrated signal: numRxCoils x numRxs
%                     solution.Sx   x component of signal
%                     solution.Sy   y component of signal
%                    (Output) Residual signal (Sref - S): numRxCoils x numRxs
%                     solution.Rx   x component
%                     solution.Ry   y component
%                    EMPTY: (Output) Gradients: numIsochromats  x numRxCoils
%                     solution.GPr  gradient for real part of Proton Density
%                     solution.GPi  gradient for imag part of Proton Density
%                     solution.GR1  gradient for R1 (1/T1)
%                     solution.GR2  gradient for R2 (1/T2)
%                    EMPTY: (Output) Jacobians: numIsochromats x numRxCoils x numRxs
%                     solution.JPDx  x jacobian for Proton Density, from d(Mx)/d(PD)
%                     solution.JPDy  y jacobian for Proton Density, from d(My)/d(PD)
%                     solution.JR1x  x jacobian for R1 (1/T1), from d(Mx)/d(R1)
%                     solution.JR1y  y jacobian for R2 (1/T2), from d(My)/d(R1)
%                     solution.JR2x  x jacobian for R2 (1/T1), from d(Mx)/d(R2)
%                     solution.JR2y  y jacobian for R2 (1/T2), from d(My)/d(R2)
%                     
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'spinTwin:sbrBloch:runJacobians';
% check args
if (nargin < 4) ...
        || isempty(solution) || isempty(simModel) ...
        || isempty(sequence) || isempty(simControl)
    ME = MException(['error:',functionName],...
        'Wrong arguments.');
    throw(ME);
end
% optional arguments
if (nargin < 5) || isempty(dbgControl)
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
    gxSignal    = gpuArray(single(sequence.gxSignal(:)));
    gySignal    = gpuArray(single(sequence.gySignal(:)));
    gzSignal    = gpuArray(single(sequence.gzSignal(:)));
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
    %
    x       = gpuArray(single(simModel.x(indexes)));
    y       = gpuArray(single(simModel.y(indexes)));
    z       = gpuArray(single(simModel.z(indexes)));
    % voxel properties
    pd      = gpuArray(single(simModel.pr(indexes)));
    pi      = gpuArray(single(simModel.pi(indexes)));
    r1      = gpuArray(single(simModel.r1(indexes)));
    r2      = gpuArray(single(simModel.r2(indexes)));
    cs      = gpuArray(single(simModel.cs(indexes)));
    bi      = gpuArray(single(simModel.bi(indexes)));
    % coil sens
    coilMapsX	= gpuArray(single(simModel.rxCoilMapsX(indexes,:)));
    coilMapsY	= gpuArray(single(simModel.rxCoilMapsY(indexes,:)));
    nCoils      = gpuArray(int32(simModel.numRxCoils));
    % voxel volume
    voxLx       = gpuArray(single(simModel.resolution(1)));
    voxLy       = gpuArray(single(simModel.resolution(2)));
    voxLz       = gpuArray(single(simModel.resolution(3)));
    
catch
    msg = sprintf( 'Failed to transfer tissue model to GPU' );
    ME = MException(['error:',functionName], '%s', msg);
    throw(ME);
end

%% convert motion model to GPU data
try
    
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
    Sz = gpuArray(single(0*solution.Sz));
    
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
        kernelFunction = 'phasorExpRF';
        rfMode = 'Analytical';
    otherwise
        % Explicit 
        kernelFunction = 'phasorExplicitODE';
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
kernelFunction = 'phasorExpGR';
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
kernelFunction = 'phasorExpRO';
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
    fprintf(fid, '\nStaring Phasor simulation using Kernel File %s',kernelPtx);
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

solution.dPx = gather(dPx);
solution.dPy = gather(dPy);
solution.dPz = gather(dPz);

solution.Sx = gather(Sx);
solution.Sy = gather(Sy);

solution.Rx = solution.Srefx - solution.Sx;
solution.Ry = solution.Srefy - solution.Sy;

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
