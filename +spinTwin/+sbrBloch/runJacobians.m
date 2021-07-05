function [solution, stats] = runJacobians( solution, ...
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
%                    (Output) Jacobians: numIsochromats x numRxCoils x numRxs
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
    pr      = gpuArray(single(simModel.pr(indexes)));
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

%% Allocate space for magnetization (initialize) and Signal solution
% get the max number of samples in any part
maxLocalRx = 0;
for ss = 1:sequence.numParts
    % part limits
    timeIni = sequence.partLimits(ss,1); 
    timeEnd = sequence.partLimits(ss,2);
    maxLocalRx = max(maxLocalRx,nnz(sequence.rxSignal(timeIni:timeEnd)));
end
try
    Mx = gpuArray(single(solution.Mx));
    My = gpuArray(single(solution.My));
    Mz = gpuArray(single(solution.Mz));
    
    dR1x = gpuArray(single(solution.dR1x));
    dR1y = gpuArray(single(solution.dR1y));
    dR1z = gpuArray(single(solution.dR1z));
    
    dR2x = gpuArray(single(solution.dR2x));
    dR2y = gpuArray(single(solution.dR2y));
    dR2z = gpuArray(single(solution.dR2z));
    
    dPx = gpuArray(single(solution.dPx));
    dPy = gpuArray(single(solution.dPy));
    dPz = gpuArray(single(solution.dPz));
    
    % generated signal
    Sx = gpuArray(single(reshape(solution.Sx,nCoils,numRxs)));
    Sy = gpuArray(single(reshape(solution.Sy,nCoils,numRxs)));
    
    % reference signal
    Srefx = gpuArray(single(reshape(solution.Srefx,nCoils,numRxs)));
    Srefy = gpuArray(single(reshape(solution.Srefy,nCoils,numRxs)));
    
    % residual function
    Rx = 0*Srefx;
    Ry = 0*Srefy;
    
    % allocate for local Jacobians
    JPDx = gpuArray(single(zeros(numIso,nCoils,maxLocalRx)));
    JPDy = gpuArray(single(zeros(numIso,nCoils,maxLocalRx)));
    
    JR1x = gpuArray(single(zeros(numIso,nCoils,maxLocalRx)));
    JR1y = gpuArray(single(zeros(numIso,nCoils,maxLocalRx)));
    
    JR2x = gpuArray(single(zeros(numIso,nCoils,maxLocalRx)));
    JR2y = gpuArray(single(zeros(numIso,nCoils,maxLocalRx)));
    
    % allocate for Jacobians if needed
    try 
        solution.JPDx = reshape(solution.JPDx,numIso,nCoils,numRxs);
    catch
        solution.JPDx = zeros(numIso,nCoils,numRxs);
    end
    try 
        solution.JPDy = reshape(solution.JPDy,numIso,nCoils,numRxs);
    catch
        solution.JPDy = zeros(numIso,nCoils,numRxs);
    end
    try 
        solution.JR1x = reshape(solution.JR1x,numIso,nCoils,numRxs);
    catch
        solution.JR1x = zeros(numIso,nCoils,numRxs);
    end
    try 
        solution.JR1y = reshape(solution.JR1y,numIso,nCoils,numRxs);
    catch
        solution.JR1y = zeros(numIso,nCoils,numRxs);
    end
    try 
        solution.JR2x = reshape(solution.JR2x,numIso,nCoils,numRxs);
    catch
        solution.JR2x = zeros(numIso,nCoils,numRxs);
    end
    try 
        solution.JR2y = reshape(solution.JR2y,numIso,nCoils,numRxs);
    catch
        solution.JR2y = zeros(numIso,nCoils,numRxs);
    end
   
catch
    msg = sprintf( 'Failed to transfer solution to GPU' );
    ME = MException(['error:',functionName], '%s', msg);
    throw(ME);
end

%% prepare to run the kernels
kernelPtx  = simControl.ptxSBR;
numThreads = simControl.threads;
numBlocks  = ceil(single(numIso)/numThreads);
% kernel proto is common to all cases
kernelProto = [ 'float *, float *, float *,' ... % Mx, My, Mz
                'float *, float *, float *,' ... % dR1x,dR1y,dR1z
                'float *, float *, float *,' ... % dR2z,dR2y,dR2z
                'float *, float *, float *,' ... % dPx,dPy,dPz
                'float *, float *,' ... % Sx, Sy
                'float *, float *,' ... % JPDx,JPDy
                'float *, float *,' ... % JR1x,JR1y
                'float *, float *,' ... % JR2z,JR2y
                'const float *,' ...    % timeDiff
                'const unsigned int *, const unsigned int *,' ... % RX, SWC
                'const float *, const float *, const float *,' ... % pulse Gx, Gy and Gz
                'const float *, const float *, const float *,' ... % pulse RFmag, RFphase and RFfreq
                'const float *, const float *, const float *,' ... % x, y, z
                'const float *, const float *,' ... % pr, pi
                'const float *, const float *,' ... % r1, r2
                'const float *, const float *,' ... % cs, bi
                'const int, const float *, const float *,' ... % nCoils, Cx, Cy
                'const unsigned int, const unsigned int, const unsigned int,' ... % nIso, tStart, tEnd
                'const float, const float, const float,' ... % gamma, b0, mu
                'const float, const float, const float' ]; % dx, dy, dz
            
%% RF kernel -- kernel for RF parts
% Explicit 
kernelFunction = 'sbrExplicitODE';
rfMode = 'Explicit';
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
kernelFunction = 'sbrExpGR';
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
kernelFunction = 'sbrExpRO';
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
stats.tGrads    = 0.0;
stats.tJacs     = 0.0;
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
    fprintf(fid, '\nStaring Gradient Bloch simulation using Kernel File %s',kernelPtx);
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
    localRxSignal   = sequence.rxSignal(timeIni:timeEnd);
    idxRx           = localRxSignal(:) > 0;
    idxRx           = localRxSignal(idxRx); % get indexes of the samples in Signal
    localNumRx      = length(idxRx(:));
    switch lower(partType)
        case 'rf'
            % use standard RF kernel
            K = KRF;
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
    [Mx, My, Mz, dR1x, dR1y, dR1z, dR2x, dR2y, dR2z, dPx, dPy, dPz, ...
        Sx, Sy, JPDx, JPDy, JR1x, JR1y, JR2x, JR2y ]...
        = feval( K, ...
        Mx, My, Mz, dR1x, dR1y, dR1z, dR2x, dR2y, dR2z, dPx, dPy, dPz, ...
        Sx, Sy, JPDx, JPDy, JR1x, JR1y, JR2x, JR2y, ...
        timeDiff, rxSignal, swcSignal, ...
        gxSignal, gySignal, gzSignal, ...
        rfmSignal, rfpSignal, rffSignal, ...
        x, y, z, pr, pi, r1, r2, cs, bi, ...
        nCoils, coilMapsX, coilMapsY,...
        numIso, timeIni, timeEnd, gamma, b0, mu, ...
        voxLx, voxLy, voxLz );
    
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
    
    %% store the jacobians
    tJacs = tic();
    % residual for the indexes of the local readouts
    Rx(:,idxRx) = Srefx(:,idxRx) - Sx(:,idxRx);
    Ry(:,idxRx) = Srefy(:,idxRx) - Sy(:,idxRx);

    % jacobians
    solution.JPDx(:,:,idxRx) = gather(JPDx(:,:,1:localNumRx));
    solution.JPDy(:,:,idxRx) = gather(JPDy(:,:,1:localNumRx));
    solution.JR1x(:,:,idxRx) = gather(JR1x(:,:,1:localNumRx));
    solution.JR1y(:,:,idxRx) = gather(JR1y(:,:,1:localNumRx));
    solution.JR2x(:,:,idxRx) = gather(JR2x(:,:,1:localNumRx));
    solution.JR2y(:,:,idxRx) = gather(JR2y(:,:,1:localNumRx));
	tJacs = toc(tJacs);
    % data collection
    stats.tJacs = stats.tJacs + tJacs;
    
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

solution.dR1x = gather(dR1x);
solution.dR1y = gather(dR1y);
solution.dR1z = gather(dR1z);

solution.dR2x = gather(dR2x);
solution.dR2y = gather(dR2y);
solution.dR2z = gather(dR2z);

solution.Sx = gather(Sx);
solution.Sy = gather(Sy);

solution.Rx = gather(Rx);
solution.Ry = gather(Ry);

solution.GPr = [];
solution.GPi = [];
solution.GR1 = [];
solution.GR2 = [];


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
    fprintf(fid, '\n  Jacobian Assembly %.3fs', stats.tJacs);
    fprintf(fid, '\n');
    if fid ~=1
        fclose(fid);
    end
end
