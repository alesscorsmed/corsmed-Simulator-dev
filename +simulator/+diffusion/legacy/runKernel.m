function [Sx,Sy,Sz,Mm,Mp,Mz] = runDiffusionKernel(...
                   initialMagnitude,initialPhase,initialFA,seq,model)


numThreads = 256;
numBlocks  = ceil(single(model.niso)/256);                                    
                                    
%% convert sequence to GPU data
timeDiff    = gpuArray(single(seq.tdiff(:)));
gxSignal    = gpuArray(single(seq.gx(:)));
gySignal    = gpuArray(single(seq.gy(:)));
gzSignal    = gpuArray(single(seq.gz(:)));
rfmSignal   = gpuArray(single(seq.rfm(:))); 
rfpSignal   = gpuArray(single(seq.rfp(:)));
rffSignal   = gpuArray(single(seq.rff(:)));
rxSignal    = gpuArray(uint32(seq.rx(:)));
swcSignal   = gpuArray(uint32(seq.swc(:)));
gamma       = gpuArray(single(seq.gamma));
nTimes      = int32(seq.nt);
nRx         = int32(seq.nrx);

%% convert model to GPU data
nIso            = uint32(model.niso);
x               = gpuArray(single(model.x(:)));
y               = gpuArray(single(model.y(:)));
z               = gpuArray(single(model.z(:)));
bInhom          = gpuArray(single(model.bi(:)));
pDInhom         = gpuArray(single(model.pd(:)));
tissueValues    = gpuArray(single(model.tissueValues(:)));  
voxelTissue     = gpuArray(uint32(model.voxelTissue(:)));  
coilMapsX       = gpuArray(single(model.coilx(:).'));  
coilMapsY       = gpuArray(single(model.coily(:).'));  
nCoils          = gpuArray(int32(model.ncoil));  
voxLx           = gpuArray(single(model.resolution(1)));
voxLy           = gpuArray(single(model.resolution(2)));
voxLz           = gpuArray(single(model.resolution(3)));
Dx              = gpuArray(single(model.diffx(:)));
Dy              = gpuArray(single(model.diffy(:)));
Dz              = gpuArray(single(model.diffz(:)));

%% Allocate space for magnetization (initialize) and Signal solution
Mm = gpuArray(single(initialMagnitude*sin(initialFA*pi/180)*ones(1,nIso)));
Mp = gpuArray(single(initialPhase*ones(1,nIso)));
Mz = gpuArray(single(initialMagnitude*cos(initialFA*pi/180)*ones(1,nIso)));

Sx = gpuArray(single(zeros(1,nRx*nCoils)));
Sy = gpuArray(single(zeros(1,nRx*nCoils)));
Sz = gpuArray(single(zeros(1,nRx*nCoils)));

%% run kernel
fprintf(1, '\n Running Diffusion Kernel\n\n');

kernel_ptx = '/efs-mount-point-MATLAB/EduTool-Jorge/EduTool-CUDA/cuda_kernel_21_phasor_selfdwi.ptx';
kernel_function = 'add3';
cproto = 'float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, unsigned int *, unsigned int *, unsigned int *, float *, float, int, float *, float *, float *, int, unsigned int, float, float, float, float *, float *, float *';
K = parallel.gpu.CUDAKernel(kernel_ptx, cproto, kernel_function);
K.ThreadBlockSize = numThreads;
K.GridSize = numBlocks;

[Mm, Mp, Mz, Sx, Sy, Sz] = feval(K, Mm, Mp, Mz, Sx, Sy, Sz,...
    gxSignal, gySignal, gzSignal, rfmSignal, rfpSignal, timeDiff, ...
    x, y, z, bInhom, tissueValues, voxelTissue, ...
    rxSignal, swcSignal, rffSignal, gamma, ...
    nCoils, coilMapsX, coilMapsY, pDInhom, ...
    nTimes, nIso, voxLx, voxLy, voxLz,...
    Dx, Dy, Dz);

%% collect data from GPU
Mm = gather(Mm);
Mp = gather(Mp);
Mz = gather(Mz);
Sx = gather(Sx);
Sy = gather(Sy);
Sz = gather(Sz);
