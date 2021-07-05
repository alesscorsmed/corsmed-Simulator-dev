function [Sx,Sy,Sz,Mm,Mp,Mz] = runPhasorKernel(M0m,M0p,M0z,seq,model)


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

%% Allocate space for magnetization (initialize) and Signal solution
Mm = gpuArray(single(M0m));
Mp = gpuArray(single(M0p));
Mz = gpuArray(single(M0z));

Sx = gpuArray(single(zeros(1,nRx*nCoils)));
Sy = gpuArray(single(zeros(1,nRx*nCoils)));
Sz = gpuArray(single(zeros(1,nRx*nCoils)));

%% run kernel
fprintf(1, '\n Running Phasor Kernel\n\n');

kernel_ptx = '/efs-mount-point-MATLAB/EduTool-Jorge/EduTool-CUDA/cuda_kernel_21_phasor.ptx';
kernel_function = 'add3';
cproto = 'float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, unsigned int *, unsigned int *, unsigned int *, float *, float, int, float *, float *, float *, int, unsigned int, float, float, float';
K = parallel.gpu.CUDAKernel(kernel_ptx, cproto, kernel_function);
K.ThreadBlockSize = numThreads;
K.GridSize = numBlocks;

[Mm, Mp, Mz, Sx, Sy, Sz] = feval(K, Mm, Mp, Mz, Sx, Sy, Sz,...
    gxSignal, gySignal, gzSignal, rfmSignal, rfpSignal, timeDiff, ...
    x, y, z, bInhom, tissueValues, voxelTissue, ...
    rxSignal, swcSignal, rffSignal, gamma, ...
    nCoils, coilMapsX, coilMapsY, pDInhom, ...
    nTimes, nIso, voxLx, voxLy, voxLz );

%% collect data from GPU
Mm = gather(Mm);
Mp = gather(Mp);
Mz = gather(Mz);
Sx = gather(Sx);
Sy = gather(Sy);
Sz = gather(Sz);
