function [Sx,Sy,Sz,Mx,My,Mz] = runStandardKernel(M0x,M0y,M0z,seq,model)


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
% if sizes are smaller than required, fill with zeros
if numThreads*numBlocks > model.niso
    model.x(numThreads*numBlocks) = 0.0;
    model.y(numThreads*numBlocks) = 0.0;
    model.z(numThreads*numBlocks) = 0.0;
    model.pd(numThreads*numBlocks) = 0.0;
    model.r1(numThreads*numBlocks) = 0.0;
    model.r2(numThreads*numBlocks) = 0.0;
    model.cs(numThreads*numBlocks) = 0.0;
    model.bi(numThreads*numBlocks) = 0.0;
    model.voxelTissue(numThreads*numBlocks) = 0;
    model.coilx(:,numThreads*numBlocks) = 0.0;
    model.coily(:,numThreads*numBlocks) = 0.0; 
end
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
Mx = zeros(1,numThreads*numBlocks);
My = zeros(1,numThreads*numBlocks);
Mz = zeros(1,numThreads*numBlocks);

Mx(1:nIso) = M0x(:);
My(1:nIso) = M0y(:);
Mz(1:nIso) = M0z(:);

Mx = gpuArray(single(Mx));
My = gpuArray(single(My));
Mz = gpuArray(single(Mz));

Sx = gpuArray(single(zeros(1,nRx*nCoils)));
Sy = gpuArray(single(zeros(1,nRx*nCoils)));
Sz = gpuArray(single(zeros(1,nRx*nCoils)));

%% run kernel
fprintf(1, '\n Running XYZ Kernel\n\n');

kernel_ptx = '/efs-mount-point-MATLAB/EduTool-Jorge/EduTool-CUDA/cuda_kernel_20_PDinhom.ptx';
kernel_function = 'add3';
cproto = 'float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, unsigned int *, unsigned int *, unsigned int *, float *, float, int, float *, float *, float *, int';
K = parallel.gpu.CUDAKernel(kernel_ptx, cproto, kernel_function);
K.ThreadBlockSize = numThreads;
K.GridSize = numBlocks;

[Mx, My, Mz, Sx, Sy, Sz] = feval(K, Mx, My, Mz, Sx, Sy, Sz,...
    gxSignal, gySignal, gzSignal, rfmSignal, rfpSignal, timeDiff, ...
    x, y, z, bInhom, tissueValues, voxelTissue, ...
    rxSignal, swcSignal, rffSignal, gamma, ...
    nCoils, coilMapsX, coilMapsY, pDInhom, nTimes );

%% collect data from GPU
Mx = gather(Mx(1:nIso));
My = gather(My(1:nIso));
Mz = gather(Mz(1:nIso));
Sx = gather(voxLx*voxLy*voxLz*Sx);
Sy = gather(voxLx*voxLy*voxLz*Sy);
Sz = gather(voxLx*voxLy*voxLz*Sz);
