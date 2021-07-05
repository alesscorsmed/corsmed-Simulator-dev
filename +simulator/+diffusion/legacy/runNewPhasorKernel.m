function [Sx,Sy,Sz,Mm,Mp,Mz,dMpDx, dMpDy, dMpDz] = runNewPhasorKernel(M0m,M0p,M0z,seq,model,phaseAcc)


numThreads = 256;
numBlocks  = ceil(single(model.niso)/256);                                    
                                    
%% convert sequence to GPU data
time        = gpuArray(single(seq.time(:)));
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
nTimes      = uint32(seq.nt);
nRx         = int32(seq.nrx);
seqType     = seq.type;

%% convert model to GPU data
nIso            = uint32(model.niso);
b0              = single(model.b0);
mu              = single(model.mu);
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

dMpDx = gpuArray(single(0*M0p));
dMpDy = gpuArray(single(0*M0p));
dMpDz = gpuArray(single(M0p./z));

Sx = gpuArray(single(zeros(1,nRx*nCoils)));
Sy = gpuArray(single(zeros(1,nRx*nCoils)));
Sz = gpuArray(single(zeros(1,nRx*nCoils)));

%% run kernel


kernel_ptx = '/efs-mount-point-MATLAB/EduTool-Jorge/EduTool-CUDA/cuda_kernel_22_phasor.ptx';

switch seqType

    
        case 'RFAPC'
        fprintf(1, '\n Running APC RF Phasor Kernel\n\n');
        
        kernel_function = 'abmRF';
        cproto = ['float *, float *, float *, ',... % Mm, Mp, Mz
            'float *, float *, float *, ',... % Sx,Sy,Sz
            'float *, unsigned int *, unsigned int *, ',... % Tdiff, RX, SWC
            'float *, float *, float *, ',... % pulse Gx, Gy and Gz
            'float *, float *, float *, ',... % pulse RFmag, RFphase and RFfreq
            'float *, float *, float *, '... % x, y, z
            'float *, float *, ',... % Bi, PD
            'float *, unsigned int *, '... % propertyValues, tissueType
            'int, float *, float *, ',... % nCoils, Cx, Cy
            'unsigned int, unsigned int, unsigned int, ',... % nIso, tStart, tEnd
            'float, float, float, ',... % gamma, b0, mu
            'float, float, float' ]; % dx, dy, dz
        K = parallel.gpu.CUDAKernel(kernel_ptx, cproto, kernel_function);
        K.ThreadBlockSize = numThreads;
        K.GridSize = numBlocks;
        
        [Mm, Mp, Mz, Sx, Sy, Sz] = feval(K, Mm, Mp, Mz, Sx, Sy, Sz,...
            timeDiff, rxSignal, swcSignal, ...
            gxSignal, gySignal, gzSignal, ...
            rfmSignal, rfpSignal, rffSignal, ...
            x, y, z, bInhom, pDInhom, tissueValues, voxelTissue, ...
            nCoils, coilMapsX, coilMapsY,...
            nIso, 0, nTimes, gamma, b0, mu, ...
            voxLx, voxLy, voxLz );
        
%         [Mm, Mp, Mz] = diffusion.kernel.abmRF(Mm, Mp, Mz, ...
%             timeDiff, rxSignal, swcSignal, ...
%             gxSignal, gySignal, gzSignal, ...
%             rfmSignal, rfpSignal, rffSignal, ...
%             x, y, z, bInhom, pDInhom, tissueValues, voxelTissue,...
%             nIso, 0, nTimes, gamma, b0, mu);
        
% Mm = gather(Mm);
% Mp = gather(Mp);
% Mz = gather(Mz);        
% Mx = Mm.*cos(Mp);
% My = Mm.*sin(Mp);

% figure(1);
% subplot(5,1,1);  hold on;
% plot(z, Mx, 'LineWidth', 2 ); ylabel('Mx');
% subplot(5,1,2);  hold on;
% plot(z, My, 'LineWidth', 2 ); ylabel('My');
% subplot(5,1,3);  hold on;
% plot(z, Mz, 'LineWidth', 2 ); ylabel('Mz');
% subplot(5,1,4);  hold on;
% plot(z, abs(Mm), 'LineWidth', 2 ); ylabel('|Mxy|');
% subplot(5,1,5);  hold on;
% plot(z, Mp*180/pi, 'LineWidth', 2 ); ylabel('Phase(Mxy)');
%         pause();
%     
    case 'RFEXP'
        fprintf(1, '\n Running EXP RF Phasor Kernel\n\n');
        
        kernel_function = 'expRF';
        cproto = ['float *, float *, float *, ',... % Mm, Mp, Mz
            'float *, float *, float *, ',... % Sx,Sy,Sz
            'float *, unsigned int *, unsigned int *, ',... % Tdiff, RX, SWC
            'float *, float *, float *, ',... % pulse Gx, Gy and Gz
            'float *, float *, float *, ',... % pulse RFmag, RFphase and RFfreq
            'float *, float *, float *, '... % x, y, z
            'float *, float *, ',... % Bi, PD
            'float *, unsigned int *, '... % propertyValues, tissueType
            'int, float *, float *, ',... % nCoils, Cx, Cy
            'unsigned int, unsigned int, unsigned int, ',... % nIso, tStart, tEnd
            'float, float, float, ',... % gamma, b0, mu
            'float, float, float' ]; % dx, dy, dz
        K = parallel.gpu.CUDAKernel(kernel_ptx, cproto, kernel_function);
        K.ThreadBlockSize = numThreads;
        K.GridSize = numBlocks;
        
        [Mm, Mp, Mz, Sx, Sy, Sz] = feval(K, Mm, Mp, Mz, Sx, Sy, Sz,...
            timeDiff, rxSignal, swcSignal, ...
            gxSignal, gySignal, gzSignal, ...
            rfmSignal, rfpSignal, rffSignal, ...
            x, y, z, bInhom, pDInhom, tissueValues, voxelTissue, ...
            nCoils, coilMapsX, coilMapsY,...
            nIso, 0, nTimes, gamma, b0, mu, ...
            voxLx, voxLy, voxLz );
        

      case 'RFEXPD'
        fprintf(1, '\n Running EXP Deriv RF Phasor Kernel\n\n');
        
        kernel_function = 'expRFD';
        cproto = ['float *, float *, float *, ',... % Mm, Mp, Mz
            'float *, float *, float *, ',... % dMpdX,dMpdY,dMpdZ
            'float *, unsigned int *, unsigned int *, ',... % Tdiff, RX, SWC
            'float *, float *, float *, ',... % pulse Gx, Gy and Gz
            'float *, float *, float *, ',... % pulse RFmag, RFphase and RFfreq
            'float *, float *, float *, '... % x, y, z
            'float *, float *, ',... % Bi, PD
            'float *, unsigned int *, '... % propertyValues, tissueType
            'int, float *, float *, ',... % nCoils, Cx, Cy
            'unsigned int, unsigned int, unsigned int, ',... % nIso, tStart, tEnd
            'float, float, float, ',... % gamma, b0, mu
            'float, float, float, float' ]; % dx, dy, dz, derivRes
        K = parallel.gpu.CUDAKernel(kernel_ptx, cproto, kernel_function);
        K.ThreadBlockSize = numThreads;
        K.GridSize = numBlocks;
        
        derivResolution = single(0.5/phaseAcc/voxLz);
        
        [Mm, Mp, Mz, dMpDx, dMpDy, dMpDz] = feval(K, Mm, Mp, Mz, dMpDx, dMpDy, dMpDz,...
            timeDiff, rxSignal, swcSignal, ...
            gxSignal, gySignal, gzSignal, ...
            rfmSignal, rfpSignal, rffSignal, ...
            x, y, z, bInhom, pDInhom, tissueValues, voxelTissue, ...
            nCoils, coilMapsX, coilMapsY,...
            nIso, 0, nTimes, gamma, b0, mu, ...
            voxLx, voxLy, voxLz, derivResolution );
         
        
    case 'GR'
        
        fprintf(1, '\n Running GR Phasor Kernel\n\n');
        
        kernel_function = 'expGR';
        cproto = ['float *, float *, float *, ',... % Mm, Mp, Mz
                  'float *, float *, float *, ',... % Sx,Sy,Sz
                  'float *, unsigned int *, unsigned int *, ',... % Tdiff, RX, SWC
                  'float *, float *, float *, ',... % pulse Gx, Gy and Gz
                  'float *, float *, float *, ',... % pulse RFmag, RFphase and RFfreq                  
                  'float *, float *, float *, '... % x, y, z
                  'float *, float *, ',... % Bi, PD
                  'float *, unsigned int *, '... % propertyValues, tissueType
                  'int, float *, float *, ',... % nCoils, Cx, Cy
                  'unsigned int, unsigned int, unsigned int, ',... % nIso, tStart, tEnd
                  'float, float, float, ',... % gamma, b0, mu
                  'float, float, float' ]; % dx, dy, dz
        K = parallel.gpu.CUDAKernel(kernel_ptx, cproto, kernel_function);
        K.ThreadBlockSize = numThreads;
        K.GridSize = numBlocks;
        
        [Mm, Mp, Mz, Sx, Sy, Sz] = feval(K, Mm, Mp, Mz, Sx, Sy, Sz,...
            timeDiff, rxSignal, swcSignal, ...
            gxSignal, gySignal, gzSignal, ...
            rfmSignal, rfpSignal, rffSignal, ...
            x, y, z, bInhom, pDInhom, tissueValues, voxelTissue, ...
            nCoils, coilMapsX, coilMapsY,...
            nIso, 0, nTimes, gamma, b0, mu, ...
            voxLx, voxLy, voxLz );
        
    otherwise
        
        fprintf(1, '\n Running DEF Phasor Kernel\n\n');
        
        kernel_function = 'add3';
        cproto = ['float *, float *, float *, ',... % Mm, Mp, Mz
                  'float *, float *, float *, ',... % Sx,Sy,Sz
                  'float *, unsigned int *, unsigned int *, ',... % Tdiff, RX, SWC
                  'float *, float *, float *, ',... % pulse Gx, Gy and Gz
                  'float *, float *, float *, ',... % pulse RFmag, RFphase and RFfreq                  
                  'float *, float *, float *, '... % x, y, z
                  'float *, float *, ',... % Bi, PD
                  'float *, unsigned int *, '... % propertyValues, tissueType
                  'int, float *, float *, ',... % nCoils, Cx, Cy
                  'unsigned int, unsigned int, unsigned int, ',... % nIso, tStart, tEnd
                  'float, float, float, ',... % gamma, b0, mu
                  'float, float, float' ]; % dx, dy, dz
        K = parallel.gpu.CUDAKernel(kernel_ptx, cproto, kernel_function);
        K.ThreadBlockSize = numThreads;
        K.GridSize = numBlocks;
        
        [Mm, Mp, Mz, Sx, Sy, Sz] = feval(K, Mm, Mp, Mz, Sx, Sy, Sz,...
            timeDiff, rxSignal, swcSignal, ...
            gxSignal, gySignal, gzSignal, ...
            rfmSignal, rfpSignal, rffSignal, ...
            x, y, z, bInhom, pDInhom, tissueValues, voxelTissue, ...
            nCoils, coilMapsX, coilMapsY,...
            nIso, 0, nTimes, gamma, b0, mu, ...
            voxLx, voxLy, voxLz );
        
        
end


%% collect data from GPU
Mm = gather(Mm);
Mp = gather(Mp);
Mz = gather(Mz);
dMpDx = gather(dMpDx);
dMpDy = gather(dMpDy);
dMpDz = gather(dMpDz);
Sx = gather(Sx);
Sy = gather(Sy);
Sz = gather(Sz);
