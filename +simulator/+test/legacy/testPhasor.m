%% Script to test phasor kernels
%
% -------------------------------------------------------------------------
%%  CORSMED AB © 2020
%   Jorge Fernandez Villena - jorge.villena@corsmed.com
% _________________________________________________________________________

clear all;
clc;
numberVoxels = [50, 25, 10, 5, 1];

%% initial conditions
initialMagnitude = 1;
initialPhase = 0;
initialFA = 90;

%% domain
fovX = 0.005;  fovY = 0.005;  fovZ = 0.005;
t1Value = 0.500;
t2Value = 0.250;

%% sequence
TAU = 50e-3;
gxScale = 0.3; gyScale = 0.3; gzScale = 0.3;
tStep = 1e-6;

% generate a PGSE waveform
gStrength = 0.05; gSlewRate = 150;
fov = 0.100; nSamples = 128;
[tDiff,gSignal,AG,TG,TAU] = diffusion.sequence.generatePGSE(TAU,fov,nSamples,gStrength,gSlewRate,tStep);

% generate the sequence structure
seq.nt    = int32(length(tDiff));
seq.dt    = single(tStep);
seq.tdiff = single(tDiff);
seq.time  = single(cumsum(tDiff(:)));
seq.rfm   = single(zeros(size(gSignal)));
seq.rfp	  = single(zeros(size(gSignal)));
seq.rff   = single(zeros(size(gSignal)));
seq.gx    = single(gSignal*gxScale);
seq.gy	  = single(gSignal*gyScale);
seq.gz    = single(gSignal*gzScale);
seq.rx    = uint32(1:seq.nt);
seq.swc   = uint32(zeros(size(gSignal)));
seq.nrx   = int32(nnz(seq.rx));
seq.gamma = single(42.577478518e6);
seq.type  = 'GR';


figure(1); clf;
figure(2); clf;
plot(cumsum(tDiff),gSignal, 'LineWidth', 2);
xlabel('time(s)'); ylabel('Gradients Signal');

for nVox = numberVoxels
    % generate a cube with problem size and dimensions
    dx = fovX/nVox;
    dy = fovY/nVox;
    dz = fovZ/nVox;
    [r,x,y,z] = diffusion.domain.generateCube(fovX,fovY,fovZ,dx,dy,dz);
    [nx,ny,nz,~] = size(r);
    if nx > 1
        dx = x(2)-x(1);
    else
        dx = fovX;
    end
    if ny > 1
        dy = y(2)-y(1);
    else
        dy = fovY;
    end
    if nz > 1
        dz = z(2)-z(1);
    else
        dz = fovZ;
    end
    
    % shift positions off center
    r = r + fovX;
    model_spatial = reshape(r,nx*ny*nz,3);
    
    model.name = 'Homogen Cube';
    model.dim = [nx,ny,nz];
    model.niso = nx*ny*nz;
    model.resolution = [dx, dy, dz];
    model.x = single(model_spatial(:,1));
    model.y = single(model_spatial(:,2));
    model.z = single(model_spatial(:,3));
    
    model.tissueValues = zeros(1,6);
    model.tissueValues(1,1) = t1Value;
    model.tissueValues(1,2) = t2Value;
    model.voxelTissue = ones(model.niso,1);
    
    model.t1 = single(t1Value*ones(model.niso,1));
    model.t2 = single(t2Value*ones(model.niso,1));
    model.pd = single(ones(model.niso,1));
    
    model.bi = single(zeros(model.niso,1));
    model.cs = single(zeros(model.niso,1));
    
    model.coilx = single(ones(1,model.niso));
    model.coily = single(zeros(1,model.niso));
    nCoils = size(model.coilx,1);
    model.ncoil = nCoils;
    
    model.b0 = 1.5;
    model.mu = 1e-6;
    
    
    % call
    [Sx,Sy,Sz,Mx,My,Mz] = diffusion.kernel.runStandardKernel(initialMagnitude,initialPhase,initialFA,seq,model);
    figure(1); hold on;
    plot(cumsum(tDiff), Sx, 'LineWidth', 2 );
    
    if nVox == 1
        [Sx,Sy,Sz,Mm,Mp,Mz] = diffusion.kernel.runPhasorKernel(initialMagnitude,initialPhase,initialFA,seq,model);
        Mx = Mm.*cos(Mp);
        My = Mm.*sin(Mp);
        
        figure(1); hold on;
        plot(cumsum(tDiff), Sx, 'k--', 'LineWidth',2 );
    end
    
end


legend('50x50x50, 0.1mm', '25x25x25, 0.2mm', '10x10x10, 0.5mm', '5x5x5, 1.0mm', '1x1x1, 5.0mm', 'Phasor 1x1x1, 5.0mm');
xlabel('time (s)'); ylabel('Sx');



