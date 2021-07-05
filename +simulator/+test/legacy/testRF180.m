%% Script to test FA
%
% -------------------------------------------------------------------------
%%  CORSMED AB © 2020
%   Jorge Fernandez Villena - jorge.villena@corsmed.com
% _________________________________________________________________________

%% define RF conditions
gamma           = 42.577478518e6;
b0              = 1.0;
mu              = 1;
f0              = gamma*b0;
dt              = 1e-5;
sliceThickness  = 6e-3;
receiverBW      = 200e3;
rfDuration      = 3e-3;
rfFlipAngle     = 180;
rfPhase         = pi/2;
rfCycles        = 1;
maxGStrenght    = 40e-3;
rwWeight        = 0*1.022;
doSliceSelect   = 0;

% initial magnetization state
initialMagnitude    = b0*mu;
initialPhase        = 0.0;
initialFA           = 90;
initialPhaseAcc     = 2*pi*gamma*maxGStrenght*5e-3;


%% domain
fovX = 0.001;
fovY = 0.001;  
fovZ = 2*sliceThickness;
dx   = 1e-3;
dy   = 1e-3;
dz   = 1e-4;
t1Value = 5.000;
t2Value = 3.000;

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
r = r;
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

model.b0 = b0;
model.mu = mu;


%% sequence
[RFblock] = diffusion.sequence.generateRFBlock(gamma,dt,sliceThickness,...
                                                rfFlipAngle,rfPhase,...
                                                rfDuration,rfCycles,...
                                                maxGStrenght,...
                                                doSliceSelect,rwWeight);

% generate the sequence structure
nt        = size(RFblock,2);
seq.nt    = int32(nt);
seq.dt    = single(dt);
seq.time  = single(dt*(1:nt));
seq.tdiff = single(dt*ones(1,nt));
seq.rfm   = single(RFblock(1,:));
seq.rfp	  = single(RFblock(2,:));
seq.rff   = single(RFblock(3,:));
seq.gx    = single(RFblock(4,:));
seq.gy	  = single(RFblock(5,:));
seq.gz    = single(RFblock(6,:));
seq.rx    = uint32(zeros(1,nt)); seq.rx(end) = 1; % unique readout at end
seq.swc   = uint32(zeros(1,nt));
seq.nrx   = int32(nnz(seq.rx));
seq.gamma = single(gamma);
seq.type  = 'RFEXPD';

%% initial magnetizations
M0m = initialMagnitude*sin(initialFA*pi/180)*ones(model.niso,1);
M0p = initialPhase*ones(model.niso,1) + initialPhaseAcc*reshape(z,model.niso,1);
M0z = initialMagnitude*cos(initialFA*pi/180)*ones(model.niso,1);
M0x = M0m.*cos(M0p(:));
M0y = M0m.*sin(M0p(:));

phaseAcc = initialPhaseAcc;

%% call default kernel
%[Sx,Sy,Sz,Mx,My,Mz] = diffusion.kernel.runStandardKernel(M0x,M0y,M0z,seq,model);

%% call new phasor kernel
[Sx,Sy,Sz,Mm,Mp,Mz,dMpDx,dMpDy,dMpDz] = diffusion.kernel.runNewPhasorKernel(M0m,M0p,M0z,seq,model,phaseAcc);
Mx = Mm.*cos(Mp);
My = Mm.*sin(Mp);
Mpnew = Mp;

figure(1);
subplot(5,1,1);  hold on;
plot(z, Mx, 'LineWidth', 2 ); ylabel('Mx');
subplot(5,1,2);  hold on;
plot(z, My, 'LineWidth', 2 ); ylabel('My');
subplot(5,1,3);  hold on;
plot(z, Mz, 'LineWidth', 2 ); ylabel('Mz');
subplot(5,1,4);  hold on;
plot(z, abs(Mm), 'LineWidth', 2 ); ylabel('|Mxy|');
subplot(5,1,5);  hold on;
plot(z, Mpnew*180/pi, 'LineWidth', 2 ); ylabel('Phase(Mxy)');


%% call phasor kernel
[Sx,Sy,Sz,Mm,Mp,Mz] = diffusion.kernel.runPhasorKernel(M0m,M0p,M0z,seq,model);
Mx = Mm.*cos(Mp);
My = Mm.*sin(Mp);

figure(1);
subplot(5,1,1);  hold on;
plot(z, Mx, 'LineWidth', 2 ); ylabel('Mx');
legend('Default', 'Phasor');
subplot(5,1,2);  hold on;
plot(z, My, 'LineWidth', 2 ); ylabel('My');
legend('Default', 'Phasor');
subplot(5,1,3);  hold on;
plot(z, Mz, 'LineWidth', 2 ); ylabel('Mz');
legend('Default', 'Phasor');
subplot(5,1,4);  hold on;
plot(z, Mm, 'LineWidth', 2 ); ylabel('|Mxy|');
legend('Default', 'Phasor');
subplot(5,1,5);  hold on;
plot(z, Mp*180/pi, 'LineWidth', 2 ); ylabel('Phase(Mxy)');
legend('Default', 'Phasor');
xlabel('z position (m)');


figure(2); hold on;
plot(z, Mpnew(:), 'LineWidth', 2 ); ylabel('Phase(Mxy)');
plot(z, Mp(:), 'LineWidth', 2 ); ylabel('Phase(Mxy)');
plot(z, M0p(:), ':', 'LineWidth', 2 ); ylabel('Phase(Mxy)');
legend('New', 'Phasor', 'Iinitial');
xlabel('z position (m)');


figure(3); hold on;
rw = initialPhaseAcc*ones(size(z));
%plot(z, [0; diff(unwrap(Mpnew))/dz], 'LineWidth', 2 ); ylabel('dPhase(Mxy)');
plot(z, dMpDz, 'LineWidth', 2 ); ylabel('dPhase(Mxy)');
plot(z, -rw, ':', 'LineWidth', 2 ); ylabel('dPhase(Mxy)');
%legend('Default', 'Phasor', 'RW');
xlabel('z position (m)')


