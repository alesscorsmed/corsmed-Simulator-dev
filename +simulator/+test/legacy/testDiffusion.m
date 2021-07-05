%% Script to test diffusion kernel
%
% -------------------------------------------------------------------------
%%  CORSMED AB © 2020
%   Jorge Fernandez Villena - jorge.villena@corsmed.com
% _________________________________________________________________________

valuesTAU = [25, 75, 150]*1e-3;
valuesGStrength = [0.03, 0.10, 0.30];
valuesGxScale = [-1.0, -0.2, 0.0, 0.4, 1.0]; 
valuesGyScale = [-1.0, -0.05, 0.0, 0.1  1.0];  
valuesGzScale = [-1.0, -0.6, 0.0, 0.01, 1.0]; 
gSlewRate = 150e3;
valuesDx = [5.0, 1.0]*1e-9;
valuesDy = [6.0, 3.0]*1e-9;
valuesDz = [4.0, 2.0]*1e-9;


%% initial conditions
initialMagnitude = 1;
initialPhase = 0;
initialFA = 90;

%% domain
numberVoxels = 25;
fovX = 0.005;  fovY = 0.005;  fovZ = 0.005;
t1Value = 0.500;
t2Value = 0.250;

% generate a cube with problem size and dimensions
dx = fovX/numberVoxels;
dy = fovY/numberVoxels;
dz = fovZ/numberVoxels;
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

% loop on cases
testcase = 0;
Bcoeff = zeros(9000,1);
ADCthe = zeros(9000,1);
ADCsim = zeros(9000,1);
for gStrength = valuesGStrength
    for gxScale = valuesGxScale
        for gyScale = valuesGyScale
            for gzScale = valuesGzScale
                for TAU = valuesTAU
                    
                    % generate a PGSE waveform
                    fov = 0.100; nSamples = 128; tStep = 1e-6;
                    [tDiff,gSignal,AG,TG,~] = diffusion.sequence.generatePGSE(TAU,fov,nSamples,gStrength,gSlewRate,tStep);
                    
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
                    seq.rx    = uint32(zeros(size(gSignal)));   seq.rx(end) = 1; % unique readout at end
                    seq.swc   = uint32(zeros(size(gSignal)));
                    seq.nrx   = int32(nnz(seq.rx));
                    seq.gamma = single(42.577478518e6);
                    seq.type  = 'GR';

                    
                    % NO diffusion as reference
                    [Sxref,Syref,Szref,Mmref,Mpref,Mzref] = diffusion.kernel.runPhasorKernel(initialMagnitude,initialPhase,initialFA,seq,model);
                    Mxref = Mmref.*cos(Mpref);
                    Myref = Mmref.*sin(Mpref);
                    
                    for Dx = valuesDx
                        for Dy = valuesDy
                            for Dz = valuesDz
                                
                                % update model with diffusion values
                                model.diffx = single(Dx*ones(model.niso,1));
                                model.diffy = single(Dy*ones(model.niso,1));
                                model.diffz = single(Dz*ones(model.niso,1));
                                
                                % Run Diffusion
                                [Sx,Sy,Sz,Mm,Mp,Mz] = diffusion.kernel.runDiffusionKernel(initialMagnitude,initialPhase,initialFA,seq,model);
                                Mx = Mm.*cos(Mp);
                                My = Mm.*sin(Mp);
                                
                                % Apparent diffusion coefficient
                                ADC = Sx/Sxref;
                                
                                % Theoretical ADC
                                WGAMMA2 = (2*pi*seq.gamma)^2;
                                AGx = AG*gxScale;
                                Bxcoeff = WGAMMA2 * (AGx^2) * (TG^2) * (TAU-TG/3);
                                AGy = AG*gyScale;
                                Bycoeff = WGAMMA2 * (AGy^2) * (TG^2) * (TAU-TG/3);
                                AGz = AG*gzScale;
                                Bzcoeff = WGAMMA2 * (AGz^2) * (TG^2) * (TAU-TG/3);
                                ADCref = exp( -Bxcoeff*abs(Dx) - Bycoeff*abs(Dy) - Bzcoeff*abs(Dz) );

                                testcase = testcase +1
                                Bcoeff(testcase) = Bxcoeff + Bycoeff + Bzcoeff;
                                ADCthe(testcase) = ADCref;
                                ADCsim(testcase) = ADC;
                                
                            end
                        end
                    end
                end
            end
        end
    end
end

figure(1); hold on;
plot(Bcoeff, ADCthe, 'go', 'LineWidth', 2, 'MarkerSize',5 );
plot(Bcoeff, ADCsim, 'b+', 'LineWidth', 1,'MarkerSize',5 );
LogErr = log10(abs(ADCsim-ADCthe));
LogErr(LogErr < -6) = -6;
plot(Bcoeff, LogErr, 'r.', 'MarkerSize',10 );
xlabel('B coefficient'); ylabel('ADC'); grid on;
ylim([-6.0001,1.0001]);
legend('Theoretical', 'Simulated', 'Abs. Error (Logscale)')
title(sprintf('Diffusion test: %d cases',testcase));

