function [iSpace] = sosFFT( kSpace, encodingData, expControl )
%
% RECONSTRUCTION.CARTESIANFFT.SOSFFT
%
%	Image with Sum-of-Squares FFT based reconstruction.
%
% INPUT
%
%
% OUTPUT
%
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'reconstructor.cartesianFFT.sosFFT';
if (nargin < 3)
    ME = MException('Reconstructor:wrongArgCount',...
        '%s : wrong argument count',functionName);
    throw(ME);
end

%% info for debugging
if expControl.debug.debugMode
    try % open file if possible, otherwise dump to stdout
        fid = fopen(expControl.debug.debugFile,'a');
    catch
        fid = 1;
    end
    % time it
    tTotal = tic();
    fprintf(fid, '\n%s : start', functionName);
end

%% extract encoding plan and FT operators
encodingPlan = encodingData.plan;
try
    iFTZ    = encodingData.operator.iFTZ;
    iFT     = encodingData.operator.iFT;    
catch
    iFTZ    = @(kSpace) fftshift(ifftn(ifftshift(kSpace),[],3));
    iFT     = @(kSpace) fftshift(ifftn(ifftshift(kSpace)));
end

%% get info about size
[numX, numY, numZ, numCH, numCI] = size(kSpace);

%% K-space padding
% factor
padX        = max(round(encodingPlan.xPadFactor), 1);
padY        = max(round(encodingPlan.yPadFactor), 1);
padZ        = max(round(encodingPlan.zPadFactor), 1);
% start positions of actual data
xPadStart   = 1 + floor((padX-1)*numX/2);
yPadStart   = 1 + floor((padY-1)*numY/2);
zPadStart   = 1 + floor((padZ-1)*numZ/2);
% incidence of padding
xPadIdx     = xPadStart:xPadStart+numX-1; 
yPadIdx     = yPadStart:yPadStart+numY-1;
zPadIdx     = zPadStart:zPadStart+numZ-1;
% allocate with padded size
numX = padX*numX;
numY = padY*numY;
numZ = padZ*numZ;
kSpacePad   = zeros(numX,numY,numZ);

%% generate the images
% allocate the image space
iSpace = zeros([numX,numY,numZ,numCH,numCI]);
% loop on Chanels (CH) and Contrast Images (CI) and apply IFT
for ic = 1:numCI
    for cc = 1:numCH
        kSpacePad(xPadIdx,yPadIdx,zPadIdx) = kSpace(:,:,:,cc,ic);
        if encodingPlan.is3D
            % 3D recon
            iSpace(:,:,:,cc,ic) = iFT(kSpacePad);
        else
            % series of 2D recon
            for zz = 1:numZ
                iSpace(:,:,zz,cc,ic) = iFT(kSpacePad(:,:,zz));
            end
        end
    end
end

%% Sum-of-Squares on the coils
if numCH>1
    iSpace = sqrt(sum(conj(iSpace).*iSpace,4));
end

%% crop the image if over sampling is applied
if encodingPlan.factorX > 1
    voxelStart = 1+floor(numX/(2*encodingPlan.factorX));
    voxelIndex = voxelStart:voxelStart+numX/encodingPlan.factorX-1;
    iSpace = iSpace(voxelIndex,:,:,:,:);
end
%% crop the image if fold over suppresion is applied
if encodingPlan.factorY > 1
    voxelStart = 1+floor(numY/(2*encodingPlan.factorY));
    voxelIndex = voxelStart:voxelStart+numY/encodingPlan.factorY-1;
    iSpace = iSpace(:,voxelIndex,:,:,:);
end

%% final message
if expControl.debug.debugMode
    % correct sizes
    [numX,numY,numZ,~,numCI] = size(iSpace);
    fprintf(fid, ...
        '\n%s : FFT-based recon done, elapsed time %.3fs',...
        functionName, toc(tTotal));
    fprintf(fid, '\n  I-Space size   %3d x %d x %d x %d',...
        numX, numY, numZ, numCI);
    fprintf(fid, '\n');
    if fid ~=1
        fclose(fid);
    end   
end
